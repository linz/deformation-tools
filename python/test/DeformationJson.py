#!/usr/bin/python3
#
# DeformationModelJson: Script to load and calculate a deformation model formatted
# with a JSON master file and GeoTIFF grid files.

from datetime import datetime
import os
import os.path
import json
import re
import math
import hashlib
import numpy as np
from DictObject import DictObject, Field, FormatError
from DeformationGrid import DeformationGridGeoTIFF

DisplacementFields = {
    "none": [],
    "horizontal": ["east_offset", "north_offset"],
    "vertical": ["vertical_offset"],
    "3d": ["east_offset", "north_offset", "vertical_offset"],
}

UncertaintyFields = {
    "none": [],
    "horizontal": ["horizontal_uncertainty"],
    "vertical": ["vertical_uncertainty"],
    "3d": ["horizontal_uncertainty", "vertical_uncertainty"],
}

DegreeFields = ["east_offset", "north_offset"]
BilinearMethod = "bilinear"
BilinearGeocentricMethod = "geocentric_bilinear"
BilinearGeocentricTypes = ["horizontal", "3d"]

E_OFFSET_ORDINATE = 0
N_OFFSET_ORDINATE = 1
U_OFFSET_ORDINATE = 2
H_UNCERTAINTY_ORDINATE = 3
V_UNCERTAINTY_ORDINATE = 4


class Context:
    def __init__(self, sourcefile, check=False):
        self.sourcefile = sourcefile
        self.basedir = os.path.dirname(sourcefile)
        self.check = check
        self.errors = []

    def raiseError(self, error):
        if self.check:
            self.errors.append(error)
        else:
            raise ValueError(error)


class TimeFunction:

    types = dict()

    @staticmethod
    def factory(value, context=None):
        funcdef = DictObject([Field("type", TimeFunction.types), Field("parameters", dict)], value, context=context)
        return TimeFunction.types[funcdef.types](funcdef.parameters, context=context)

    @staticmethod
    def parseTime(value):
        return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ")

    @staticmethod
    def decimalYear(value):
        if type(value) == str:
            value = TimeFunction.parseTime(value)
        if type(value) == datetime:
            r0 = datetime(value.year, 1, 1)
            r1 = datetime(value.year + 1, 1, 1)
            value = value.year + (value - r0).total_seconds() / (r1 - r0).total_seconds()
        if type(value) != float:
            raise ValueError("Invalid value {0} for date/time".format(value))
        return value

    def valueAt(self, epoch):
        raise RuntimeError("TimeFunction.valueAt must be overridden")


class VelocityFunction(DictObject, TimeFunction):

    TimeFunction.type["velocity"] = VelocityFunction
    fields = [Field("reference_epoch", TimeFunction.decimalYear)]

    def __init__(self, value, context=None):
        DictObject.init(self, self.fields, value, context=context)

    def valueAt(self, epoch):
        return epoch - self.reference_epoch


class PiecewiseFunctionPoint(DictObject):

    fields = [Field("epoch", TimeFunction.decimalYear), Field("scale_factor", float)]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)


class PiecewiseFunction(DictObject):

    TimeFunction.types["piecewise"] = PiecewiseFunction

    ZERO = 0
    CONSTANT = 1
    LINEAR = 2
    Extrapolation = {
        "zero": PiecewiseFunction.ZERO,
        "constant": PiecewiseFunction.CONSTANT,
        "linear": PiecewiseFunction.LINEAR,
    }

    fields = [
        Field("model", [PiecewiseFunctionPoint]),
        Field("before_first", lambda v: PiecewiseFunction.Extrapolation[v]),
        Field("after_last", lambda v: PiecewiseFunction.Extrapolation[v]),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)
        if len(self.model) < 1:
            raise ValueError("Piecewise function must contain at least one function point")
        for p0, p1 in zip(self.model[:-1], self.model[1:]):
            if p1.epoch < p0.epoch:
                raise ValueError("Piecewise function epoch {0} cannot be less than {1}".format(p1.epoch, p0.epoch))
        model = self.model
        if self.before_first == self.LINEAR and (len(model) < 2 or model[0].epoch == model[1].epoch):
            raise ValueError("Cannot use linear extrapolation on first two points of piecewise function")
        if self.after_last == self.LINEAR and (len(model) < 2 or model[-2].epoch == model[-1].epoch):
            raise ValueError("Cannot use linear extrapolation on last two points of piecewise function")

    def linearExtrapolation(self, r0, r1, epoch):
        if r0.epoch == r1.epoch:
            return r1.value
        return ((epoch - r0.epoch) * r1.value + (r1.epoch - epoch) * r0.value) / (r1.epoch - r0.epoch)

    def valueAt(self, epoch):
        if epoch <= self.model[0].epoch:
            if self.before_first == self.ZERO:
                return 0.0
            elif self.before_first == self.CONSTANT:
                return self.model[0].value
            else:
                return self.linearExtrapolation(self.refpoint[0], self.refpoint[1], epoch)
        elif epoch > self.model[1].epoch:
            if self.after_last == self.ZERO:
                return 0.0
            elif self.after_last == self.CONSTANT:
                return self.model[-1].value
            else:
                return self.linearExtrapolation(self.refpoint[-2], self.refpoint[-1], epoch)
        for r0, r1 in zip(self.model[:-1], self.model[1:]):
            if epoch <= r1.epoch:
                return self.linearExtrapolation(r0, r1, epoch)
        raise RuntimeError("Failed to find epoch {0} in piecewise time function".format(epoch))


class StepFunction(DictObject, PiecewiseFunction):

    TimeFunction.types["step"] = StepFunction

    fields = [Field("step_epoch", TimeFunction.decimalYear)]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)
        PiecewiseFunction.__init__(
            self,
            {"before_first": "zero", "after_last": "constant", "model": [{"epoch": self.step_epoch, "scale_factor": 1.0}]},
        )


class ReverseStepFunction(DictObject, PiecewiseFunction):

    TimeFunction.types["reverse_step"] = ReverseStepFunction

    fields = [Field("step_epoch", TimeFunction.decimalYear)]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)
        PiecewiseFunction.__init__(
            self,
            {"before_first": "constant", "after_last": "zero", "model": [{"epoch": self.step_epoch, "scale_factor": 1.0}]},
            context=context,
        )


class ExponentialFunction(DictObject, TimeFunction):

    TimeFunction.types["exponential"] = ExponentialFunction

    fields = [
        Field("reference_epoch", TimeFunction.decimalYear),
        Field("end_epoch", TimeFunction.decimalYear, optional=True),
        Field("relaxation_constant", float, min=0.0),
        Field("before_scale_factor", float),
        Field("initial_scale_factor", float),
        Field("final_scale_factor", float),
    ]


def __init__(self, value, context=None):
    DictObject.__init__(self, self.fields, value, context=context)


def valueAt(self, epoch):
    if self.end_epoch is not None and epoch > self.end_epoch:
        epoch = self.end_epoch
    if epoch < self.start_epoch:
        return self.before_scale_factor
    return self.initial_scale_factor + (
        (self.final_scale_factor - self.initial_scale_factor)
        * (1.0 - math.exp((self.ref_epoch - epoch) / self.relaxation_constant))
    )


class TimeExtent(DictObject):

    fields = [Field("first", TimeFunction.decimalYear), Field("last", TimeFunction.decimalYear)]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)

    def contains(self, epoch):
        return epoch >= self.first and epoch < self.last


class Extent:
    def __init__(self, value, context=None):
        if "type" not in value or value["type"] != "bbox":
            raise ValueError("Invalid extent type {0} - must be bbox".format(value.get(type, "undefined")))
        if type(value.get("parameters")) != dict or "bbox" not in type(value["parameters"]):
            raise ValueError("Invalid extent - requires bbox parameter")
        try:
            extents = value["parameters"]["bbox"]
            xmin, ymin, xmax, ymax = (float(x) for x in extents)
            if xmax < xmin or ymax < ymin:
                raise ValueError("Error in bounding box - max value less than min value")
        except:
            raise ValueError("Invalid bounding box - must be [xmin,ymin,xmax,ymax]")

    def contains(self, x, y):
        if x < xmin or x >= xmax:
            return False
        if y < ymin or y >= ymax:
            return False


class SourceFile:
    def __init__(self, filename, md5=None, context=None):
        self.filename = filename
        basedir = context.basedir if context is not None else ""
        filepath = os.path.join(basedir, filename)
        if not os.path.isfile(filepath):
            raise ValueError("File {0} is missing".format(filepath))
        if context is not None and context.check and md5 is not None:
            if self.calcMd5() != md5:
                raise ValueError("File {0} md5 checksum is not correct")

    def calcMd5(self):
        md5 = hashlib.md5()
        with open(self.filepath, "rb") as gf:
            while True:
                buffer = gf.read(2048)
                if len(buffer) == 0:
                    break
                md5.update(buffer)
        return md5.hexdigest()


class SpatialModel(DictObject):

    displacementTypes = {t: t for t in DisplacementFields}
    uncertaintyTypes = {t: t for t in UncertaintyFields}

    fields = [
        Field("type", str, re="^GeoTIFF$"),
        Field("interpolation_method", str, "^(bilinear|geocentric_bilinear)$"),
        Field("filename", str),
        Field("md5_checksum", str),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)
        self.displacement_type = None
        self.uncertainty_type = None
        sourcefile = SourceFile(self.filename, self.md5_checksum, context)
        self.grid = DeformationGridGeoTIFF(sourcefile.filepath)

    def check(self, displacement_type, uncertainty_type, context):
        self.grid.load()
        grid_disptype = self.grid.displacement_type()
        grid_unctype = self.grid.uncertainty_type()

        if grid_disptype != displacement_type:
            raise ValueError("Grid displacement type {0} doesn't match expected {1}".format(grid_disptype, displacement_type))
        if grid_unctype != uncertainty_type:
            raise ValueError("Grid uncertainty type {0} doesn't match expected {1}".format(grid_unctype, uncertainty_type))


class Component(DictObject):

    fields = [
        Field("description", str),
        Field("extent", Extent),
        Field("displacement_type", displacementTypes),
        Field("uncertainty_type", uncertaintyTypes),
        Field("horizontal_uncertainty", float, min=0.0),
        Field("vertical_uncertainty", float, min=0.0),
        Field("spatial_model", SpatialModel),
        Field("time_function", TimeFunction.factory),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)
        if context is not None and context.check:
            self.spatial_model.check(self.displacment_type, self.uncertainty_type, context)
        self.epoch = None
        self.tfunc = None

    def _zero(self):
        return np.array([0.0, 0.0, 0.0, 0.0, 0.0])

    def setUncertainty(self, result):
        if self.uncertainty_type == "none":
            result[H_UNCERTAINTY_ORDINATE] = self.horizontal_uncertainty
            result[V_UNCERTAINTY_ORDINATE] = self.vertical_uncertainty
        elif self.uncertainty_type == "horizontal":
            result[V_UNCERTAINTY_ORDINATE] = self.vertical_uncertainty
        elif self.uncertainty_type == "vertical":
            result[H_UNCERTAINTY_ORDINATE] = self.horizontal_uncertainty

    def valueAt(self, x, y, epoch):
        if not self.extent.contains(x, y):
            result = self._zero()
        else:
            if epoch != self.epoch:
                self.tfunc = self.time_function.valueAt(epoch)
            if self.tfunc == 0.0:
                result = self._zero()
            else:
                result = self.spatial_model.valueAt(x, y) * self.tfunc
        self.setUncertainty(result)
        return result


class Authority(DictObject):
    fields = [
        Field("name", str, regex=r".*\S"),
        Field("url", str, regex=r".*\S"),
        Field("address", str, regex=r".*\S"),
        Field("email", str, regex=r"^.*\@.*"),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)


class Link(DictObject):

    fields = [
        Field("rel", str, regex=r"^(about|source|metadata|license)$"),
        Field("url", str, regex=r"^https?\:\/\/"),
        Field("type", str, regex=r"^(text\/html|application\/zip|application\/xml)$"),
        Field("title", str),
    ]

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)


class DeformationModel(DictObject):
    fields = [
        Field("file_type", str, regex=r"^deformation_model_master_file$"),
        Field("format_version", str, regex=r"1.0"),
        Field("name", str),
        Field("version", str),
        Field("publication_date", TimeFunction.parseTime),
        Field("license", str),
        Field("description", str),
        Field("authority", Authority),
        Field("links", [Link], optional=True),
        Field("source_crs", str, regex=r"^EPSG\:\d+$"),
        Field("target_crs", str, regex=r"^EPSG\:\d+$"),
        Field("definition_crs", str, regex=r"^EPSG\:\d+$"),
        Field("reference_epoch", TimeFunction.decimalYear),
        Field("uncertainty_reference_epoch", TimeFunction.decimalYear),
        Field("horizontal_offset_unit", str, regex=r"^metre$"),
        Field("vertical_offset_unit", str, regex=r"^metre$"),
        Field("horizontal_uncertainty_type", str, regex=r"^circular 95% confidence limit$"),
        Field("horizontal_uncertainty_unit", str, regex=r"^metre$"),
        Field("vertical_uncertainty_type", str, regex=r"^95% confidence limit$"),
        Field("vertical_uncertainty_unit", str, regex=r"^metre$"),
        Field("horizontal_offset_method", str, regex=r"^addition$"),
        Field("extent", Extent),
        Field("time_extent", TimeExtent),
        Field("components", [Component]),
    ]

    @staticmethod
    def LoadJson(sourcefile, check=False):
        with open(sourcefile, "r") as jsonf:
            value = json.load(jsonf)
        context = Context(sourcefile, check)
        model = DeformationModel(value, context)
        return model

    def __init__(self, value, context=None):
        DictObject.__init__(self, self.fields, value, context=context)


def main():
    import argparse

    parser = argparse.ArgumentParser()

    parser = argparse.ArgumentParser(description="Calculate deformation in a deformation model JSON file")
    parser.add_argument("json_file", help="Name of the deformation model JSON file")
    parser.add_argument("--point", nargs=3, help="lon,lat,epoch to evaluate")
    parser.add_argument("-i", "--point-file", help="Name of CSV file - columns id,lon,lat,epoch")
    parser.add_argument("-o", "--output-file", help="Name of CSV output file")
    parser.add_argument("-c", "--check", action="store_true", help="Check the deformation model")
    args = parser.parse_args(argv)
    model = DeformationModel.LoadJson(args.json_file)
    if args.point:
        lon = float(args.point[0])
        lat = float(args.point[1])
        epoch = TimeFunction.decimalYear(args.point[2])
        print(model.valueAt(lon, lat, epoch))
    if args.point_file:
        import csv

        csvi = csv.DictReader(open(args.point_file, "r"))
        for field in ("id", "lon", "lat", "epoch"):
            if field not in csvi.fieldnames:
                raise RuntimeError("Input file {0} is missing field {1}".format(args.point_file, field))
        if args.output_file:
            csvo = csv.writer(open(args.output_file, "w"))
        else:
            csvo = csv.writer(sys.stdout)
        csvo.writerow(["id", "lon", "lat", "epoch", "de", "dn", "du", "eh", "ev"])
        nrow = 1
        for row in csvi:
            nrow += 1
            try:
                id = row["id"]
                lon = float(row["lon"])
                lat = float(row["lat"])
                epochstr = row["epoch"]
                epoch = TimeFunction.decimalYear(epochstr)
                result = model.valueAt(lon, lat, epoch)
                output = [id, str(lon), str(lat), epochstr] + [str(r) for r in result]
                csvo.writerow(output)
            except Exeption as ex:
                print("Error in row {0}: {1}".format(nrow, ex))


if __name__ == "__main__":
    main()

