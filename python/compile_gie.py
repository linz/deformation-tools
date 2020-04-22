#!/usr/bin/python3
#
# Script to compile a PROJ gie script to validate a PROJ NZGD2000
# implementation.  This uses the published model as a source from which
# to build the gie file.  It looks at the model components to generate
# set of points and epochs that will ensure that each component of the
# model is tested.
#
# The script downloads the official NZGD2000 model www.geodesy.linz.govt.nz,
# unpacks it in a temporary directory, and then uses the LINZ NZGD2000 python
# module to interrogate it and generate the set of test points and expected
# outputs.

import sys
import os.path
import argparse
import random
import tempfile
import requests
import zipfile
import pyproj
from LINZ.DeformationModel import Model, Time

model_url_template = "https://www.geodesy.linz.govt.nz/download/nzgd2000/nzgd2000_deformation_{version}_full.zip"
gie_file_template = "nzgd2000-{version}.gie"
gie_inv_template = "nzgd2000-inv-{version}.gie"
proj_template = "+proj=defmodel +model=nzgd2000-20180701.json"
pointspergrid = 5
tolerance = 0.001
gie_header = """
============================================================
Test: NZGD2000 version {version} deformation model
============================================================
"""
gie_inv_header = """
============================================================
Test: NZGD2000 version {version} deformation model (inverse)
============================================================
"""
proj_nztm = "+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("version", help="Version of deformation model to test")
    parser.add_argument("-n", "--npoints", type=int, default=pointspergrid, help="Number of test points per source grid")
    parser.add_argument("-i", "--inverse", action="store_true", help="Build a gie file for the inverse transformation")
    parser.add_argument("-s", "--seed", action="store_true", help="Seed random number generator for repeatable output")
    args = parser.parse_args()
    compileGie(args.version, pointspergrid=args.npoints, inverse=args.inverse, useseed=args.seed)


def componentTestPoints(component, pointspergrid):
    maxzero = max(1, int(0.2 * pointspergrid))
    maxreject = pointspergrid * 10
    tf = component.timeFunction
    t0 = Time.Time.Parse(tf.time0)
    t1 = Time.Time.Parse(tf.time1)
    y0 = t0.asYear()
    y1 = t1.asYear()
    times = [y0 - 0.1]
    if y1 > y0:
        times.append(y0 * 0.75 + y1 * 0.25)
    times.append(y1 + 0.0)
    models = []
    for model in component.spatialModel.models():
        model.load()
        models.append(model)
    prioritymodels = []
    points = []
    randht = lambda: random.uniform(0.0, 1000.0)
    for m in models:
        ntest = 0
        nzero = 0
        nreject = 0
        useall = False
        randpt = lambda: [random.uniform(m.min_lon, m.max_lon), random.uniform(m.min_lat, m.max_lat)]
        while ntest < pointspergrid:
            x, y = randpt()
            if not useall:
                reject = False
                for p in prioritymodels:
                    if p.containsPoint(x, y):
                        reject = True
                        break
                if not reject:
                    denu, skip = m.calcDeformation(x, y)
                    zero = sum((d * d for d in denu)) <= 0.000001
                    if zero and nzero < maxzero:
                        reject = True
                if reject:
                    nreject += 1
                    if nreject > maxreject:
                        useall = True
                    continue
            ntest += 1
            for t in times:
                points.append([x, y, randht(), t])
        if useall:
            print("Warning: Gave up trying to find enough good test points for {0}".format(m))
    return points


def loadModel(modelurl, targetdir):
    if not os.path.isdir(targetdir):
        raise RuntimeError("Invalid target directory {0}".format(targetdir))
    modelzip = os.path.join(targetdir, "model.zip")
    response = requests.get(modelurl)
    if not response.ok:
        raise RuntimeError("Cannot retrieve {0}".format(modelurl))
    with open(modelzip, "wb") as zf:
        zf.write(response.content)
    response.close()
    with zipfile.ZipFile(modelzip, "r") as mz:
        mz.extractall(targetdir)
    return Model.Model(os.path.join(targetdir, "model"))


def writeGieFile(giefile, version, inverse, tests):
    nztm = pyproj.Proj(proj_nztm)
    projstring = proj_template.replace("{version}", version)
    inv = " +inv" if inverse else ""
    projfull = "+proj=pipeline +step +inv {0} +step{2} {1} +step {0}".format(proj_nztm, projstring, inv)
    with open(giefile, "w") as gf:
        gf.write((gie_inv_header if inverse else gie_header).replace("{version}", version))
        gf.write("<gie>\n")
        gf.write("operation {0}\n".format(projfull))
        gf.write("tolerance {0}\n\n".format(tolerance))
        for test, result in tests:
            crdin = test[:4]
            crdout = result[:3]
            crdin[:2] = nztm(crdin[0], crdin[1])
            crdout[:2] = nztm(crdout[0], crdout[1])
            crdout.append(crdin[3])
            gf.write("accept {0:.4f} {1:.4f} {2:.4f} {3:.3f}\n".format(*crdin))
            gf.write("except {0:.4f} {1:.4f} {2:.4f} {3:.3f}\n".format(*crdout))
            gf.write("roundtrip 1\n\n")
        gf.write("</gie>\n")


def compileGie(version, pointspergrid=pointspergrid, inverse=False, useseed=False):
    model_url = model_url_template.replace("{version}", version)
    gie_file = (gie_inv_template if inverse else gie_file_template).replace("{version}", version)
    if useseed:
        random.seed(42)
    testpoints = []
    with tempfile.TemporaryDirectory() as tempdir:
        model = loadModel(model_url, tempdir)
        for c in model.components():
            testpoints.extend(componentTestPoints(c, pointspergrid))
        tests = [(p, model.applyTo(*p, subtract=inverse)) for p in testpoints]
        writeGieFile(gie_file, version, inverse, tests)


if __name__ == "__main__":
    main()
