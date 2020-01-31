#!/usr/bin/python3
###############################################################################
# $Id$
#
#  Project:  PROJ-data
#  Purpose:  Convert deformation model .CSV files to GTIFF format
#  Author:   Chris Crook <ccrook@linz.govt.nz>
#
###############################################################################
#  Copyright (c) 2019, Land Information New Zealand (www.linz.govt.nz)
#
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
###############################################################################
# This software adapted from ntv2_to_gtiff.py by Even Rouault

synopsis = """
Creates a GeoTIFF formatted nest grid file based on a set of CSV files 
formatted as defined for the LINZ deformation format (url).

The files to include are defined in an input file which includes the 
interpolation crs defined as EPSG:code, the number of rows in the grid,
the grid dimensions (displacement:error), and for each grid the number of
rows and columns.

This is formatted as a JSON file with structure:

  {
  "crs": "EPSG:4959",
  "type": "3d:none",
  "copyright": "Land Information New Zealand (2013): Released under Create Commons Attribution 4.0 International",
  "description": "Darfield earthquake",
  "version": "YYYYMMDD",
  "grids": [
    {
      "filename": "patch_c1_20100904/grid_c1_L1.csv",
      "extent": [168.1,176.2,-46.75,-40.375],
      "size": [55,52]
    },
    ...
    ]
    }
"""

from osgeo import gdal
from osgeo import osr
from cloud_optimize_gtiff import generate_optimized_file
import argparse
import datetime
import csv
import json
import math
import os
import struct
import numpy as np

def get_args(args=None):
    parser = argparse.ArgumentParser(description="Convert JSON csv file list to GeoTIFF nested grid.")
    parser.add_argument("source", help="Source JSON file list")
    parser.add_argument("dest", help="Destination GeoTIFF file")
    parser.add_argument("-m", "--model-dir", help="Base directory of deformation model")
    parser.add_argument("-c","--compact-metadata",action="store_true",help="Reduce size of metadata in GeoTIFF directory")
    parser.add_argument(
        "--uint16-encoding",
        dest="uint16_encoding",
        action="store_true",
        help="Use uint16 storage with linear scaling/offseting",
    )
    return parser.parse_args()

def isSubgrid( grid1, grid2 ):
    # Test if the extents of grid1 are fully within grid2
    ext1 = grid1['extent']
    ext2 = grid2['extent']
    return ( ext1[0] >= ext2[0] and ext1[2] <=ext2[2] and ext1[1] >= ext2[1] and ext1[3] <= ext2[3])

def create_unoptimized_file(sourcefilename, tmpfilename, args):
    gridspec = json.loads(open(sourcefilename).read())
    subdatsets = gridspec['grids']

    # Determine source CSV fields and target bands for GeoTIFF file

    displacement_type, error_type = gridspec['type'].upper().split(':')
    csv_fields = []
    bands = []
    if displacement_type in ('HORIZONTAL', '3D'):
        csv_fields.append('de')
        csv_fields.append('dn')
        bands.append('east_offset')
        bands.append('north_offset')
    if displacement_type in ('VERTICAL','3D'):
        csv_fields.append('du')
        bands.append('vertical_offset')
    if error_type in ('HORIZONTAL', '3D'):
        csv_fields.append('eh')
        bands.append('horizontal_uncertainty')
    if error_type in ('VERTICAL','3D'):
        csv_fields.append('ev')
        bands.append('vertical_uncertainty')
    nbands = len(bands)
    usecols=list(range(2,len(bands)+2))

    # Prepare the subgrids for inclusion into the data set
    # Set the name for each grid, check the source file exists ...

    modeldir=args.model_dir
    subgrids = {}
    for nsub,subds in enumerate(subdatsets):
        name = subds['filename']
        name = name.replace('/', '_')
        if name.endswith('.csv'):
            name = name[:-4]
        subds['name'] = name
        filename=subds['filename']
        if modeldir is not None:
            filename = os.path.join(modeldir, 'model', filename)
        if not os.path.exists(filename):
            raise RuntimeError('Source file {0} is missing'.format(filename))
        try:
            csvr = csv.reader(open(filename))
            source_fields = next(csvr)
            csvr = None
        except:
            raise RuntimeError("Error opening CSV file {0}".format(filename))
        if source_fields[:2] != ['lon', 'lat'] or source_fields[2:] != csv_fields:
            raise RuntimeError("Invalid fields in CSV {0}: should be lon, lat, {1}"
                    .format(filename, ', '.join(csv_fields)))
        subds['sourcefile'] = filename
        subds['children'] = []
        parent=None
        while nsub > 0:
            nsub -= 1
            if isSubgrid(subds, subdatsets[nsub]):
                parent = subdatsets[nsub]
                break
        subds['parent'] = parent
        if parent is not None:
            parent['children'].append(subds)

    # Compile the grids
 
    for idx_ifd, subds in enumerate(subdatsets):

        # Read the data from the CSV file

        csvfile = subds['sourcefile']
        data = np.genfromtxt(csvfile, names=True, delimiter=",", usecols=usecols)
        # Shape is (nrows=nlat values,ncols = n lon values)
        data=data.reshape((subds['size'][1],subds['size'][0]))

        size=subds['size']
        tmp_ds = gdal.GetDriverByName("GTiff").Create(
            "/vsimem/tmp",
            size[0],
            size[1],
            nbands,
            gdal.GDT_Float32 if not args.uint16_encoding else gdal.GDT_UInt16,
        )

        if idx_ifd == 0:
            description = gridspec['description']
            copyright = gridspec['copyright']
            version=gridspec['version']
            griddate="{0}:{1}:{2} 00:00:00".format(version[:4],version[4:6],version[6:8])
            tmp_ds.SetMetadataItem("TIFFTAG_IMAGEDESCRIPTION", description)
            tmp_ds.SetMetadataItem("TIFFTAG_COPYRIGHT", copyright)
            tmp_ds.SetMetadataItem("TIFFTAG_DATETIME", griddate)

        if idx_ifd == 0 or not args.compact_metadata:
            tmp_ds.SetMetadataItem("TYPE", "DEFORMATION_MODEL")
            tmp_ds.SetMetadataItem("DISPLACEMENT_TYPE", displacement_type)
            tmp_ds.SetMetadataItem("UNCERTAINTY_TYPE", error_type)
            for i, band in enumerate(bands):
                tmp_ds.GetRasterBand(i+1).SetDescription(band)
                tmp_ds.GetRasterBand(i+1).SetUnitType("metre")

        # Calculate the GeoTransform
        ext = subds['extent']
        dx = (ext[2] - ext[0]) / (size[0] - 1)
        dy = (ext[3] - ext[1]) / (size[1] - 1)
        transform = [ext[0] - dx / 2.0, dx, 0, ext[1] - dy / 2.0, 0, dy]

        src_crs = osr.SpatialReference()
        src_crs.SetFromUserInput(gridspec['crs'])
        tmp_ds.SetSpatialRef(src_crs)
        tmp_ds.SetGeoTransform(transform)
        tmp_ds.SetMetadataItem("AREA_OR_POINT", "Point")

        tmp_ds.SetMetadataItem("grid_name", subds["name"])
        if subds['parent'] is not None:
            tmp_ds.SetMetadataItem("parent_grid_name", subds['parent']['name'])
        if len(subds['children']) > 0:
            tmp_ds.SetMetadataItem("number_of_nested_grids", str(len(subds['children'])))

        for i, field in enumerate(csv_fields):
            bdata = data[field].copy()
            # Might be useful to reconsider how to handle scaling if we 
            # decide to employ it. eg to try 1.0e-6, 1.0e-5 ... to better
            # round to actual values.  Worst case (?) 10m earth movement
            # would be rounded to 1mm.  Also check rounding of source data,
            # eg if already rounded to 0.0001.
            if args.uint16_encoding:
                bmin = bdata.min()
                bmax = bdata.max()
                scale = (bmax - bmin) / 65535
                bdata = (bdata - bmin) / scale
                tmp_ds.GetRasterBand(i+1).SetOffset(min)
                tmp_ds.GetRasterBand(i+1).SetScale(scale)           
            tmp_ds.GetRasterBand(i+1).WriteArray(data)

        options = [
            "PHOTOMETRIC=MINISBLACK",
            "COMPRESS=DEFLATE",
            "PREDICTOR=3" if not args.uint16_encoding else "PREDICTOR=2",
            "INTERLEAVE=BAND",
            "GEOTIFF_VERSION=1.1",
        ]
        if tmp_ds.RasterXSize > 256 and tmp_ds.RasterYSize > 256:
            options.append("TILED=YES")
        else:
            options.append("BLOCKYSIZE=" + str(tmp_ds.RasterYSize))
        if gdal.VSIStatL(tmpfilename) is not None:
            options.append("APPEND_SUBDATASET=YES")

        assert gdal.GetDriverByName("GTiff").CreateCopy(tmpfilename, tmp_ds, options=options)


# def check(sourcefilename, destfilename, args):
#     gridspec = gdal.Open(sourcefilename)
#     assert gridspec.GetDriver().ShortName in ("NTv2", "NTv1", "CTable2")
#     src_subdatsets = [(sourcefilename, None)]
#     src_subdatsets += gridspec.GetSubDatasets()

#     dst_ds = gdal.Open(destfilename)
#     dst_subdatsets = dst_ds.GetSubDatasets()
#     if not dst_subdatsets:
#         dst_subdatsets = [(destfilename, None)]

#     assert len(src_subdatsets) == len(dst_subdatsets)
#     for src_subds, dst_subds in zip(src_subdatsets, dst_subdatsets):
#         gridspec = gdal.Open(src_subds[0])
#         dst_ds = gdal.Open(dst_subds[0])
#         if not args.uint16_encoding:
#             for i in range(min(gridspec.RasterCount, dst_ds.RasterCount)):
#                 data = gridspec.GetRasterBand(i + 1).ReadRaster(buf_type=gdal.GDT_Float32)

#                 if gridspec.GetDriver().ShortName == "CTable2":
#                     nvalues = gridspec.RasterXSize * gridspec.RasterYSize
#                     out_data = b""
#                     # From radian to arc-seconds
#                     for v in struct.unpack("f" * nvalues, data):
#                         out_data += struct.pack("f", v / math.pi * 180.0 * 3600)
#                     data = out_data

#                 if i + 1 == 2 and args.positive_longitude_shift_value == "east":
#                     nvalues = gridspec.RasterXSize * gridspec.RasterYSize
#                     out_data = b""
#                     for v in struct.unpack("f" * nvalues, data):
#                         out_data += struct.pack("f", -v)
#                     data = out_data
#                 assert dst_ds.GetRasterBand(i + 1).ReadRaster() == data
#         else:
#             import numpy as np

#             for i in range(min(gridspec.RasterCount, dst_ds.RasterCount)):
#                 src_data = gridspec.GetRasterBand(i + 1).ReadAsArray()
#                 dst_data = dst_ds.GetRasterBand(i + 1).ReadAsArray()
#                 offset = dst_ds.GetRasterBand(i + 1).GetOffset()
#                 scale = dst_ds.GetRasterBand(i + 1).GetScale()
#                 dst_data = dst_data * scale + offset
#                 if i + 1 == 2 and args.positive_longitude_shift_value == "east":
#                     dst_data = -dst_data
#                 max_error = np.max(abs(dst_data - src_data))
#                 assert max_error <= 1.01 * scale / 2, (max_error, scale / 2)

def build_deformation_gtiff(source,target,args):
    tmpfilename = target + ".tmp"
    gdal.Unlink(tmpfilename)
    create_unoptimized_file(source, tmpfilename, args)
    generate_optimized_file(tmpfilename, target)
    #check(source, target, args)
    gdal.Unlink(tmpfilename)

if __name__ == "__main__":
    args = get_args()
    build_deformation_gtiff(args.source,args.dest,args)
