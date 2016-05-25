#!/usr/bin/python

from osgeo import ogr
from shapely import wkt, wkb
from shapely.ops import unary_union, transform
from shapely.geometry import MultiPolygon, Polygon
import math
import sys

# Polygons defining area over which model is required to meet standards...
# Tolerance is offset in degrees used to buffer and simplify model in order to 
# simplify calculations
# Land areas is a multipolygon defining the required coverage.  If it is empty then
# land is treated as infinite!

land_area_buffer=20000
land_area_tolerance=5000

def load_land_areas( polygon_shapefile, extents_wktfile, buffer=land_area_buffer, tolerance=land_area_tolerance ):
    areas=[]
    driver=ogr.GetDriverByName('ESRI Shapefile')
    print("Loading land area definition from "+polygon_shapefile)
    datasource=driver.Open(polygon_shapefile,0)
    if datasource is None:
        raise RuntimeError('Cannot open land areas file '+polygon_file)
    layer=datasource.GetLayer()
    for feature in layer:
        mp=wkb.loads(feature.GetGeometryRef().ExportToWkb())
        if type(mp) != MultiPolygon:
            mp=[mp]
        for p in mp:
            print("Polygon: {0} points".format(len(list(p.exterior.coords))))
            p=Polygon(p.exterior)

            minx,miny,maxx,maxy=p.bounds
            midx=(minx+maxx)/2.0
            midy=(miny+maxy)/2.0
            yscale=100000
            xscale=100000*math.cos(math.radians(midy))
            p=transform(lambda x,y,z=None: ((x-midx)*xscale,(y-midy)*yscale,0.0),p)
            p=p.buffer(buffer,3).simplify(tolerance)
            p=transform(lambda x,y,z=None: (x/xscale+midx,y/yscale+midy,0.0),p)
            areas.append(p)
    if areas:
        land_areas=unary_union(areas)
        try:
            from shapely.wkt import dumps
            with open(extents_wktfile,"w") as laf:
                laf.write(dumps(land_areas))
        except:
            pass

if __name__=="__main__":
    import argparse
    import os.path
    parser=argparse.ArgumentParser("Build buffered land extents wkt file")
    parser.add_argument("nzpoly_shapefile",help="NZ island polygons shape file")
    parser.add_argument("output_wktfile",help="Output buffered island WKT file")
    parser.add_argument("-b","--buffer",type=float,default=land_area_buffer,help="Buffer applied to island polygons")
    parser.add_argument("-t","--tolerance",type=float,default=land_area_tolerance,help="Tolerance for shape simplification")
    parser.add_argument("-f","--force",action="store_true",help="Force building even if current version already built")

    args=parser.parse_args()
    sfile=args.nzpoly_shapefile
    ofile=args.output_wktfile

    if not os.path.exists(sfile):
        print "NZ polygon shape file "+sfile+" does not exist - aborting!"
        sys.exit()

    if not args.force and os.path.exists(ofile) and os.path.getmtime(ofile) > os.path.getmtime(sfile):
        try:
            with open(ofile) as laf:
                pgnwkt=laf.read()
                land_areas=wkt.loads(pgnwkt)
                print( "Using existing wkt land extents in "+ofile)
            sys.exit()
        except Exception as ex:
            print("Cannot read existing polygon WKT in "+ofile)
    print("Building coastline polygon extents")
    print("Buffering by {0} degrees".format(args.buffer))
    print("Simplification tolerance {0} degrees".format(args.tolerance))
    load_land_areas(sfile,ofile,buffer=args.buffer,tolerance=args.tolerance)
