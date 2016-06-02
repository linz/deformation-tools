#!/usr/bin/python

from osgeo import ogr
from shapely import wkt, wkb
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import cascaded_union
from shapely import affinity
import math
import sys

# Polygons defining area over which model is required to meet standards...
# Tolerance is offset in degrees used to buffer and simplify model in order to 
# simplify calculations
# Land areas is a multipolygon defining the required coverage.  If it is empty then
# land is treated as infinite!

land_area_buffer=20000
land_area_tolerance=5000

def buffered_polygon( polygon, buffer, simplify=0.0 ):
    '''
    Buffer a polygon by an extents in metres.  Approximately scales to metres, applies buffer,
    and scales back. Returns a multipolygon
    '''

    if buffer > 0:
        minx,miny,maxx,maxy=polygon.bounds
        midx=(minx+maxx)/2.0
        midy=(miny+maxy)/2.0
        origin=[midx,midy,0.0]
        scale_factor=(1.0e-5/math.cos(math.radians(midy)),1.0e-5)
        polygon=affinity.scale(polygon,xfact=1.0/scale_factor[0],yfact=1.0/scale_factor[1],origin=origin)
        # Simplifying first to speed up process
        polygon=polygon.simplify(buffer/10.0)
        polygon=polygon.buffer(buffer)
        if simplify > 0:
            polygon=polygon.simplify(simplify)
        polygon=affinity.scale(polygon,xfact=scale_factor[0],yfact=scale_factor[1],origin=origin)
    return polygon

def create_land_areas( polygon_shapefile, extents_wktfile, buffer=land_area_buffer, tolerance=land_area_tolerance, min_points=0, verbose=False ):
    areas=[]
    driver=ogr.GetDriverByName('ESRI Shapefile')
    print("Loading land area definition from "+polygon_shapefile)
    datasource=driver.Open(polygon_shapefile,0)
    if datasource is None:
        raise RuntimeError('Cannot open land areas file '+polygon_file)
    layer=datasource.GetLayer()
    npoints=0
    nskip=0
    areas=[]
    for feature in layer:
        mp=wkb.loads(feature.GetGeometryRef().ExportToWkb())
        if type(mp) == Polygon:
            mp=[mp]
        for p in mp:
            npoints1=len(p.exterior.coords)
            if min_points and len(p.exterior.coords) < min_points:
                nskip += 1
                continue
            p=Polygon(p.exterior)
            p=buffered_polygon(p,buffer,tolerance)
            npoints2=len(p.exterior.coords)
            if type(p) == Polygon:
                p=[p]
            areas.extend(p)
            npoints += npoints2
            if verbose:
                print("Polygon: {0} points reduced to {1} points".format(npoints1,npoints2))

    if verbose:
        print("Skipped {0} polygons < {1} points".format(nskip,min_points))
    if areas:
        if verbose:
            print("Forming union of areas - total of {0} points in {1} polygons"
              .format(npoints,len(areas)))
        areas=MultiPolygon(areas)
        areas=areas.buffer(0)
        try:
            if verbose:
                print("Writing wkt file {0}".format(extents_wktfile))
            from shapely.wkt import dumps
            with open(extents_wktfile,"w") as laf:
                laf.write(dumps(areas))
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
    parser.add_argument("-p","--min-points",type=int,default=0,help="Tolerance for shape simplification")
    parser.add_argument("-f","--force",action="store_true",help="Force building even if current version already built")
    parser.add_argument("-v","--verbose",action="store_true",help="Print output during calculations")

    args=parser.parse_args()
    sfile=args.nzpoly_shapefile
    ofile=args.output_wktfile
    verbose=args.verbose

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
    if verbose:
        print("Building coastline polygon extents")
        print("Buffering by {0} degrees".format(args.buffer))
        print("Simplification tolerance {0} degrees".format(args.tolerance))
    create_land_areas(sfile,ofile,buffer=args.buffer,tolerance=args.tolerance,
                      min_points=args.min_points,verbose=verbose)
