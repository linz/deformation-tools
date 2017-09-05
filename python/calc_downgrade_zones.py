#!/usr/bin/python
#
# Script to calculate the zones in which orders are downgraded.  Takes as input
# a file defining zones with records of the format
# zone_code "grid" gridfile column tolerance [buffer=#] [simplify=#] [error_factor=#]
# zone_code "wkt" wkt_file [buffer=#] [simplify=#]
#
# Values in ".." are literal text (entered without the "s)
# zone_code is a code defining the zone being defined, which in turn defines the orders
# that are downgraded (see downgrade definition file)
# 
# error_factor is a factor by which the value from the grid is multiplied in order to 
# determine the value that will be compared with the tolerance.  The assumption is that 
# the error of the calculated value (ie the error of the deformation model which impacts
# on the order of the mark) is roughly proportional to the value itself.  The values used
# are the relative errors, the relative vector error for the horizontal component and the 
# tilt for the vertical component.
#
# Buffer and simplification values are used to transform shapes to simpler objects
# to improve efficiency in using the resulting zones.
#
# The defined shapes are all assumed to be polygons.  For each zone code the components are 
# combined using a spatial union (likely to generate a multipolygon), then split into individual
# polygons for the final WKT definition.  Only the exterior of each polygon is used - holes
# are ignored.  Splitting into individual polygons is assumed to be more efficient in determining
# the value at any location - depends on the efficiency of intersecting a point with a multipolygon
# compared with the efficiency of first selecting the polygons to consider using the spatial 
# index.

import argparse
import math
import csv
import re
from shapely.geometry import Polygon, MultiPolygon
from shapely import affinity
from shapely import wkt
import defgrid

def exteriors( shape ):
    if type(shape) is Polygon:
        shape=[shape]
    elif type(shape) is not MultiPolygon:
        raise RuntimeError('Can only deal with Polygons and MultiPolygons')
    exteriors=[]
    for p in shape:
        exteriors.append(Polygon(p.exterior))
    return exteriors

def bufferPoly( pgn, buffer, simplify ):
    if buffer <=0 and simplify <= 0:
        return
    c=pgn.centroid
    yscl=100000.0
    xscl=yscl*math.cos(math.radians(c.y))
    pgn=affinity.scale(pgn,xfact=xscl,yfact=yscl,origin=c)
    if buffer > 0:
        pgn=pgn.buffer(buffer,resolution=4)
    if simplify > 0:
        pgn=pgn.simplify(simplify)
    pgn=affinity.scale(pgn,xfact=1.0/xscl,yfact=1.0/yscl,origin=c)
    return pgn

def gridPolys( gridfile, column, tolerance, buffer=0.0, simplify=0.0 ):
    g=defgrid.DeformationGrid(gridfile)
    pgns=g.regionsExceedingLevel( column, tolerance )
    pgns=[bufferPoly(p,buffer,simplify) for p in pgns]
    return pgns

def wktPolys( wktfile, buffer=0.0, simplify=0.0 ):
    with open( wktfile ) as wktf:
        csvf=csv.DictReader(wktf)
        geoms=[]
        for r in csvf:
            wkts=r["WKT"]
            geom=wkt.loads(wkts)
            geom=bufferPoly(geom,buffer,simplify)
            geom=exteriors(geom)
            geoms.extend(geom)
        return compilePolys(geoms)

def compilePolys( geoms ):
    union=None
    ntries=5
    
    while ntries > 0:
        ntries -= 1
        failed=[]
        for g in geoms:
            if union is None:
                union=g
            else:
                try: 
                    u=union.union(g)
                    if not u.is_valid:
                        failed.append(g)
                    else:
                        union=u
                except:
                    failed.append(g)
        if not failed:
            break

    result=exteriors(union)
    for p in failed:
        result.extend(exteriors(p))
    return result

def main():
    parser=argparse.ArgumentParser(description='Compile downgrade polygon zones')
    parser.add_argument('zone_def_file',help='File defining components to include in zones')
    parser.add_argument('zone_file',help='Output zone definition file')
    parser.add_argument('parameters',nargs='*',help='Parameters substituted into file')
    parser.add_argument('-b','--buffer',type=float,default=2000.0,help="Buffer applied to final polygon")
    parser.add_argument('-s','--simplify',type=float,default=1000.0,help="Simplification tolerance applied to final polygon")
    args=parser.parse_args()

    params={}
    for p in args.parameters:
        if '=' not in p:
            raise RuntimeError("Invalid command line parameter "+p)
        k,v=p.split('=',1)
        params[k]=v

    def replace_param(m):
        key=m.group(1)
        if k in params:
            return params[k]
        if k in sys.environ:
            return sys.environ[k]
        raise RuntimeError('Undefined parameter $\{{0}\} in configuration'.format(k))

    zones={}
    with open(args.zone_def_file) as zdf:
        for l in zdf:
            if re.match(r'^\s*(\#|$)',l):
                continue
            l=re.sub(r'\$\{(\w+)\}',replace_param,l)
            try:
                parts=l.split()
                if len(parts) >= 5 and parts[1] == "grid":
                    zone_code,skip,gridfile,column,tolerance=parts[:5]
                    tolerance=float(tolerance) 
                    buffer=2000.0
                    simplify=1000.0
                    factor=1.0
                    for p in parts[5:]:
                        if '=' not in p:
                            raise RuntimeError("Invalid component "+p)
                        item,value=p.split('=',1)
                        value=float(value)
                        if item=='buffer':
                            buffer=value
                        elif item=='simplify':
                            simplify=value
                        elif item=='error_factor':
                            factor=value
                        else:
                            raise RuntimeError("Invalid item "+item)
                    geoms=gridPolys(gridfile,column,tolerance/factor,buffer,simplify)
                elif len(parts) >= 3 and parts[1] == "wkt":
                    zone_code,skip,wktfile=parts[:3]
                    for p in parts[3:]:
                        if '=' not in p:
                            raise RuntimeError("Invalid component "+p)
                        item,value=p.split('=',1)
                        value=float(value)
                        if item=='buffer':
                            buffer=value
                        elif item=='simplify':
                            simplify=value
                        else:
                            raise RuntimeError("Invalid item "+item)
                    geoms=wktPolys(wktfile,buffer,simplify)
                else:
                    raise RuntimeError('Invalid data, not grid or shape definition')
                if zone_code not in zones:
                    zones[zone_code]=[]
                zones[zone_code].extend(geoms)
            except Exception as ex:
                raise
                msg=ex.message
                raise RuntimeError("Error in line: {0}\n{1}".format(l,ex.message))

    buffer=args.buffer
    simplify=args.simplify
    with open(args.zone_file,'w') as zf:
        for z in zones:
            polys=compilePolys(zones[z])
            polys=[bufferPoly(p,buffer,simplify) for p in polys]
            for p in compilePolys(polys):
                if not p.is_valid:
                    print "Skipping invalid zone {0} poly".format(zone_code)
                    continue
                zf.write("{0}|Zone {0}|{1}|\n".format(z,p.wkt))

if __name__=="__main__":
    main()
