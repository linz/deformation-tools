#!/usr/bin/python

description='''
Script to calculate the zones in which orders are downgraded.  Takes as input
a file defining zones with records of the format
zone_code "grid" gridfile column tolerance [buffer=#] [simplify=#] [error_factor=#]
zone_code "wkt" wkt_file [buffer=#] [simplify=#]
zone_code "fault_model" fault_wkt_file [buffer=#] [minwidth=#]

Values in ".." are literal text (entered without the "s)
zone_code is a code defining the zone being defined, which in turn defines the orders
that are downgraded (see downgrade definition file)

error_factor is a factor by which the value from the grid is multiplied in order to 
determine the value that will be compared with the tolerance.  The assumption is that 
the error of the calculated value (ie the error of the deformation model which impacts
on the order of the mark) is roughly proportional to the value itself.  The values used
are the relative errors, the relative vector error for the horizontal component and the 
tilt for the vertical component.

Buffer and simplification values are used to transform shapes to simpler objects
to improve efficiency in using the resulting zones.
'''

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
from shapely.geometry import Polygon, MultiPolygon, MultiPoint
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

def bufferPoly( pgn, buffer, simplify=0.0 ):
    if buffer <=0 and simplify <= 0:
        return pgn
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

def wktPolys( wktfile, buffer=0.0, simplify=0.0, wktonly=False ):
    pgnwkt=[]
    with open( wktfile ) as wktf:
        if wktonly:
            pgnwkt=wktf.readlines()
        else:
            csvf=csv.DictReader(wktf)
            for r in csvf:
                wkts=r["WKT"]
                pgnwkt.append(wkts)
    geoms=[]
    for wkts in pgnwkt:
        geom=wkt.loads(wkts)
        if not geom:
            continue
        geom=bufferPoly(geom,buffer,simplify)
        geom=exteriors(geom)
        geoms.extend(geom)
    return compilePolys(geoms)

def faultPolys( wktfile, buffer, minwidth=None ):
    crdre=r'(?:\-?\d+(?:\.\d+)?\s+\-?\d+(?:\.\d+)\s+\-?\d+(?:\.\d+))'
    rectre=re.compile(r'\(('+crdre+r'(\s*\,\s*'+crdre+r'){4})\)')
    minwidth=buffer/10
    geoms=[]
    with open(wktfile) as wktf:
        header=next(wktf)
        for segdef in wktf:
            segdef=segdef.strip()
            if segdef == '':
                continue
            match=rectre.search(segdef)
            if not match:
                print("Invalid fault segment: {0}".format(segdef))
                continue
            crddef=match.group(1)
            crds=[]
            for crd in crddef.split(','):
                crds.append([float(x) for x in crd.split()])

            line1=[crds[0],crds[1]]
            line2=[crds[3],crds[2]]
            z1=-crds[0][2]
            z2=-crds[3][2]
            if z1 >= buffer and z2 >= buffer:
                continue
            if z1 > buffer or z2 > buffer:
                if z1 > z2:
                    z1,z2=z2,z1
                    line1,line2=line2,line1
                intp=lambda x1, x2: ((z1-buffer)*x2+(buffer-z2)*x1)/(z1-z2)
                line2=[
                    [intp(line1[0][0],line2[0][0]),intp(line1[0][1],line2[0][1])],
                    [intp(line1[1][0],line2[1][0]),intp(line1[1][1],line2[1][1])]
                    ]
                z2=buffer
            w1=math.sqrt(max(0.0,buffer*buffer-z1*z1))
            w2=math.sqrt(max(0.0,buffer*buffer-z2*z2))
            x0=(line1[0][0]+line1[1][0]+line2[0][0]+line2[1][0])/4.0
            y0=(line1[0][1]+line1[1][1]+line2[0][1]+line2[1][1])/4.0
            xf=100000.0*math.cos(math.radians(y0))
            yf=100000.0
            toxy=lambda crd: [(crd[0]-x0)*xf,(crd[1]-y0)*yf]
            toll=lambda crd: [crd[0]/xf+x0,crd[1]/yf+y0]
            line1=[toxy(c) for c in line1]
            line2=[toxy(c) for c in line2]
            dw=math.hypot(line1[0][0]-line2[0][0],line1[0][1]-line2[0][1])
            if max(w1*2,w2*2,w1+w2+dw) < minwidth:
                continue
            dl=math.hypot(line1[0][0]-line1[1][0],line1[0][1]-line1[1][1])
            dx=(line1[1][0]-line1[0][0])/dl
            dy=(line1[1][1]-line1[0][1])/dl
            crds=[]
            for w,l in zip((w1,w2),(line1,line2)):
                for f,p in zip((-1,1),l):
                    crds.append(toll([p[0]+f*w*dx-w*dy,p[1]+f*w*dy+w*dx]))
                    crds.append(toll([p[0]+f*w*dx+w*dy,p[1]+f*w*dy-w*dx]))
            geoms.append(MultiPoint(crds).convex_hull)
    return compilePolys(geoms, True)

def compilePolys( geoms, keep_holes=False ):
    union=None
    ntries=5

    if keep_holes:
        ext=lambda x: [x] if type(x) == Polygon else x
    else:
        ext=exteriors
    
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

    result=[]
    if union != None:
        result=ext(union)
        for p in failed:
            result.extend(ext(p))
    return result

def main():
    parser=argparse.ArgumentParser(description=description)
    parser.add_argument('zone_def_file',help='File defining components to include in zones')
    parser.add_argument('zone_file',help='Output zone definition file')
    parser.add_argument('parameters',nargs='*',help='Parameters substituted into file')
    parser.add_argument('-b','--buffer',type=float,default=2000.0,help="Buffer applied to final polygon")
    parser.add_argument('-s','--simplify',type=float,default=1000.0,help="Simplification tolerance applied to final polygon")
    parser.add_argument('-d','--decimal-places',type=int,default=5,help="Number of decimal places in WKT")
    parser.add_argument('-w','--wkt-only',action='store_true',help="Output file is simple list of wkt shapes")
    parser.add_argument('-m','--multipolygon',action='store_true',help="Polygons for each zone are compiled to multipolygon")
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
                elif len(parts) >= 3 and parts[1] == "wkt" or parts[1] == "fault_model":
                    zone_code,skip,wktfile=parts[:3]
                    buffer=2000.0
                    simplify=1000.0
                    minwidth=None
                    for p in parts[3:]:
                        if '=' not in p:
                            raise RuntimeError("Invalid component "+p)
                        item,value=p.split('=',1)
                        value=float(value)
                        if item=='buffer':
                            buffer=value
                        elif item=='simplify':
                            simplify=value
                        elif item=='minwidth':
                            minwidth=value
                        else:
                            raise RuntimeError("Invalid item "+item)
                    if parts[1] == "wkt" or parts[1] == "wktonly":
                        geoms=wktPolys(wktfile,buffer,simplify,parts[1] == "wktonly")
                    else:
                        geoms=faultPolys(wktfile,buffer,minwidth)
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
    ndp=args.decimal_places
    template="{1}\n" if args.wkt_only else "{0}|Zone {0}|{1}\n"
    with open(args.zone_file,'w') as zf:
        for z in zones:
            polys=zones[z]
            polys=[bufferPoly(p,buffer,simplify) for p in polys]
            if len(polys) > 1:
                polys=compilePolys(zones[z])
            if args.multipolygon:
                mp=MultiPolygon(polys)
                if not mp.is_valid:
                    mp=bufferPoly(mp,0.0)
                if not mp.is_valid:
                    mp=bufferPoly(mp,1.0)
                polys=[mp]
            for p in polys:
                if not p.is_valid:
                    print("Skipping invalid zone {0} poly".format(zone_code))
                    continue
                pwkt=p.wkt
                pwkt=re.sub(r'(\.\d{{{0}}})\d+'.format(ndp),r'\1',pwkt)
                pgn=wkt.loads(pwkt)
                if not pgn.is_valid:
                    pgn=pgn.buffer(0.0)
                    pwkt=pgn.wkt
                zf.write(template.format(z,pwkt))

if __name__=="__main__":
    main()
