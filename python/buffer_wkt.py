#!/usr/bin/python
#
# Approximately buffer a lat/lon WKT by a specified distance.

import sys
import shapely
import math
from shapely.wkt import loads
from shapely.ops import transform
from shapely.geos import ReadingError

if len(sys.argv) != 4:
    print "Syntax: buffer_wkt lat_lon_wkt_file buffer_distance_m output_wkt_file"
    sys.exit()

src_wkt, buffer, output_wkt=sys.argv[1:]
buffer=float(buffer)

with open(src_wkt) as srcf, open(output_wkt,'w') as outputf:
    for l in srcf:
        try:
            geom=loads(l)
            minx,miny,maxx,maxy=geom.bounds
            midx=(minx+maxx)/2.0
            midy=(miny+maxy)/2.0
            yscale=100000
            xscale=100000*math.cos(math.radians(midy))
            geom=transform(lambda x,y,z=None: ((x-midx)*xscale,(y-midy)*yscale,0.0),geom)
            geom=geom.buffer(buffer,resolution=3)
            geom=transform(lambda x,y,z=None: (x/xscale+midx,y/yscale+midy,0.0),geom)
            outputf.write(geom.wkt)
            outputf.write("\n")
        except ReadingError as e:
            outputf.write(l)





