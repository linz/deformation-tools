import sys
import math

infile=sys.argv[1]
cols=(int(x) for x in sys.argv[2:])

with open(infile) as inf:
    headers=inf.readline().split()
    for l in inf:
        values=l.split()
        for i in cols:
            v1,v2=(float(x) for x in values[i:i+2])
            mag=math.hypot(v1,v2)
            angle=math.degrees(math.atan2(v2,v1))
            print headers[i],mag,angle


