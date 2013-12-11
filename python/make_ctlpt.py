
import argparse
import os
import os.path
import re
from collections import namedtuple
import subprocess
from math import *
from datetime import date
import struct

class ellipsoid( object ):

    convergence=1.0e-10

    @staticmethod
    def _cossin( angle ):
        angle=radians(angle)
        return cos(angle),sin(angle)

    def __init__( self, a, rf ):
        '''
        Initiallize an ellipsoid based on semi major axis and inverse flattening
        '''
        self.a=float(a)
        self.rf=float(rf)
        self.b=a-a/rf if rf else a
        self.a2=a*a
        self.b2=self.b*self.b
        self.a2b2=self.a2-self.b2

    def metres_per_degree( self, lat, hgt=0 ):
        '''
        Calculate the number of metres per degree east and
        north
        '''
        clt,slt = ellipsoid._cossin(lat)
        bsac=hypot(self.b*slt,self.a*clt)
        p = self.a2*clt/bsac + hgt*clt
        dedln=radians(p)
        dndlt=radians((self.a2*self.b2)/(bsac*bsac*bsac)+hgt)
        return dedln,dndlt

grs80 = ellipsoid(6378137.0,298.257222101)

def make_ll_grid( f ):
    if not f.endswith(".csv"):
        raise RuntimeError( "File"+f+" is not a .csv file - ignored")
    with open(f) as inf:
        l=inf.readline().strip()
        if l != "lon,lat,de,dn,du":
            raise RuntimeError( "File"+f+"is not a shift csv file - incorrect header line")
        f1=f[:-4]+"_ll.csv"
        dedln=1.0
        dndlt=1.0
        glat=-99.0
        npt=0
        with open(f1,"w") as outf:
            # print "Converting file",f,"to",f1
            outf.write("lon,lat,dlon,dlat,du\n")
            for l in inf:
                data=[float(x) for x in l.strip().split(',')]
                if len(data) != 5:
                    continue
                lat=data[1]
                if glat != lat:
                    glat=lat
                    dedln,dndlt=grs80.metres_per_degree(data[1])
                data[2] /= dedln
                data[3] /= dndlt
                outf.write(','.join((str(x) for x in data)))
                outf.write('\n')
                npt += 1
            # print npt,"grid points processed"
    return f1


parser=argparse.ArgumentParser('Convert patch to ESRI control point file')
parser.add_argument('-b','--base-directory',default='.',help='Directory in which patch CSV files are located')
parser.add_argument('-k','--keep-files',action='store_true',help='Keep working files')
parser.add_argument('-o','--include-offsets',action='store_true',help='Add dln,dlt offset columns')
parser.add_argument('ctlpt_filename',default='ctlpts.dat',help='Control point file name')

args=parser.parse_args()

basedir=args.base_directory
cptfile=args.ctlpt_filename

if not os.path.isdir(basedir):
    raise RuntimeError('Base directory '+basedir+' is not a directory')

print "Creating control point file as",cptfile

gridfile=namedtuple('gridfile','filename parent children')
gridlevel=lambda x: int(re.search(r'L(\d+)\.csv$',x).group(1))

filenames = [x for x in os.listdir(basedir) if re.search(r'L\d+\.csv$',x)]
filenames.sort(key=gridlevel)

parent=dict()
llgrid=dict()
children={f:[] for f in filenames}
working_files=[]

for f in filenames:
    level=gridlevel(f)
    for pf in reversed(filenames):
        if gridlevel(pf) >= level:
            continue
        pfs = re.sub(r'\_S?L\d+\.csv','',pf)
        if f[:len(pfs)] == pfs:
            parent[f]=pf
            children[pf].append(f)
            # print pf,"=>",f
            break

for f in filenames:
    lgf=make_ll_grid(basedir+'/'+f)
    lgf=lgf[len(basedir)+1:]
    working_files.append(lgf)
    # If there is a parent, then add the parent offsets
    # This will cumulate as we are processing files from top to bottom
    if f in parent:
        lgfp=llgrid[parent[f]]
        lgfn=re.sub(r'\.csv$','_sum.csv',lgf)
        subprocess.call([
            'gridtool',
            'read','csv',
            basedir+'/'+lgf,
            'add','csv',
            basedir+'/'+lgfp,
            'write','csv',
            basedir+'/'+lgfn
        ])
        lgf=lgfn
        working_files.append(lgf)
    llgrid[f]=lgf

# ogr2ogr nested grid bug (proj issue #177) workaround - construct the grids started with deepest 
# level and ignoring nesting.  Looks like proj/ogr2ogr just takes the first grid that applies at 
# a point, so this should just pick the most detailed grid.

griddef={}

def ongrid(x,g0,dg):
    return fabs(fmod(fabs((x-g0))/dg+0.5,1.0)-0.5) < 0.001

cptformat=(
    "{0} {1:.8f} {2:.8f} {3:.8f} {4:.8f} {5:.5f} {6:.5f}\n"
    if args.include_offsets else
    "{0} {1:.8f} {2:.8f} {3:.8f} {4:.8f}\n"
)

with open(cptfile,'w') as cpt:
    gptid=0
    nskip=0
    for i,f in enumerate(filenames):
        ischild = f in parent
        if ischild:
            plon0,pdlon,plat0,pdlat=griddef[parent[f]]

        points=[]
        with open(basedir+'/'+llgrid[f]) as lgf:
            lgf.readline()
            nlon=0
            for line in lgf:
                if ',' not in line:
                    break
                lon,lat,dln,dlt=[float(fv) for fv in line.split(',')][:4]
                isdone=False
                if ischild:
                    isdone = ongrid(lon,plon0,pdlon) and ongrid(lat,plat0,pdlat)
                    if isdone:
                        nskip += 1 
                if not isdone:
                    gptid += 1
                    cpt.write(cptformat.format(
                        gptid,lon,lat,lon+dln,lat+dlt,dln*3600,dlt*3600))
                if len(points) == 0:
                    lat0=lat
                    lon0=lon
                else:
                    if nlon==0 and lon < lon1:
                        nlon=len(points)
                lat1=lat
                lon1=lon
                points.append([dln,dlt])
            nlat=int(len(points)/nlon)
            loninc=(lon1-lon0)/(nlon-1)
            latinc=(lat1-lat0)/(nlat-1)
            griddef[f]=[lon0,loninc,lat0,latinc]

print nskip,"duplicated points skipped"

if not args.keep_files:
    for f in working_files:
        os.unlink(basedir+'/'+f)

