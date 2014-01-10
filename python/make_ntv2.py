
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


parser=argparse.ArgumentParser('Convert patch to NTV2 format')
parser.add_argument('-b','--base-directory',default='.',help='Directory in which patch CSV files are located')
parser.add_argument('-v','--version',default='20130801',help='Distortion grid version')
parser.add_argument('-c','--created',default=date.today().strftime("%Y%m%d"),help='Created date')
parser.add_argument('-A','--australian',action='store_true',help='Build australian format binary')
parser.add_argument('-B','--big-endian',action='store_true',help='File is big endian (default little endian)')
parser.add_argument('-k','--keep-files',action='store_true',help='Keep working files')
parser.add_argument('--ogr2ogr-bug-workaround',action='store_true',help='Workaround for ogr2ogr nested grid bug (#177) - carefully ordered non-nested grids!')
parser.add_argument('ntv2_filename',default='grid.asc',help='Base name of the NTv2 grid file')

args=parser.parse_args()

basedir=args.base_directory
ntfile=args.ntv2_filename
version=args.version
created=args.created

if not os.path.isdir(basedir):
    raise RuntimeError('Base directory '+basedir+' is not a directory')

print "Creating NTv2 files as",ntfile

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
    # Also ensure that the 
    # This will cumulate as we are processing files from top to bottom
    if f in parent:
        lgfp=llgrid[parent[f]]
        lgfn=re.sub(r'\.csv$','_sum.csv',lgf)
        subprocess.call([
            'gridtool',
            'read','csv',
            basedir+'/'+lgf,
            'alignto','csv',
            basedir+'/'+lgfp,
            'trimto','csv',
            basedir+'/'+lgfp,
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

if args.ogr2ogr_bug_workaround:
    print "Applying ogr2ogr bug (proj issue #177) bug fix"
    filenames.reverse()
    parent=dict()

endian = '>' if args.big_endian else '<'
sformat=endian+'8s8s'
iformat=endian+'8si' if args.australian else endian+'8si4x'
dformat=endian+'8sd'
gsformat=endian+'ffff'

created=(created+'       ')[:8]
version=(version+'       ')[:8]

with open(ntfile+'.asc','w') as nt,open(ntfile+'.gsb','wb') as nb:
    nt.write("{0:8s}{1:3d}\n".format('NUM_OREC',11))
    nt.write("{0:8s}{1:3d}\n".format('NUM_SREC',11))
    nt.write("{0:8s}{1:3d}\n".format('NUM_FILE',len(filenames)))
    nt.write("{0:8s}{1:8s}\n".format('GS_TYPE','SECONDS'))
    nt.write("{0:8s}{1:8s}\n".format('VERSION',version))
    nt.write("{0:8s}{1:8s}\n".format('SYSTEM_F','GRS80'))
    nt.write("{0:8s}{1:8s}\n".format('SYSTEM_T','GRS80'))
    nt.write("{0:8s}{1:12.3f}\n".format('MAJOR_F',grs80.a))
    nt.write("{0:8s}{1:12.3f}\n".format('MINOR_F',grs80.b))
    nt.write("{0:8s}{1:12.3f}\n".format('MAJOR_T',grs80.a))
    nt.write("{0:8s}{1:12.3f}\n".format('MINOR_T',grs80.b))

    nb.write(struct.pack(iformat,'NUM_OREC',11))
    nb.write(struct.pack(iformat,'NUM_SREC',11))
    nb.write(struct.pack(iformat,'NUM_FILE',len(filenames)))
    nb.write(struct.pack(sformat,'GS_TYPE ','SECONDS '))
    nb.write(struct.pack(sformat,'VERSION ',version))
    nb.write(struct.pack(sformat,'SYSTEM_F','GRS80   '))
    nb.write(struct.pack(sformat,'SYSTEM_T','GRS80   '))
    nb.write(struct.pack(dformat,'MAJOR_F ',grs80.a))
    nb.write(struct.pack(dformat,'MINOR_F ',grs80.b))
    nb.write(struct.pack(dformat,'MAJOR_T ',grs80.a))
    nb.write(struct.pack(dformat,'MINOR_T ',grs80.b))

    ntgridname={}
    for i,f in enumerate(filenames):
        name="GRID{0:02d}  ".format(i)
        ntgridname[f]=name
        pname=ntgridname[parent[f]] if f in parent else 'NONE    '

        points=[]
        with open(basedir+'/'+llgrid[f]) as lgf:
            lgf.readline()
            nlon=0
            for line in lgf:
                if ',' not in line:
                    break
                lon,lat,dln,dlt=[float(f)*3600 for f in line.split(',')][:4]
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
            if nlat*nlon != len(points):
                raise RuntimeError("Error reading grid "+llgrid[f])
            loninc=(lon1-lon0)/(nlon-1)
            latinc=(lat1-lat0)/(nlat-1)
            nt.write("{0:8s}{1:8s}\n".format('SUB_NAME',name))
            nt.write("{0:8s}{1:8s}\n".format('PARENT',pname))
            nt.write("{0:8s}{1:8s}\n".format('CREATED',created))
            nt.write("{0:8s}{1:8s}\n".format('UPDATED',created))
            nt.write("{0:8s}{1:15.6f}\n".format('S_LAT',lat0))
            nt.write("{0:8s}{1:15.6f}\n".format('N_LAT',lat1))
            nt.write("{0:8s}{1:15.6f}\n".format('E_LONG',-lon1))
            nt.write("{0:8s}{1:15.6f}\n".format('W_LONG',-lon0))
            nt.write("{0:8s}{1:15.6f}\n".format('LAT_INC',latinc))
            nt.write("{0:8s}{1:15.6f}\n".format('LONG_INC',loninc))
            nt.write("{0:8s}{1:6d}\n".format('GS_COUNT',len(points)))

            nb.write(struct.pack(sformat,'SUB_NAME',name))
            nb.write(struct.pack(sformat,'PARENT  ',pname))
            nb.write(struct.pack(sformat,'CREATED ',created))
            nb.write(struct.pack(sformat,'UPDATED ',created))
            nb.write(struct.pack(dformat,'S_LAT   ',lat0))
            nb.write(struct.pack(dformat,'N_LAT   ',lat1))
            nb.write(struct.pack(dformat,'E_LONG  ',-lon1))
            nb.write(struct.pack(dformat,'W_LONG  ',-lon0))
            nb.write(struct.pack(dformat,'LAT_INC ',latinc))
            nb.write(struct.pack(dformat,'LONG_INC',loninc))
            nb.write(struct.pack(iformat,'GS_COUNT',len(points)))

            for ilat in range(nlat):
                ilatn=ilat*nlon
                for ilon in reversed(range(nlon)):
                    dln,dlt = points[ilatn+ilon]
                    nt.write("{0:10.6f}{1:10.6f}{2:10.6f}{3:10.6f}\n".format(dlt,-dln,-1.0,-1.0))
                    nb.write(struct.pack(gsformat,dlt,-dln,-1.0,-1.0))


if not args.keep_files:
    for f in working_files:
        os.unlink(basedir+'/'+f)

