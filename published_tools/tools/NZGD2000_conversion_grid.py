#!/usr/bin/python
import sys
import argparse
import numpy as np

from LINZ.DeformationModel.Time import Time
from LINZ.DeformationModel import Model
from LINZ.DeformationModel.Error import ModelDefinitionError, OutOfRangeError, UndefinedValueError
from LINZ.Geodetic import ITRF_NZGD2000
from LINZ.Geodetic.ellipsoid import grs80

parser=argparse.ArgumentParser(description='''
Calculate the corrections to convert an ITRF coordinate to or from 
an NZGD2000 coordinate on a longitude/latitude grid.
The corrections are generated in a comma separated (CSV) file.''')
parser.add_argument('min_lon',type=float,help="Minimum longitude")
parser.add_argument('min_lat',type=float,help="Minimum latitude")
parser.add_argument('max_lon',type=float,help="Maximum longitude")
parser.add_argument('max_lat',type=float,help="Maximum latitude")
parser.add_argument('nlon',type=float,help="Number of longitude values")
parser.add_argument('nlat',type=float,help="Number of latitude values")
parser.add_argument('output_file',type=str,help="Output file name")
parser.add_argument('-d','--date',type=str,default='now',
                    help="Date of transformation (default current date)")
parser.add_argument('-i','--itrf',type=str, default='ITRF2008',
                    help="ITRF reference frame version (default ITRF2008)")
parser.add_argument('-m','--model-dir',type=str,default='../model',
                    help="Deformation model base directory (default ../model)")
parser.add_argument('-v','--version',type=str,
                    help="Version of deformation model to use (default current version)")
parser.add_argument('-t','--type',choices=('llh','enu','xyz'),default='llh',
                    help="Type of correction to calculate")
parser.add_argument('-s', '--size', action='store_true',
                    help="Use grid cell size (dlon,dlat) instead of number of values (nlon,nlat)")
parser.add_argument('-r','--reverse',action='store_true',
                    help="Generate grid to convert NZGD2000 to ITRF")
parser.add_argument('-o','--order',choices='es en ws wn se sw ne nw'.split(),
                    default='es',
                    help="Grid order (start corner and first axis to increment)")
parser.add_argument('-q','--quiet',action='store_true',
                    help="Suppress output")
parser.add_argument('--cache',choices=('ignore','clear','reset','use'),
                    default='use',
                    help='Deformation model cache option (requires pytables)')
parser.add_argument('--logging',action='store_true',help="Enable trace logging")

args=parser.parse_args()


modeldir=args.model_dir
version=args.version
date=args.date
try:
    date=Time.Parse(args.date) 
except:
    print "Invalid date "+v+" requested, must be formatted YYYY-MM-DD"
    sys.exit()
outputfile=args.output_file
itrf=args.itrf
increment=args.size
order=args.order
reverse=args.reverse
corrtype=args.type
quiet=args.quiet
usecache = args.cache in ('use','reset')
clearcache = args.cache in ('clear','reset')

if args.logging:
    logging.basicConfig(level=logging.INFO)

if not modeldir:
    from os.path import dirname, abspath, join
    modeldir = join(dirname(dirname(abspath(__file__))),'model')

# Use a loop to make exiting easy...

for loop in [1]:
    # Setup the transformation
    transform=None
    try:
        transform=ITRF_NZGD2000.Transformation(
            itrf,
            toNZGD2000=not args.reverse,
            modeldir=modeldir,
            version=version,
            usecache=usecache,
            clearcache=clearcache )
    except ModelDefinitionError:
        print "Error loading model:"
        print str(sys.exc_info()[1])
        break
    except RuntimeError:
        print str(sys.exc_info()[1])
        break

    # Determine the source for input

    coords=[]
    try:
        min_lon=args.min_lon
        min_lat=args.min_lat
        max_lon=args.max_lon
        max_lat=args.max_lat
        if max_lon<=min_lon or max_lat<=min_lat:
            raise ValueError('Minimum latitude or longitude larger than maximum')
        lonval=None
        latval=None
        if increment:
            dlon=args.nlon
            dlat=args.nlat
            if dlon <= 0 or dlat <= 0:
                raise ValueError('')
            lonval=np.arange(min_lon,max_lon+dlon*0.99,dlon)
            latval=np.arange(min_lat,max_lat+dlat*0.99,dlat)
        else:
            nlon=int(args.nlon)
            nlat=int(args.nlat)
            if nlon<2 or nlat < 2:
                raise ValueError('Must be at least two longitude and latitude values')
            lonval=np.linspace(min_lon,max_lon,nlon)
            latval=np.linspace(min_lat,max_lat,nlat)

        if 'w' in order:
            lonval=lonval[::-1]
        if 'n' in order:
            latval=latval[::-1]
        if order[0] in 'ew':
            for lat in latval:
                for lon in lonval:
                    coords.append((lon,lat))
        else:
            for lon in lonval:
                for lat in latval:
                    coords.append((lon,lat))

    except ValueError as e:
        print "Invalid grid definition: "+e.message
        break

    # Create the output file

    if not quiet:
        if reverse:
            print "Calculating NZGD2000 to "+itrf+" corrections at "+str(date)
        else:
            print "Calculating "+itrf+" to NZGD2000 corrections at "+str(date)
        print "Deformation model "+transform.model.name() + " version "+transform.version

    try:
        outstream = open(outputfile,"wb")
    except:
        print "Cannot open output file",outputfile
        break


    if corrtype=="llh":
        outstream.write("lon,lat,dlon,dlat,dhgt\n")
    elif corrtype=="enu":
        outstream.write("lon,lat,de,dn,dh\n")
    elif corrtype=="xyz":
        outstream.write("lon,lat,dx,dy,dz\n")


    ncalc=0
    nrngerr=0
    nmissing=0
    hgt=0.0
    rvs=-1.0 if reverse else 1.0
    for lon,lat in coords:
        try:
            llh = transform(lon,lat,hgt,date)
            if corrtype=="llh":
                outstream.write("{0:.5f},{1:.5f},{2:.9f},{3:.9f},{4:.4f}\n".format(
                    lon,lat,rvs*(llh[0]-lon),rvs*(llh[1]-lat),rvs*llh[2]))
            elif corrtype=="enu":
                dedln,dndlt=grs80.metres_per_degree(lon,lat)
                outstream.write("{0:.5f},{1:.5f},{2:.4f},{3:.4f},{4:.4f}\n".format(
                    lon,lat,rvs*dedln*(llh[0]-lon),rvs*dndlt*(llh[1]-lat),rvs*llh[2]))
            elif corrtype=="xyz":
                xyz0=grs80.xyz(lon,lat,hgt)
                xyz1=grs80.xyz(llh[0],llh[1],llh[2])
                outstream.write("{0:.5f},{1:.5f},{2:.4f},{3:.4f},{4:.4f}\n".format(
                    lon,lat,
                    rvs*(xyz1[0]-xyz0[0]),
                    rvs*(xyz1[1]-xyz0[1]),
                    rvs*(xyz1[2]-xyz0[2])))
            ncalc += 1
        except OutOfRangeError:
            nrngerr += 1
        except UndefinedValueError:
            nmissing += 1
        except:
            raise
            print str(sys.exc_info()[1])
            nerror += 1

    outstream.close()

    if not quiet:
        print "{0} corrections calculated".format(ncalc)
    if nrngerr > 0:
        print "{0} points were outside the valid range of the model".format(nrngerr)
    if nmissing > 0:
        print "{0} deformation values were undefined in the model".format(nmissing)
