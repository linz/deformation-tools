#!/usr/bin/python
import sys
import getopt
import csv

from LINZ.DeformationModel.Time import Time
from LINZ.DeformationModel import Model
from LINZ.DeformationModel.Error import ModelDefinitionError, OutOfRangeError, UndefinedValueError
from LINZ.Geodetic import ITRF_NZGD2000

syntax='''
ITRF_NZGD2000.py: program to calculate deformation at a specific time and place using
the LINZ deformation model.

Syntax:
    python ITRF_NZGD2000.py [options] input_file output_file
    python ITRF_NZGD2000.py [options] -a longitude latitude

Options are:
  -d date           The date at which to calculate the transformation 
    --date=..       (default current date), or ":col_name"

  -i itrf           The ITRF reference frame to transform from/to. Eg 2008,
    --itrf=         ITRF2008, 97, etc. (default ITRF2008)
  -r                Reverse transformation (ie NZGD2000->ITRF)
    --reverse
  -c lon:lat:h      Column names for the longitude, latitude, and height columns           
    --columns=      in input_file.  Default is lon, lat, hgt (optional)
  -f format         Format of the input and output files - format can be one of 
    --format=       csv (excel compatible csv), tab (tab delimited), 
                    w (whitespace delimited).
  -g grid           Use a grid for input rather than an input file. The grid is
    --grid=         entered as "min_lon:min_lat:max_lon:max_lat:nlon:nlat"
  -x                Evaluate at longitude, latitude specified as arguments, rather
    --atpoint       using input/output files
  -m dir            Model base directory (default ../model)
    --model-dir=
  -v version        Version of model to calculate (default latest version)
    --version=      
  -q                Suppress optional output
    --quiet
  --cache=..        Cache options, ignore, clear, reset, use (default is use)
  --logging         Enable trace logging 
'''

def help():
    print syntax
    sys.exit()


modeldir=None
version=None
reverse=False
columns="lon lat hgt".split()
date_column=None
format="c"
date=None
griddef=None
inputfile=None
outputfile=None
quiet=False
atpoint=False
ptlon=None
ptlat=None
pthgt=None
itrf='ITRF2008'

if len(sys.argv) < 2:
    help()

optlist=None
args=None
try:
    optlist, args = getopt.getopt( sys.argv[1:], 'hd:c:f:g:i:m:v:rqlx', 
         ['help', 'date=', 'columns=','format=',
          'grid=','itrf=', 'model-dir=','version=',
          'quiet','cache=','reverse','logging','atpoint'])
except getopt.GetoptError:
    print str(sys.exc_info()[1])
    sys.exit()

nargs = 2
nargsmax=2
usecache=True
clearcache=False
for o,v in optlist:
    if o in ('-h','--help'):
        help()
    elif o in ('-d','--date'):
        if v.startswith(':'):
            date_column=v[1:]
        else:
            try:
                date=Time.Parse(v) 
            except:
                print "Invalid date "+v+" requested, must be formatted YYYY-MM-DD"
                sys.exit()
    elif o in ('-c','--columns'):
        columns=v.split(':')
        if len(columns) != 3:
            print "Invalid columns specified - must be 3 colon separated column names"
            sys.exit()
    elif o in ('-f','--format'):
        v = v.lower()
        if v in ('csv','tab','whitespace','c','t','w'):
            format=v[:1]
        else:
            print "Invalid format specified, must be one of csv, tab, or whitespace"
            sys.exit()
    elif o in ('-g','--grid'):
        griddef=v
        nargs=1
        nargsmax=1
    elif o in ('-x','--atpoint'):
        atpoint=True
        nargs=2
        nargsmax=3
    elif o in ('-i','--itrf'):
        itrf=v
    elif o in ('-m','--model-dir'):
        modeldir = v
    elif o in ('-v','--version'):
        version = v
    elif o in ('-r','--reverse'):
        reverse=True
    elif o in ('-q','--quiet'):
        quiet = True
    elif o in ('--cache'):
        if v in ('use','clear','ignore','reset'):
            usecache = v in ('use','reset')
            clearcache = v in ('clear','reset')
        else:
            print "Invalid cache option - must be one of use, clear, reset, ignore"
            sys.exit()
    elif o in ('--logging'):
        logging.basicConfig(level=logging.INFO)
    else:
        print "Invalid parameter "+o+" specified"

if len(args) > nargsmax:
    print "Too many arguments specified: " + " ".join(args[nargs:])
    sys.exit()
elif len(args) < nargs:
    if atpoint:
        print "Require longitude, latitude and optional height coordinate"
    elif nargs - len(args) == 2:
        print "Require input and output filename arguments"
    else:
        print "Require output filename argument"
    sys.exit()

if atpoint:
    try:
        ptlon=float(args[0])
        ptlat=float(args[1])
        pthgt=0.0
        if len(args) > 2:
            pthgt=float(args[2])
    except:
        print "Invalid longitude/latitude "+args[0]+" "+args[1]
        sys.exit()
else:
    if nargs == 2:
        inputfile=args[0]
    if nargs > 0:
        outputfile=args[-1]

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
            toNZGD2000=not reverse,
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


    if date is None and date_column is None:
        date = Time.Now()

    # Determine the source for input

    reader = None
    headers = None
    colnos=None
    date_colno = None

    ncols = 2
    dialect = csv.excel_tab if format =='t' else csv.excel;
    if atpoint:
        pass

    elif griddef:
        # Grid format
        try:
            parts=griddef.split(':')
            if len(parts) != 6:
                raise ValueError('')
            min_lon=float(parts[0])
            min_lat=float(parts[1])
            max_lon=float(parts[2])
            max_lat=float(parts[3])
            nlon=int(parts[4])
            nlat=int(parts[5])
            if max_lon<=min_lon or max_lat<=min_lat or nlon<2 or nlat < 2:
                raise ValueError('')
            dlon = (max_lon-min_lon)/(nlon-1)
            dlat = (max_lat-min_lat)/(nlat-1)
            datestr=date.strftime()
            def readf():
                lat = min_lat-dlat
                for ilat in range(nlat):
                    lat += dlat
                    lon = min_lon-dlon
                    for ilon in range(nlon):
                        lon += dlon
                        yield [str(lon),str(lat),'0.0',datestr,str(lon),str(lat),'0.0']
            reader=readf
        except:
            raise
            print "Invalid grid definition",griddef
            break
        ncols=7
        colnos=[4,5,6]
        header_itrf=('lon_'+itrf.lower(),'lat_'+itrf.lower(),'hgt_'+itrf.lower())
        header_nzgd2000=('lon_nzgd2000','lat_nzgd2000','hgt_nzgd2000')
        headers=[]
        if reverse:
            headers.extend(header_nzgd2000)
            headers.append('date')
            headers.extend(header_itrf)
        else:
            headers.extend(header_itrf)
            headers.append('date')
            headers.extend(header_nzgd2000)
        
    else:
        try:
            instream = open(inputfile,"rb")
        except:
            print "Cannot open input file "+inputfile
            break
        # Whitespace
        if format == 'w':
            headers=instream.readline().strip().split()
            def readf():
                for line in instream:
                    yield line.strip().split()
            reader = readf
        # CSV format
        else:
            csvrdr = csv.reader(instream,dialect=dialect)
            headers=csvrdr.next()
            def readf():
                for r in csvrdr:
                    yield r
            reader = readf
        ncols = len(headers)
        colnos=[]
        for c in columns:
            if c in headers:
                colnos.append(headers.index(c))
            else:
                break
        if len(colnos) < 3:
            print "Column",c,"missing in",inputfile
            break
        if date_column:
            if date_column in headers:
                date_colno = headers.index(date_column)
            else:
                print "Column",date_column,"missing in",inputfile
                break

            date_colno = colno

    # Create the output file

    if not quiet:
        action = "Converting "
        action = action+"NZGD2000 to "+itrf if reverse else action+itrf+" to NZGD2000" 
        datopt = "the date in column "+date_column if date_column else str(date)
        action = action + ' at '+datopt
        print action
        print "Deformation model "+transform.model.name() + " version "+transform.version
        print "for "+transform.model.datumName()

    if atpoint:
            llh = transform(ptlon,ptlat,pthgt, date )
            print "{0:.9f} {1:.9f} {2:.4f}".format(*llh)
            break

    try:
        outstream = open(outputfile,"wb")
    except:
        print "Cannot open output file",outputfile
        break

    writefunc = None
    if format=='w':
        def writef(cols):
            outstream.write(' '.join(cols))
            outstream.write("\n")
        writefunc = writef
    else:
        csvwrt = csv.writer(outstream,dialect=dialect)
        writefunc=csvwrt.writerow

    writefunc(headers)

    ncalc=0
    nrngerr=0
    nmissing=0
    ncolnos=len(colnos)
    for data in reader():
        if len(data) < ncols:
            continue
        else:
            data=data[:ncols]
        try:
            lon = float(data[colnos[0]])
            lat = float(data[colnos[1]])
            hgt = float(data[colnos[2]]) if ncolnos > 2 else 0.0
            if date_colno is not None:
                date = data[date_colno]
            llh = transform(lon,lat,hgt, date)
            data[colnos[0]]="{0:.9f}".format(llh[0])
            data[colnos[1]]="{0:.9f}".format(llh[1])
            if ncolnos > 2:
                data[colnos[2]]="{0:.4f}".format(llh[2])
            writefunc(data)
            ncalc += 1
        except OutOfRangeError:
            nrngerr += 1
        except UndefinedValueError:
            nmissing += 1
        except:
            raise
            print str(sys.exc_info()[1])
            nerror += 1

    if not quiet:
        print ncalc," coordinates values converted"
        if nrngerr > 0:
            print nrngerr,"points were outside the valid range of the model"
        if nmissing > 0:
            print nmissing,"deformation values were undefined in the model"
