#!/usr/bin/python
import sys
import os
import datetime

__version__='1.0b'

if sys.version_info.major != 2:
    print "This program requires python 2"
    sys.exit()
if sys.version_info.minor < 5:
    print "This program requires python 2.5 or more recent"
    sys.exit()

try:
    import numpy
except ImportError:
    print "This program requires the python numpy module is installed"
    sys.exit()

import getopt
import re
import csv

from LINZ.DeformationModel.Time import Time
from LINZ.DeformationModel import Model
from LINZ.DeformationModel.Error import ModelDefinitionError, OutOfRangeError, UndefinedValueError
from LINZ.Geodetic.ellipsoid import grs80

syntax='''
CalcDeformation.py: program to calculate deformation at a specific time and place using
the LINZ deformation model.

Syntax:
    python CalcDeformation.py [options] input_file output_file
    python CalcDeformation.py [options] -a longitude latitude

Options are:
  -d date           The date at which to calculate the deformation 
    --date=..       (default current date), or ":col_name"
  -b date           The reference date at which to calculate deformation.  The 
    --base_date=..   default is to calculate relative to the reference coordinates
  -u                Update the coordinates (default action is to calculate the 
    --update        displacement components) 
  -c lon:lat:h      Column names for the longitude, latitude, and height columns           
    --columns=      in input_file.  Default is lon, lat, hgt (optional)
  -f format         Format of the input and output files - format can be one of 
    --format=       csv (excel compatible csv), tab (tab delimited), 
                    w (whitespace delimited).
  -x de:dn:du       Displacement components to calculate, any of de, dn, du,
  --calculate=      eh, ev, separated by colon characters (default is de,dn,du)
  -g grid           Use a grid for input rather than an input file. The grid is
    --grid=         entered as "min_lon:min_lat:max_lon:max_lat:nlon:nlat"
  -a                Evaluate at longitude, latitude specified as arguments, rather
    --atpoint       using input/output files
  -m dir            Model base directory (default ../model)
    --model-dir=
  -v version        Version of model to calculate (default latest version)
    --version=      
  -r version        Calculate change relative to previous (base) version
    --baseversion=
  -p                Calculate reverse patch corrections 
    --patch 
  -k                Check that the model is correctly formatted - do not do any 
    --check         calculations
  -o submodel       Only calculate for specified submodel (ndm or patch directory name)
    --only=
  -l                Just list details of the model (input and output file 
    --list          are ignored)
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
base_version=None
subtract=False
reverse_patch=False
update=False
columns="lon lat hgt".split()
usercolumns=False
date_column=None
format="c"
date=None
base_date=None
listonly=False
check=False
griddef=None
inputfile=None
outputfile=None
quiet=False
atpoint=False
ptlon=None
ptlat=None
calcfields='de dn du eh ev'.split()
calculate=[0,1,2]

ell_a=6378137.0
ell_rf=298.257222101
ell_b=ell_a*(1-1.0/ell_rf)
ell_a2=ell_a*ell_a
ell_b2=ell_b*ell_b

if len(sys.argv) < 2:
    help()

optlist=None
args=None
try:
    optlist, args = getopt.getopt( sys.argv[1:], 'hd:b:uc:f:e:g:m:v:po:kqlxas', 
         ['help', 'date=', 'base_date=', 'apply','subtract','columns=','format=','elements=',
          'grid=','model-dir=','version=','patch','check',
          'only=','list','quiet','cache=','logging','atpoint'])
except getopt.GetoptError:
    print str(sys.exc_info()[1])
    sys.exit()

nargs = 2
maxargs = 2
usecache=True
clearcache=False
submodel=None
for o,v in optlist:
    if o in ('-h','--help'):
        help()
    if o in ('-l','--list'):
        listonly = True
        nargs=0
        maxargs=0
    elif o in ('-k','--check'):
        check = True
        nargs=0
        maxargs=0
    elif o in ('-d','--date'):
        if v.startswith(':'):
            date_column=v[1:]
        else:
            try:
                date=Time.Parse(v) 
            except:
                print "Invalid date "+v+" requested, must be formatted YYYY-MM-DD"
                sys.exit()
    elif o in ('-b','--base_date'):
       try:
            base_date=Time.Parse(v) 
       except:
            print "Invalid base date "+v+" requested, must be formatted YYYY-MM-DD"
            sys.exit()
    elif o in ('-a','--apply'):
        update=True
    elif o in ('-s','--subtract'):
        update=True
        subtract=True
    elif o in ('-c','--columns'):
        usercolumns=True
        columns=v.split(':')
        if len(columns) not in (2,3,4):
            print "Invalid columns specified - must be 2 or 3 colon separated column names"
            sys.exit()
    elif o in ('-f','--format'):
        v = v.lower()
        if v in ('csv','tab','whitespace','c','t','w'):
            format=v[:1]
        else:
            print "Invalid format specified, must be one of csv, tab, or whitespace"
            sys.exit()
    elif o in ('-e','--elements'):
        cols = v.lower().split(':')
        for c in cols:
            if c not in calcfields:
                print "Invalid calculated value "+c+" requested, must be one of "+' '.join(calcfields)
                sys.exit()
        calculate = [i for i,c in enumerate(calcfields) if c in cols]
    elif o in ('-g','--grid'):
        griddef=v
        nargs=1
        maxargs=1
    elif o in ('-x','--atpoint'):
        atpoint=True
        nargs=2
        maxargs=3
    elif o in ('-m','--model-dir'):
        modeldir = v
    elif o in ('-v','--version'):
        m=re.match(r'^(\d{8})?(?:\-(\d{8}))?$',v)
        if not v or not m:
            print "Invalid model version "+v+" selected"
            sys.exit()
        if m.group(1):
            version=m.group(1)
            if m.group(2):
                base_version=m.group(2)
        else:
            version=m.group(2)
            subtract=True
    elif o in ('-p','--patch'):
        reverse_patch = True
    elif o in ('-q','--quiet'):
        quiet = True
    elif o in ('-o', '--only'):
        submodel=v
    elif o in ('--cache'):
        if v in ('use','clear','ignore','reset'):
            usecache = v in ('use','reset')
            clearcache = v in ('clear','reset')
        else:
            print "Invalid cache option - must be one of use, clear, reset, ignore"
            sys.exit()
    elif o in ('--logging'):
        import logging
        logging.basicConfig(level=logging.INFO)
    else:
        print "Invalid parameter "+o+" specified"

if len(args) > maxargs:
    print "Too many arguments specified: " + " ".join(args[nargs:])
    sys.exit()
elif len(args) < nargs:
    if atpoint:
        print "Require longitude and latitude coordinate"
    elif nargs - len(args) == 2:
        print "Require input and output filename arguments"
    else:
        print "Require output filename argument"
    sys.exit()

if atpoint:
    try:
        ptlon=float(args[0])
        ptlat=float(args[1])
        pthgt=float(args[2]) if len(args)==3 else 0.0
    except:
        print "Invalid longitude/latitude "+args[0]+" "+args[1]
        sys.exit()
else:
    if nargs == 2:
        inputfile=args[0]
    if nargs > 0:
        outputfile=args[-1]

if not modeldir:
    modeldir=os.environ.get('NZGD2000_DEF_MODEL')

if not modeldir:
    from os.path import dirname, abspath, join, isdir, exists
    modeldir = join(dirname(dirname(abspath(__file__))),'model')
    modelcsv=join(modeldir,'model.csv')
    if not isdir(modeldir) or not exists(modelcsv):
        modeldir='model'

# Load the model, print its description if listonly is requested

model=None

# Use a loop to make exiting easy...
try:
    for loop in [1]:
        try:
            model = Model.Model(modeldir,loadAll=check,
                                useCache=usecache,clearCache=clearcache,loadSubmodel=submodel )
        except ModelDefinitionError:
            print "Error loading model:"
            print str(sys.exc_info()[1])
            break
        if check:
            print "The deformation model is correctly formatted"
            break

        # Set the model version

        if version == None:
            version = model.currentVersion()

        if reverse_patch:
            if base_version == None:
                base_version = model.versions()[0]
            model.setVersion( base_version, version )

            if date == None and date_column==None:
                date = model.datumEpoch()
            else:
                if not quiet:
                    print "Using a date or date column with a patch option - are you sure?"
        else:
            if date == None and date_column == None:
                date = Time.Now()
            model.setVersion( version, base_version )

        if listonly:
            print model.description()
            break

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
                def readf():
                    lat = min_lat-dlat
                    for ilat in range(nlat):
                        lat += dlat
                        lon = min_lon-dlon
                        for ilon in range(nlon):
                            lon += dlon
                            yield [str(lon),str(lat)]
                reader=readf
            except:
                print "Invalid grid definition",griddef
                break
            colnos=[0,1]
            headers=columns[0:2]
            
        else:
            try:
                instream = sys.stdin if inputfile=='-' else open(inputfile,"rb")
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
            columnsvalid=True
            if len(columns) > 3:
                date_column=date_column or columns[3]
                if date_column != columns[3]:
                    print "Inconsistent names specified for date column"
                    break
                columns=columns[:3]
                if columns[2] == '':
                    columns=columns[:2]
            for c in columns:
                if c in headers:
                    colnos.append(headers.index(c))
                elif c=='hgt' and not usercolumns:
                    break
                else:
                    print "Column",c,"missing in",inputfile
                    columnsvalid=False
            if not columnsvalid:
                break
            if date_column:
                if date_column in headers:
                    date_colno = headers.index(date_column)
                else:
                    print "Column",date_column,"missing in",inputfile
                    break

        # Create the output file

        if not quiet:
            action = "Updating with" if update else "Calculating"
            value = "patch correction" if reverse_patch else "deformation"
            vsnopt = "between versions "+base_version+" and "+version if base_version else "for version "+version
            datopt = "the date in column "+date_column if date_column else str(date)
            if base_date:
                datopt = "between "+str(base_date)+" and "+datopt
            else:
                datopt = "at "+datopt
            print "Deformation model "+model.name()
            print "for datum "+model.datumName()
            print action + " " + value + " " + vsnopt + " " + datopt

        if atpoint:
                defm = model.calcDeformation(ptlon,ptlat,date,base_date)
                if subtract:
                    for i in range(3):
                        defm[i]=-defm[i]
                if update:
                    dedln,dndlt=grs80.metres_per_degree(ptlon,ptlat)
                    ptlon += defm[0]/dedln
                    ptlat += defm[1]/dndlt
                    pthgt += defm[2]
                    print "{0:.8f} {1:.8f} {2:.4f}".format(ptlon,ptlat,pthgt)
                elif quiet:
                    print "{0:.4f} {1:.4f} {2:.4f}".format(defm[0],defm[1],defm[2])
                else:
                    print "Deformation at {0:.6f} {1:.6f}: {2:.4f} {3:.4f} {4:.4f}".format(
                        ptlon,ptlat,defm[0],defm[1],defm[2])
                break

        try:
            outstream = sys.stdout if outputfile=='-' else open(outputfile,"wb")
        except:
            print "Cannot open output file",outputfile
            break

        if not update:
            for c in calculate:
                headers.append(calcfields[c])

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

        latcalc=None
        dedln=None
        dndlt=None
        nerror=0
        nrngerr=0
        nmissing=0
        ncalc=0
        nrec=0

        for data in reader():
            nrec+=1
            if len(data) < ncols:
                if not quiet:
                    print "Skipping record",nrec,"as too few columns"
                    continue
            else:
                data=data[:ncols]
            try:
                lon = float(data[colnos[0]])
                lat = float(data[colnos[1]])
                if date_colno != None:
                    date = data[date_colno]
                defm = model.calcDeformation(lon,lat,date,base_date)
                if subtract:
                    for i in range(3):
                        defm[i] = -defm[i]
                if update:
                    dedln,dndlt=grs80.metres_per_degree(lon,lat)
                    lon += defm[0]/dedln
                    lat += defm[1]/dndlt
                    data[colnos[0]]="%.8lf"%(lon,)
                    data[colnos[1]]="%.8lf"%(lat,)
                    if len(colnos) > 2:
                        hgt = float(data[colnos[2]])
                        hgt += defm[2]
                        data[colnos[2]] = "%.3lf"%(hgt,)
                else:
                    for c in calculate:
                        data.append("%.4lf"%defm[c])
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
            print ncalc,"deformation values calculated"
            if nrngerr > 0:
                print nrngerr,"points were outside the valid range of the model"
            if nmissing > 0:
                print nmissing,"deformation values were undefined in the model"

except:
    print str(sys.exc_info()[1])


if model:
    model.close()
