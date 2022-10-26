#!/usr/bin/python
#
# Script to test the reverse patch binary blob against the
# published model and the source model files, and to test
# the snap binary against the published model
#


import sys
import os.path
sys.path.append(os.path.join(os.path.dirname(__file__),'..','tools','python'))
import os
import re
import csv
import math
import numpy as np
import numbers
import ellipsoid
import atexit
from random import uniform, seed
from subprocess import call
from datetime import datetime, timedelta
from collections import namedtuple
from subprocess import call, check_output
from ellipsoid import grs80
from distutils.spawn import find_executable as which
import euler

tooldir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
calc_okada_image=os.path.join(tooldir,'okada','calc_okada')
if not os.path.exists(calc_okada_image):
    raise RuntimeError('Cannot find calc_okada program at '+calc_okada_image)

ntestpergrid=20
tmpfiles=[]

def tempfilename():
    tfn=os.tempnam()
    if tfn not in tmpfiles:
        tmpfiles.append(tfn)
    return tfn

def remove_temp_files():
    for tfn in tmpfiles:
        try:
            os.remove(tfn)
        except:
            pass

atexit.register(remove_temp_files)

def load_area( polygon_file, area_type, qclog ):
    from shapely.wkt import loads
    try:
        with open(polygon_file) as laf:
            wkt=laf.read()
            land_area=loads(wkt)
        qclog.printlog( "Using {0} definition from {1}".format(area_type,polygon_file))
        return land_area
    except Exception as ex:
        raise RuntimeError("Cannot load {0} definition from {1}\n{2}"
                           .format(area_type,polygon_file,ex.message))

def test_points( ldm, include_areas=None, points_per_grid=10, version=None ):
    ranges=[]
    rangestr=set()
    for c in sorted(ldm.components(),key=lambda x: x.name):
        if version is not None and c.versionAdded != version:
            continue
        sm=c.spatialModel
        for m in sm.models():
            grange=[m.min_lon,m.min_lat,m.max_lon,m.max_lat]
            if str(grange) not in rangestr:
                rangestr.add(str(grange))
                ranges.append(grange)

    testpts=[]
    for r in ranges:
        xmin,ymin,xmax,ymax=r
        if xmax > 180:
            xmax=180
        for i in range(points_per_grid):
            tries=20
            while True:
                lon= round(uniform(xmin,xmax),5)
                lat= round(uniform(ymin,ymax),5)
                if include_areas:
                    from shapely.geometry import Point
                    pt=Point(lon,lat)
                    ok=True
                    for use_inside,area in include_areas:
                        inside=area.contains(pt)
                        if (inside and not use_inside) or (not inside and use_inside):
                            ok=False
                            break
                    if not ok:
                        tries -= 1
                        if tries > 0:
                            continue
                        break
                testpts.append(( lon, lat, round(uniform(0.0,1000.0),3)))
                break
    return np.array(testpts)

def load_test_points( tpfilename ):
    testpts=[]
    with open(tpfilename,'rb') as tpf:
        csvr=csv.reader(tpf)
        for r in csvr:
            try:
                lon=float(r[1])
                lat=float(r[2])
            except:
                continue
            try:
                hgt=float(r[3])
            except:
                hgt=0.0
            testpts.append((lon,lat,hgt))
    return np.array(testpts)

def points_file( pts, format="{0} {1} {2}", header=None ):
    pfn=tempfilename()
    with open(pfn,'w') as pf:
        if header is not None:
            pf.write(header)
            pf.write('\n')
        for i,p in enumerate(pts):
            pf.write(format.format(i+1,p[0],p[1],p[2]))
            pf.write('\n')
    return pfn

def read_denu( filename=None,source=None,cols=(2,5),header=False,expected=None,delim=None,fillfrom=None,factor=1.0,relative=None):
    denu=[]
    ncol=cols[1]
    hline=''
    dline=''
    ndata=cols[1]-cols[0]
    if filename:
        source=open(filename).readlines()
    start=1 if header else 0
    hline=source[0].strip() if source and header else ''
    if relative is not None:
        expected=len(relative)

    for ipt,l in enumerate(source[start:]):
        if not dline:
            dline=l
        parts=l.strip().split(delim)
        if len(parts) >= ncol:
            if fillfrom is not None:
                id=int(parts[0])
                while fillfrom < id:
                    denu.append([0.0]*ndata)
                    fillfrom += 1
                fillfrom=id+1
            ldata=[float(x)*factor for x in parts[cols[0]:cols[1]]]
            if relative is not None:
                if ipt >= len(relative):
                    raise RuntimeError('Too many points in output data')
                relpt=relative[ipt]
                dedln,dndlt=grs80.metres_per_degree(relpt[0],relpt[1])
                ldata[0]=(ldata[0]-relpt[0])*dedln
                ldata[1]=(ldata[1]-relpt[1])*dndlt
                if len(ldata) > 2:
                    ldata[2]=ldata[2]-relpt[2]
            while len(ldata) < 3:
                ldata.append(0.0)
            denu.append(ldata)
    if expected is not None and len(denu) != expected:
        if fillfrom is not None and len(denu) < expected:
            denu.extend([[0.0]*ndata]*(expected-len(denu)))
        else:
            raise RuntimeError('Wrong number of dislocations returned: {0} of {1}\n{2}{3}'.
                           format(len(denu),expected,hline,dline))
    return np.array(denu)


def date_as_year( dt ):
    dy0=datetime(dt.year,1,1)
    dy1=datetime(dt.year+1,1,1)
    return dt.year+float((dt-dy0).days)/float((dy1-dy0).days)

def parse_date( datestr ):
    datestr=datestr.strip()
    for format in ('%d %b %Y','%d %B %Y','%Y-%m-%d'):
        try:
            return datetime.strptime(datestr,format)
        except:
            pass
    raise ValueError("Cannot parse date \"{0}\"".format(datestr))

def fault_models( modeldir, points, calcmodel=True ):

    class FaultModel(namedtuple('FaultModel','file version ramp denu')):

        def factor(self, t=None, reverse=False):
            if isinstance(t,datetime):
                t=date_as_year(t)
            ramp=self.ramp
            factor=ramp[-1].factor
            if t is None or t >= ramp[-1].year:
                pass
            elif t < ramp[0].year:
                 factor=0.0
            else:
                for n,rp in enumerate(ramp):
                    if n < rp.year:
                        rn=ramp[n+1]
                        factor= (rn.factor*(rp.year-t)+rp.factor*(t-rn.year))/(rp.year-rn.year)
                        break
            if reverse:
                factor -= ramp[-1].factor
            return factor

        def calc_denu(self,t=None,t0=None,reverse=False,qclog=None):
            factor=self.factor(t,reverse)
            if t0 is not None:
                factor -= self.factor(t0,reverse)
            if qclog:
                qclog.printlog("Using {0} times {1:.4f}".format(os.path.basename(self.file),factor),
                           verbose=1)
            return self.denu*factor

    class TimeStep(namedtuple('TimeStep','time year factor')):

        @staticmethod
        def parse( tsstr ):
            tsstr=tsstr.strip()
            if tsstr == '':
                return None
            (datestr,valuestr)=tsstr.split()[:2]
            tsdate=parse_date(datestr)
            tsvalue=float(valuestr)
            tsyear=date_as_year(tsdate)
            return TimeStep(tsdate,tsyear,tsvalue)

    models = []

    for f in sorted(os.listdir(modeldir)):
        if f.endswith('.model'):
            mfile=os.path.join(modeldir,f)
            pfile=os.path.splitext(mfile)[0]+'.patch'
            timemodel=[]
            inmodel=False
            version=None
            with open(pfile) as rf:
                for l in rf:
                    if re.match(r'^\s*(\#|$)',l):
                        continue
                    cmdmatch=re.match(r'^\s*(\w+)\s*\:\s*(.*?)\s*$',l)
                    if not cmdmatch: 
                        if inmodel:
                            ts=TimeStep.parse(l)
                            if ts:
                                timemodel.append(ts)
                        continue
                    inmodel=False
                    command=cmdmatch.group(1).lower()
                    value=cmdmatch.group(2)
                    if command == 'date':
                        mdate=parse_date(value)
                        timemodel=[TimeStep(mdate,date_as_year(mdate),1.0)]
                    elif command == 'timemodel':
                        timemodel=[]
                        ts=TimeStep.parse(value)
                        if ts is not None:
                            timemodel.append(ts)
                        inmodel=True
                    elif command == 'version':
                        version=value
            denu=None
            if calcmodel:
                denu=calc_okada(mfile,points)
            models.append(FaultModel(mfile,version,timemodel,denu))

    return models


def test_times( models, version=None ):
    transitions=set()
    for m in models:
        if version is not None and m.version != version:
            continue
        d0=None
        for step in m.ramp:
                tdate=step.time
                if d0:
                    nd=(tdate-d0).days/3
                    d1=d0+timedelta(days=nd)
                    transitions.add(d1)
                    d1=d0+timedelta(days=2*nd)
                    transitions.add(d1)
                d0=tdate
                transitions.add(d0-timedelta(days=1))
        if d0:
            transitions.add(d0+timedelta(days=1))
    times=sorted(list(transitions))
    return times

def calc_okada( m, pts ):
    global calc_okada_image
    print("Calculating fault model file",m)
    tn1=points_file(pts,format="{1} {2}")
    tn2=tempfilename()
    call((calc_okada_image,'-f',m,tn1,tn2))
    os.remove(tn1)
    disloc=[]
    try:
        disloc=read_denu(tn2,header=True,expected=len(pts))
    finally:
        os.remove(tn2)
    return disloc

def find_patches( pdir, pts ):
    pfiles=[]
    disloc=[]
    for base,dirs,files in os.walk(pdir):
        for f in files:
            if f.lower().endswith('.bin'):
                pfiles.append(os.path.join(base,f))
    Patch=namedtuple('Patch','file denu')
    tn1=points_file(pts)
    tn2=tempfilename()
    try:
        for pf in pfiles:
            call(('runshift','-p',pf,tn1,tn2))
            try:
                denu=read_denu(tn2,cols=(3,6),expected=len(pts))
                yield Patch(pf,np.array(denu))
            finally:
                os.remove(tn2)
    finally:
        os.remove(tn1)

def calc_velocities( gnsvel, eulerdef, points, qclog ):
    with open(eulerdef) as f:
        eparams=f.readline().split()
        evals=[float(x) for x in eparams]
    pfile=tempfilename()
    with open(pfile,'w') as pf:
        pf.write('{0}\n'.format(len(points)))
        for i,p in enumerate(points):
            pf.write('{0} {1} {2}\n'.format(i+1,p[1],p[0]))
    ofile=tempfilename()
    call(('gns_velocity_linz',gnsvel,pfile,ofile,'0','0','0','0'))
    venu=np.zeros((len(points),3))
    nbad=0
    erot=euler.euler_rotation(grs80,evals[1],evals[0],evals[2])
    with open(ofile) as of:
        of.readline()
        for l in of:
            parts=l.split()
            if len(parts) < 5:
                continue
            try:
                id=int(parts[0])
                (lat,lon,ve,vn)=(float(x) for x in parts[1:5])
                de,dn=erot.velocity(lon,lat)
                venu[id-1]=[ve*0.001-de,vn*0.001-dn,0.0]
            except:
                nbad += 1
    qclog.printlog("Calculating velocities from ",gnsvel)
    qclog.printlog(" {0} points could not be calculated".format(nbad))
    return venu

class QcLog( object ):

    def __init__( self, qcdir, prefix='', verbose=0 ):
        if not os.path.isdir(qcdir):
            os.makedirs(qcdir)
        prefix=prefix or ''
        if prefix:
            prefix=re.sub(r'_?$','_',prefix)
        self.dir=qcdir
        self.prefix=prefix
        self.verbose=verbose if type(verbose) == int else 1 if verbose else 0
        qclog=open(self.qcfilename('qc.log'),'w')
        self.log=qclog
        self.log.write("Running deformation QC tests\n")
        self.log.write("Runtime: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'\n')

    def qcfilename( self, filename ):
        return os.path.join(self.dir,self.prefix+filename)

    def printlog( self, *messages, **options ):
        verbose=options.get('verbose',0)
        verbose=verbose if type(verbose) == int else 1 if verbose else 0
        if verbose <= self.verbose:
            self.log.write(' '.join(messages))
            self.log.write('\n')

    def start_test(self,message):
        self.log.write("\n")
        self.log.write(message)
        self.log.write("\n")

    def write_test_points(self,points):
        tpfilename=self.qcfilename('test_points.csv')
        csvw=csv.writer(open(tpfilename,'w'))
        csvw.writerow(('id','lon','lat'))
        for i,p in enumerate(points):
            csvw.writerow((str(i),str(p[0]),str(p[1])))

    def compare(self,
                description,
                points,
                datasets,
                csvfile=None,
                csvsummary=None,
                tolerance=0.001,
                skip=None,
                listall=False,
                date0=None,
                date=None,
                append=False
               ):

        if len(datasets) > 2:
            for i,d0 in enumerate(datasets):
                for d1 in datasets[i+1:]:
                    self.compare(
                        description,
                        points,
                        (d0,d1),
                        csvfile,
                        csvsummary,
                        tolerance,
                        skip,
                        listall,
                        date0,
                        date,
                        append
                    )
                    append=True
            return

        if isinstance(tolerance,numbers.Number):
            tolerances=[tolerance]
        else:
            tolerances=list(tolerance)
        tolerance=tolerances[0]

        self.start_test(description)
        self.log.write('Comparing:\n')
        for i,d in enumerate(datasets):
            self.log.write("    {0}: {1}\n".format(d[0],d[1]))

        # Find rows exceeding tolerance

        error=np.array([0.0]*len(points))
        d0,d1=datasets
        diff=d1[2]-d0[2]
        error=np.sqrt(np.sum(np.square(diff),axis=1))
        ds=np.sqrt(np.sum(np.square(diff[:,:2]),axis=1))
        dh=np.abs(diff[:,2])

        skiprows=[]
        if skip is not None:
            skiprows=np.nonzero(np.sum(np.square(d0[2]),axis=1) < skip*skip)[0]
            error[skiprows]=0.0

        badrows=np.nonzero(error>tolerance)[0]

        nbad=len(badrows)
        ntest=len(points)-len(skiprows)

        if nbad:
            self.log.write('**** {0} of {1} differences exceed tolerance {2}\n'.format(nbad,ntest,tolerance))
            for tol in tolerances[1:]:
                nfail=len(np.nonzero(error>tol)[0])
                self.log.write('    {0} of {1} differences exceed tolerance {2}\n'.format(nfail,ntest,tol))

            self.log.write('     Worst difference {0:.4f}\n'.format(np.max(error)))
            print('**** {0} of {1} differences exceed tolerance {2} (max {3:.4f})'.format(nbad,ntest,tolerance,np.max(error)))
        else: 
            self.log.write('Maximum of {0} differences is {1:.4f}\n'.format(ntest,np.max(error)))
            print('** All {1} differences within tolerance {2} (max {3:.4f})'.format(nbad,ntest,tolerance,np.max(error)))

        date=date or ''
        date0=date0 or ''
        comparison=d1[0]+'-'+d0[0]

        csvf=None
        haveheader=False
        if csvsummary:
            csvf=self.qcfilename(csvsummary+'.csv')
            if not append and os.path.exists(csvf):
                os.remove(csvf)
            elif os.path.exists(csvf):
                haveheader=True
            csvw=csv.writer(open(csvf,'ab'))
            header=['date','date0','comparison','ntest']
            header.extend('tol_{0:.3f}'.format(t) for t in tolerances)
            if not haveheader:
                csvw.writerow(header)
            row=[date,date0,comparison,str(ntest)]
            for tol in tolerances[1:]:
                nfail=len(np.nonzero(error>tol)[0])
                row.append(str(nfail))
            csvw.writerow(row)
            csvw=None

        csvf=None
        haveheader=False
        if csvfile:
            csvf=self.qcfilename(csvfile+'.csv')
            if os.path.exists(csvf) and not append:
                os.remove(csvf)
            elif os.path.exists(csvf):
                haveheader=True
        date=date or ''
        date0=date0 or ''
        if (nbad or listall) and csvf:
            if listall:
                badrows=list(range(ntest))
            csvw=csv.writer(open(csvf,'ab'))
            header=['id','date','date0','comparison','lon','lat','diff','dh','dv','de','dn','du','de0','dn0','du0','de1','dn1','du1']
            if not haveheader:
                csvw.writerow(header)
            for ib in badrows:
                row=[str(ib+1),date,date0,comparison,str(points[ib][0]),str(points[ib][1])]
                row.append("{0:.4f}".format(error[ib]))
                row.append("{0:.4f}".format(ds[ib]))
                row.append("{0:.4f}".format(dh[ib]))
                row.extend(("{0:.4f}".format(x) for x in diff[ib]))
                row.extend(("{0:.4f}".format(x) for x in d0[2][ib]))
                row.extend(("{0:.4f}".format(x) for x in d1[2][ib]))
                csvw.writerow(row)
            csvw=None
            self.log.write('     Bad residuals in {0}.csv\n'.format(csvfile))

def calc_total_fault_model( fmodels, points, t=None, t0=None, reverse=False, qclog=None ):
    pfault=np.array([[0.0,0.0,0.0]]*len(points))
    for fm in fmodels:
        pfault += fm.calc_denu(t=t,t0=t0,reverse=reverse,qclog=qclog)
    return pfault

from argparse import ArgumentParser
parser=ArgumentParser('QC checking of deformation models and products - confirm everything matches!')
# parser.add_argument('def_files',nargs='*',help='List of shift .def files from which to construct zip files')
parser.add_argument('--seed','-s',default='Consistent random seed',help='String to seed "random" number generator')
parser.add_argument('--qc-dir','-q',default='qc',help='Directory to write qc results to')
parser.add_argument('--qc-prefix',default='',help='Prefix applied to qc results files')
parser.add_argument('--patch-dir','-p',default='reverse_patch',help='Directory to scan for .def files')
parser.add_argument('--linzdef-file','-d',help='Binary LINZDEF deformation file (SNAP/Landonline)')
parser.add_argument('--binshift-file',help='Binary reverse patch shift model file (Landonline)')
parser.add_argument('--binshift-component',default='HV',help='Components of shift to test (H, V, or HV)')
parser.add_argument('--ntv2-file',help='NTv2 reverse patch file')
parser.add_argument('--model-dir','-m',default='published',help='Directory from which to read the deformation model')
parser.add_argument('--fault-model-dir','-f',default='model',help='Directory from which to read fault models')
parser.add_argument('--land-area','-l',help='WKT file from which to load land areas - used to restrict test points')
parser.add_argument('--fault-zone','-z',help='WKT file from which to load fault zone within which no points are created')
parser.add_argument('--gns-velocities','-v',default='NDM/GNS_2011_V4_velocity_model/solution.gns',help='GNS velocity model file')
parser.add_argument('--ndm-euler-rotation-file','-e',default='ndm_eulerdef',help='File containing NDM euler rotation pole')
parser.add_argument('--n-test-points','-n',type=int,default=ntestpergrid,help='Number of test points per grid area')
parser.add_argument('--test-version',help='Version of model to test (tests components created by the version, default is latest')
parser.add_argument('--test-points',help='Test point file (default is model based semi-random points)')
parser.add_argument('--test-patches',action='store_true',help='Test reverse patch .def files')
parser.add_argument('--test-csvmodel',action='store_true',help='Test published model deformation')
parser.add_argument('--test-binmodel-gns',action='store_true',help='Compare binary file with GNS source as well as with CSV model')
parser.add_argument('--test-velocities',action='store_true',help='Test secular velocity model')
parser.add_argument('--use-ndm-velocities',action='store_true',help='Use velocity model from NDM for patch/dislocation testing')
parser.add_argument('--test-linzdef',action='store_true',help='Test SNAP/Landonline binary model reverse patch')
parser.add_argument('--test-ntv2',action='store_true',help='Test NTv2 format reverse patch')
parser.add_argument('--from-date',help="Start date in format YYYY-MM- D")
parser.add_argument('--to-date',help="End date in format YYYY-MM-DD")
parser.add_argument('--verbose',action='store_true',help="More verbose output")
parser.add_argument('--list-all',action='store_true',help="Include all points in CSV outputs")

args=parser.parse_args()
modeldir=args.model_dir
patchdir=args.patch_dir
lnzdef=args.linzdef_file
binshift=args.binshift_file
verbose=args.verbose
qclog=QcLog(args.qc_dir,verbose=verbose, prefix=args.qc_prefix)

seed(args.seed)

test_patches=args.test_patches
test_csvmodel=args.test_csvmodel
test_velocities=args.test_velocities
test_binmodel=lnzdef is not None
test_binmodel_gns=args.test_binmodel_gns and test_binmodel
test_binshift=binshift is not None
binshift_components=args.binshift_component.upper()
test_ntv2=args.test_ntv2 and args.ntv2_file is not None
ntv2_file=args.ntv2_file

include_areas=[]
if args.land_area is not None:
    include_areas.append((True,load_area( args.land_area, 'land area', qclog)))
fault_zone=None
if args.fault_zone is not None:
    include_areas.append((False,load_area( args.fault_zone,'fault zone',qclog)))

#sys.path.append(os.path.join(modeldir,'tools'))
from LINZ.DeformationModel.Model import Model
ldm=Model(os.path.join(modeldir,'model'))
test_version=ldm.version()
if args.test_version is not None:
    test_version=args.test_version
    if test_version == 'all':
        test_version=None
if test_version is not None:
    qclog.printlog('Testing components added at version {0}'.format(test_version))

# Set up calcdeformation program
cdfm=None
if which('calcdeformation'):
    cdfm=[which('calcdeformation')]
if cdfm is None:
    cdfmf=os.path.join(modeldir,'tools','calcdeformation.py')
    if os.path.exists(cdfmf):
        cdfm=['python',cdfmf]
if cdfm is None:
    raise RuntimeError('Cannot find calcdeformation executable from python-linz-deformationmodel project')
cdfm.extend(('--ndp=6','-m',os.path.join(modeldir,'model')))
cdfm=tuple(cdfm)

if args.test_points:
    if verbose:
        qclog.printlog('Loading test points from',args.test_points)
    points=load_test_points(args.test_points)
else:
    if verbose:
        qclog.printlog('Calculating {0} test points per deformation grid file'.format(args.n_test_points))
    points=test_points(ldm,include_areas=include_areas,points_per_grid=args.n_test_points,version=test_version)
qclog.write_test_points( points )

tolerances=(0.001,0.002,0.005,0.01,0.02,0.05,0.10,0.20,0.50,1.0)

fmodels=None
times=None
time0=None
time1=None
times=[]
calcmodel=test_patches or test_binmodel_gns or test_csvmodel
needtimes=calcmodel or test_binmodel
if needtimes or calcmodel:
    fmodels=fault_models(args.fault_model_dir,points,calcmodel)
    times=test_times(fmodels,version=test_version)
    if args.from_date is not None or args.to_date is not None:
        times=[
            times[0] if args.from_date is None else datetime.strptime(args.from_date,"%Y-%m-%d"),
            times[-1] if args.to_date is None else datetime.strptime(args.to_date,"%Y-%m-%d")
            ]
    time0=times[0]
    time1=times[-1]

# CalcDeformation input file (to be deleted at end of run)

cdin=points_file(points,format="{0},{1},{2},{3}",header="id,lon,lat,hgt")

# Test patch against the published deformation model, and the fault source models.

if test_patches:
    raise RuntimeError('test_patches option not currently functional')
    patches=find_patches(patchdir,points)
    pgroups={}
    for p in patches:
        pdir=os.path.dirname(p.file)
        if pdir not in pgroups:
            pgroups[pdir]=[p]
        else:
            pgroups[pdir].append(p)
    for pdir in pgroups:
        plist=pgroups[pdir]
        pname=os.path.basename(pdir)
        print("Testing reverse patch group",pdir)
        pcsv=os.path.basename(pdir)+'_qc'
        skip=None
        # Don't expect parcel patches to match where less than 0.05.
        if 'parcel' in pdir:
            skip=0.05
        pdenu=plist[0].denu
        for p in plist[1:]:
            pdenu = pdenu+p.denu
        pfile=' '.join((p.file for p in plist))
        qclog.compare("Testing reverse patch file "+pdir+"\n     "+pfile,points,
                      (('flt','Calculated directly from fault models',pfault),
                       ('grd','Calculated from patch file',pdenu)),
                      csvfile=pcsv,
                      tolerance=tolerances,
                      skip=skip
                     )

if test_csvmodel:
    t0=time0
    t0str=t0.strftime('%Y-%m-%d')
    append=False
    for t in times[1:]:
        tstr=t.strftime('%Y-%m-%d')
        print("Testing calcdeformation.py patches from {0} to {1}".format(t0str,tstr))
        pfault=calc_total_fault_model( fmodels, points, t=t, t0=t0, qclog=qclog )
        outfile=tempfilename()
        command=cdfm+('-d',tstr,'-b',t0str,'-o','-ndm',cdin,outfile)
        # print command
        call(command)
        cdenu=read_denu(outfile,cols=(4,7),header=True,expected=len(points),delim=',')
        os.remove(outfile)
        qclog.compare("Testing CalcDeformation.py patches {0} - {1}".format(t0str,tstr),
                          points,
                          (('flt','Calculated directly from fault models',pfault),
                           ('pub','Calculated from published model',cdenu)),
                          csvfile='calcdef_patch',
                          tolerance=tolerances,
                          date=tstr,
                          date0=t0str,
                          listall=args.list_all,
                          append=append
                         )
        append=True

if test_velocities or (test_binmodel and not args.use_ndm_velocities):
    velocities=calc_velocities(args.gns_velocities,args.ndm_euler_rotation_file, points, qclog)

if test_velocities or (test_binmodel and args.use_ndm_velocities):
    outfile=tempfilename()
    command=cdfm+('--date=2100-01-01','--base-date=2000-01-01','--only=ndm',cdin,outfile)
    call(command)
    venu=read_denu(outfile,cols=(4,7),header=True,expected=len(points),delim=',')
    venu /= 100.0
    os.remove(outfile)

if test_velocities:
    print("Testing velocity component in CalcDeformation")
    qclog.compare("Testing CalcDeformation.py velocities",points,
                      (('vel','Calculated directly from gns_velocity+euler',velocities),
                       ('pub','Calculated from published model',venu)),
                      csvfile='ndm_velocities',
                      tolerance=tolerances,
                     )

if test_binmodel:
    if not which('runlnzdef'):
        raise RuntimeError('Require runlnzdef program from the Landonline dbl4u code')
    pfile=points_file(points)
    outfile=tempfilename()
    append=False
    testvel=venu if args.use_ndm_velocities else velocities
    datestr0=None
    mdenu0=None
    bdenu0=None
    cdenu0=None
    call(('runlnzdef',lnzdef,'2000.0',pfile,outfile))
    bdenu2k=read_denu(outfile,cols=(3,6),expected=len(points))
    command=cdfm+('--date=2000-01-01',cdin,outfile)
    call(command)
    cdenu2k=read_denu(outfile,cols=(4,7),header=True,expected=len(points),delim=',')
    os.remove(outfile)

    for dt in times:
        # csvf="dislocations_"+dt.strftime("%Y%m%d")
        csvf="binmodel"
        csvsum="binmodel_summary"
        year=date_as_year(dt)
        datestr="{0:%Y-%m-%d}".format(dt)
        print("Testing dislocations at {0} {1:.2f}".format(datestr,year))
        # Model dislocations
        mdenu=None
        if test_binmodel_gns:
            mdenu=testvel*(year-2000.0)
            for fm in fmodels:
                mdenu += fm.calc_denu(year,qclog=qclog)

        # SNAP/Landonline binary file dislocations
        outfile=tempfilename()
        call(('runlnzdef',lnzdef,str(year),pfile,outfile))
        bdenu=read_denu(outfile,cols=(3,6),expected=len(points))
        os.remove(outfile)

        # CalcDeformation dislocations

        command=cdfm+('--date='+dt.strftime('%Y-%m-%d'),cdin,outfile)
        call(command)
        cdenu=read_denu(outfile,cols=(4,7),header=True,expected=len(points),delim=',')
        os.remove(outfile)
        if test_binmodel_gns:
            comparison=[
                   ('gns','Calculated directly from gns files',mdenu),
                   ('bin','Calculated from SNAP/Landonline binary file',bdenu-bdenu2k),
                   ('pub','Calculated from published model',cdenu-cdenu2k)
                   ]
            qclog.compare("Testing dislocations at {0} relative to 2000-10-01".format(datestr),points,
                          comparison,
                          csvfile=csvf,
                          csvsummary=csvsum,
                          date=datestr,
                          date0='2000-01-01',
                          tolerance=tolerances,
                          append=append
                         )
            append=True
        qclog.compare("Testing dislocations at {0}".format(datestr),points,
                      (('bin','Calculated from SNAP/Landonline binary file',bdenu),
                      ('pub','Calculated from published model',cdenu)),
                      csvfile=csvf,
                      csvsummary=csvsum,
                      date=datestr,
                      tolerance=tolerances,
                      append=append
                     )
        append=True
        if datestr0 is not None:
            comparison=[
                       ('bin','Calculated from SNAP/Landonline binary file',bdenu-bdenu0),
                       ('pub','Calculated from published model',cdenu-cdenu0)
                       ]
            if test_binmodel_gns:
                comparison.insert(0,('gns','Calculated directly from gns files',mdenu-mdenu0))
            qclog.compare("Testing dislocations between  {0} and {1}".format(datestr0,datestr),points,
                          comparison,
                          csvfile=csvf,
                          csvsummary=csvsum,
                          date0=datestr0,
                          date=datestr,
                          tolerance=tolerances,
                          append=append
                         )
            append=True
        datestr0=datestr
        mdenu0=mdenu
        bdenu0=bdenu
        cdenu0=cdenu

# Test reverse patch binary
if test_binshift:
    if not which('runshift'):
        raise RuntimeError('Require runshift program from the Landonline dbl4u code')
    pfile=points_file(points)
    outfile=tempfilename()
    call(('runshift','-p',binshift,pfile,outfile))
    bshift=read_denu(outfile,cols=(3,6),expected=len(points))
    command=cdfm+('-p',cdin,outfile)
    call(command)
    cshift=read_denu(outfile,cols=(4,7),header=True,expected=len(points),delim=',')
    os.remove(outfile)
    if 'H' not in binshift_components:
        bshift[:,:2]=0.0
        cshift[:,:2]=0.0
    if 'V' not in binshift_components:
        bshift[:,2:]=0.0
        cshift[:,2:]=0.0
    csvf='shift'
        
    print("Testing shift model components {0}".format(binshift_components))
    qclog.compare("Testing shift model components {0}".format(binshift_components),points,
                  (('bin','Calculated from Landonline binary shift file',bshift),
                  ('pub','Calculated from published model',cshift)),
                  csvfile=csvf,
                  tolerance=tolerances,
                 )

if test_ntv2:
    ntv2_cvt=os.path.join(tooldir,'ntv2','ntv2_cvt')
    cs2cs=which('cs2cs')
    if not os.path.exists(ntv2_cvt):
        raise RuntimeError('Cannot find ESRI ntv2_cvt at {0}'.format(ntv2_cvt))
    if not cs2cs:
        raise RuntimeError('Cannot find PROJ4 cs2cs')
    # cs2cs likes an absolute path
    gsb=os.path.abspath(ntv2_file)
    if not os.path.exists(gsb):
        raise RuntimeError('NTv2 file {0} does not exists'.format(gsb))
    pfile=points_file(points,format='{1} {2}')
    outfile=tempfilename()

    # Run ESRI ntv2_cfg
    command=[ntv2_cvt,'-r','-p',pfile,gsb]
    ntv2_out=check_output(command)
    eshift=read_denu(source=ntv2_out.split('\n'),
                     cols=(0,2),header=False,
                     relative=points)
    # Run cs2cs
    command=[cs2cs,'-f','%.10f',
            '+proj=longlat', '+ellps=GRS80', '+nadgrids='+gsb,
            '+to',
            '+proj=longlat', '+ellps=GRS80', '+towgs84=0,0,0,0,0,0,0',
            pfile]
    cs2cs_out=check_output(command)
    pshift=read_denu(source=cs2cs_out.split('\n'),
                     cols=(0,2),header=False,
                     relative=points)

    # Run calcdef
    command=cdfm+('-p',cdin,outfile)
    call(command)
    cshift=read_denu(outfile,cols=(4,7),header=True,expected=len(points),delim=',')
    os.remove(outfile)
    os.remove(pfile)

    eshift[:,2:]=0.0
    pshift[:,2:]=0.0
    cshift[:,2:]=0.0
    csvf='ntv2'

    print("Testing NTv2 reverse patch")
    qclog.compare("Testing NTv2 reverse patch",points,
                  (('esr','Calculated from NTv2 with ESRI ntv2_cvt',eshift),
                  ('prj','Calculated from NTv2 with PROJ cs2cs',pshift),
                  ('pub','Calculated from published model',cshift)),
                  csvfile=csvf,
                  tolerance=tolerances,
                 )
