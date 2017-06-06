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
import numpy as np
import ellipsoid
import atexit
from random import uniform, seed
from subprocess import call
from datetime import datetime, timedelta
from collections import namedtuple
from subprocess import call
from ellipsoid import grs80
import euler

calc_okada_image=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'okada','calc_okada')
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

def load_land_area( polygon_file, qclog ):
    from shapely.wkt import loads
    try:
        with open(polygon_file) as laf:
            wkt=laf.read()
            land_area=loads(wkt)
        qclog.printlog( "Using land area definition from "+polygon_file)
        return land_area
    except Exception as ex:
        raise RuntimeError("Cannot load land area definition from "+polygon_file+"\n"+ex.message)

def test_points( ldm, land_area ):
    ranges=[]
    rangestr=set()
    for c in sorted(ldm.components(),key=lambda x: x.name):
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
        for i in range(ntestpergrid):
            tries=20
            while True:
                lon= round(uniform(xmin,xmax),5)
                lat= round(uniform(ymin,ymax),5)
                if land_area is not None:
                    from shapely.geometry import Point
                    if not land_area.contains(Point(lon,lat)):
                        tries -= 1
                        if tries > 0:
                            continue
                        break
                testpts.append(( lon, lat, round(uniform(0.0,1000.0),3)))
                break
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

def read_denu( tmpfile, cols=(2,5),header=False,expected=None,delim=None,fillfrom=None,factor=1.0):
    denu=[]
    ncol=cols[1]
    hline=''
    dline=''
    ndata=cols[1]-cols[0]
    with open(tmpfile) as pf:
        if header:
            hline=pf.readline()
        for l in pf:
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
                denu.append([float(x)*factor for x in parts[cols[0]:cols[1]]])
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

def fault_models( modeldir, points ):

    class FaultModel(namedtuple('FaultModel','file ramp denu')):

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

        def calc_denu(self,t=None,t0=None,reverse=False,verbose=False):
            factor=self.factor(t,reverse)
            if t0 is not None:
                factor -= self.factor(t0,reverse)
            if verbose:
                print "Using {0} times {1:.4f}".format(os.path.basename(self.file),factor)
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
            denu=calc_okada(mfile,points)
            models.append(FaultModel(mfile,timemodel,denu))

    return models


def test_times( models ):
    transitions=set()
    for m in models:
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

    times=[x for x in transitions]
    return sorted(times)

def calc_okada( m, pts ):
    global calc_okada_image
    print "Calculating fault model file",m
    tn1=points_file(pts,format="{1} {2}")
    tn2=tempfilename()
    call((calc_okada_image,m,tn1,tn2))
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

    def __init__( self, qcdir ):
        if not os.path.isdir(qcdir):
            os.makedirs(qcdir)
        qclog=open(os.path.join(qcdir,'qc.log'),'w')
        self.dir=qcdir
        self.log=qclog
        self.log.write("Running deformation QC tests\n")
        self.log.write("Runtime: "+datetime.now().strftime('%Y-%m-%d %H:%M:%S')+'\n')

    def printlog( self, *messages ):
        self.log.write(' '.join(messages))
        self.log.write('\n')

    def start_test(self,message):
        self.log.write("\n")
        self.log.write(message)
        self.log.write("\n")

    def write_test_points(self,points):
        csvw=csv.writer(open(os.path.join(self.dir,'test_points.csv'),'w'))
        csvw.writerow(('id','lon','lat'))
        for i,p in enumerate(points):
            csvw.writerow((str(i),str(p[0]),str(p[1])))

    def compare(self,
                description,
                points,
                datasets,
                csvfile=None,
                tolerance=0.001,
                skip=None,
                listall=False,
                date0=None,
                date=None,
                append=False
               ):
        try:
            tolerances=list(tolerance)
        except:
            tolerances=(tolerance,)
        tolerance=tolerances[0]

        self.start_test(description)
        self.log.write('Comparing:\n')
        for i,d in enumerate(datasets):
            self.log.write("    {0}: {1}\n".format(d[0],d[1]))

        # Find rows exceeding tolerance

        error=np.array([0.0]*len(points))
        for id0,d0 in enumerate(datasets[:-1]):
            for d1 in datasets[id0+1:]:
                error=np.maximum(error,np.sqrt(np.sum(np.square(d0[2]-d1[2]),axis=1)))

        skiprows=[]
        if skip is not None:
            skiprows=np.nonzero(np.sum(np.square(datasets[0][2]),axis=1) < skip*skip)[0]
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
            print '**** {0} of {1} differences exceed tolerance {2} (max {3:.4f})'.format(nbad,ntest,tolerance,np.max(error))
        else: 
            self.log.write('Maximum of {0} differences is {1:.4f}\n'.format(ntest,np.max(error)))

        csvf=None
        haveheader=False
        if csvfile:
            csvf=os.path.join(self.dir,csvfile+'.csv')
            if not append:
                try:
                    os.path.remove(csvf)
                except:
                    pass
            elif os.path.exists(csvf):
                haveheader=True
        date=date or ''
        date0=date0 or ''
        if (nbad or listall) and csvf:
            if listall:
                badrows=range(ntest)
            csvw=csv.writer(open(csvf,'ab' if append else 'wb'))
            header=['id','date','date0','lon','lat','diff']
            for i in range(len(datasets)):
                stri=datasets[i][0]
                header += ['de_'+stri,'dn_'+stri,'du_'+stri]
            if not haveheader:
                csvw.writerow(header)
            for ib in badrows:
                row=[str(ib+1),date,date0,str(points[ib][0]),str(points[ib][1])]
                data=[error[ib]]
                for ds in datasets:
                    data.extend(ds[2][ib])
                row.extend(("{0:.4f}".format(x) for x in data))
                csvw.writerow(row)
            csvw=None
            self.log.write('     Bad residuals in {0}.csv\n'.format(csvfile))

def calc_total_fault_model( fmodels, points, t=None, t0=None, reverse=False, verbose=False ):
    pfault=np.array([[0.0,0.0,0.0]]*len(points))
    for fm in fmodels:
        pfault += fm.calc_denu(t=t,t0=t0,reverse=reverse,verbose=verbose)
    return pfault

from argparse import ArgumentParser
parser=ArgumentParser('QC checking of deformation models and products - confirm everything matches!')
# parser.add_argument('def_files',nargs='*',help='List of shift .def files from which to construct zip files')
parser.add_argument('--seed','-s',default='Consistent random seed',help='String to seed "random" number generator')
parser.add_argument('--qc-dir','-q',default='qc',help='Directory to write qc results to')
parser.add_argument('--patch-dir','-p',default='reverse_patch',help='Directory to scan for .def files')
parser.add_argument('--linzdef-file','-d',help='Binary LINZDEF deformation file (SNAP/Landonline)')
parser.add_argument('--model-dir','-m',default='published',help='Directory from which to read the deformation model')
parser.add_argument('--fault-model-dir','-f',default='model',help='Directory from which to read fault models')
parser.add_argument('--land-area','-l',help='WKT file from which to load land areas - used to restrict test points')
parser.add_argument('--gns-velocities','-v',default='NDM/GNS_2011_V4_velocity_model/solution.gns',help='GNS velocity model file')
parser.add_argument('--ndm-euler-rotation-file','-e',default='ndm_eulerdef',help='File containing NDM euler rotation pole')
parser.add_argument('--n-test-points','-n',type=int,default=ntestpergrid,help='Number of test points per grid area')
parser.add_argument('--test-patches',action='store_true',help='Test reverse patch .def files')
parser.add_argument('--test-calcdef',action='store_true',help='Test published model deformation')
parser.add_argument('--test-velocities',action='store_true',help='Test secular velocity model')
parser.add_argument('--use-ndm-velocities',action='store_true',help='Use velocity model from NDM for patch/dislocation testing')
parser.add_argument('--test-linzdef',action='store_true',help='Test SNAP/Landonline binary model dislocations')
parser.add_argument('--from-date',help="Start date in format YYYY-MM-DD")
parser.add_argument('--to-date',help="End date in format YYYY-MM-DD")
parser.add_argument('--verbose',action='store_true',help="More verbose output")
parser.add_argument('--list-all',action='store_true',help="Include all points in CSV outputs")

args=parser.parse_args()
modeldir=args.model_dir
patchdir=args.patch_dir
lnzdef=args.linzdef_file
qclog=QcLog(args.qc_dir)

seed(args.seed)

test_patches=args.test_patches
test_calcdef_patch=args.test_calcdef
test_velocities=args.test_velocities
test_dislocations=lnzdef is not None

land_area=None
if args.land_area is not None:
    land_area=load_land_area( args.land_area, qclog)

#sys.path.append(os.path.join(modeldir,'tools'))
from LINZ.DeformationModel.Model import Model
ldm=Model(os.path.join(modeldir,'model'))

cdfm=None
if os.path.exists('/usr/bin/calcdeformation'):
    cdfm=['/usr/bin/calcdeformation']
else:
    cdfmf=os.path.join(modeldir,'tools','calcdeformation.py')
    if os.path.exists(cdfmf):
        cdfm=['python',cdfmf]
if cdfm is not None:
    cdfm.extend(('--ndp=6','-m',os.path.join(modeldir,'model')))
    cdfm=tuple(cdfm)

points=test_points(ldm,land_area)
qclog.write_test_points( points )

fmodels=None
times=None
time0=None
time1=None
if test_patches or test_dislocations:
    fmodels=fault_models(args.fault_model_dir,points)
    times=test_times(fmodels)
    time0=datetime(2000,1,1) if args.from_date is None else datetime.strptime(args.from_date,"%Y-%m-%d")
    time1=times[-1] if args.to_date is None else datetime.strptime(args.to_date,"%Y-%m-%d")

# CalcDeformation input file (to be deleted at end of run)

cdin=points_file(points,format="{0},{1},{2},{3}",header="id,lon,lat,hgt")

# Test patch against the published deformation model, and the fault source models.

# Compile 

if test_patches:
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
        print "Testing reverse patch group",pdir
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
                      tolerance=(0.001,0.002,0.005,0.01),
                      skip=skip
                     )

if test_calcdef_patch:
    t=time1
    t0=time0
    tstr=t.strftime('%Y-%m-%d')
    t0str=t0.strftime('%Y-%m-%d')
    print "Testing calcdeformation.py patches from {0} to {1}".format(t0str,tstr)
    pfault=calc_total_fault_model( fmodels, points, t=t, t0=t0, verbose=args.verbose )
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
                      tolerance=(0.001,0.002,0.005,0.01),
                      listall=args.list_all
                     )


if test_velocities or (test_dislocations and not args.use_ndm_velocities):
    velocities=calc_velocities(args.gns_velocities,args.ndm_euler_rotation_file, points, qclog)

if test_velocities or (test_dislocations and args.use_ndm_velocities):
    outfile=tempfilename()
    command=cdfm+('--date=2100-01-01','--base-date=2000-01-01','--only=ndm',cdin,outfile)
    call(command)
    venu=read_denu(outfile,cols=(4,7),header=True,expected=len(points),delim=',')
    venu /= 100.0
    os.remove(outfile)

if test_velocities:
    print "Testing velocity component in CalcDeformation"
    qclog.compare("Testing CalcDeformation.py velocities",points,
                      (('vel','Calculated directly from gns_velocity+euler',velocities),
                       ('pub','Calculated from published model',venu)),
                      csvfile='ndm_velocities',
                      tolerance=(0.0001,0.0002,0.0005,0.001)
                     )

if test_dislocations:
    pfile=points_file(points)
    append=False
    testvel=venu if args.use_ndm_velocities else velocities
    datestr0=None
    mdenu0=None
    bdenu0=None
    cdenu0=None
    for dt in times:
        # csvf="dislocations_"+dt.strftime("%Y%m%d")
        csvf="dislocations"
        year=date_as_year(dt)
        datestr="{0:%Y-%m-%d}".format(dt)
        print "Testing dislocations at {0} {1:.2f}".format(datestr,year)
        # Model dislocations
        mdenu=testvel*(year-2000.0)
        for fm in fmodels:
            mdenu += fm.calc_denu(year,reverse=True)

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
        qclog.compare("Testing dislocations at {0}".format(datestr),points,
                      (('gns','Calculated directly from gns files',mdenu),
                      ('bin','Calculated from SNAP/Landoline binary file',bdenu),
                      ('pub','Calculated from published model',cdenu)),
                      csvfile=csvf,
                      date=datestr,
                      tolerance=(0.001,0.002,0.005,0.01),
                      append=append
                     )
        append=True
        if datestr0 is not None:
            qclog.compare("Testing dislocations between  {0} and {1}".format(datestr0,datestr),points,
                          (('gns','Calculated directly from gns files',mdenu-mdenu0),
                          ('bin','Calculated from SNAP/Landoline binary file',bdenu-bdenu0),
                          ('pub','Calculated from published model',cdenu-cdenu0)),
                          csvfile=csvf,
                          date0=datestr0,
                          date=datestr,
                          tolerance=(0.001,0.002,0.005,0.01),
                          append=append
                         )
        datestr0=datestr
        mdenu0=mdenu
        bdenu0=bdenu
        cdenu0=cdenu
