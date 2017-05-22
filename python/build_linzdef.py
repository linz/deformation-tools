#!/usr/bin/python

# Script to build a LINZDEF deformaton model file from the published format
#
# Very messy script as new format is not a perfect fit with the format used to build the files
# currently used by SNAP and Landonline.  This script does not handle the full scope of the
# new model, but is sufficient for current requirements.  Long term propose to migrate SNAP and 
# Landonline to use the new format.

import os
import os.path
from os.path import join as joinpath
from os.path import isdir
from collections import namedtuple
import sys
import argparse
from datetime import datetime
from subprocess import call

argparser=argparse.ArgumentParser('Build a LINZDEF deformation file from a published deformation model')
argparser.add_argument('model_dir',help='Base directory of model, containing model and tools directories')
argparser.add_argument('target_dir',help='Target directory in which to build model')
argparser.add_argument('linzdef_file',help='Name of final model')

args = argparser.parse_args()

md=args.model_dir
if not isdir(md) or not isdir(joinpath(md,'model')):
    raise RuntimeError("Invalid model directory: "+md)

bd=args.target_dir
if not isdir(bd):
    if os.path.exists(bd):
        raise RuntimeError("Invalid target build directory: "+bd)
    else:
        os.makedirs(bd)

defname=args.linzdef_file
deffile=open(joinpath(bd,defname+'.def'),"w")

toolsdir=joinpath(md,'tools')
if os.path.exists(toolsdir):
    sys.path.append(joinpath(md,toolsdir))

try:
    from LINZ.DeformationModel import Model, Time
except ImportError:
    raise RuntimeError('Model tools directory not found and LINZ.DeformationModel not intalled')

m=Model.Model(joinpath(md,'model'))

header= '''
DEFORMATION_MODEL NZGD2000 deformation model
FORMAT LINZDEF2B
VERSION_NUMBER {version}
VERSION_DATE  {versiondate}
START_DATE 1-Jan-1850
END_DATE 1-Jan-2200
COORDSYS NZGD2000
DESCRIPTION
{description}
END_DESCRIPTION
'''.format(
    version=m.version(),
    versiondate=m.versionInfo(m.version()).release_date.strftime('%d-%b-%Y'),
    description=m.description(submodels=False)
    )

# print header
deffile.write(header)
deffile.write("\n")

TimeStep=namedtuple('TimeStep','mtype t0 f0 t1 f1')
DefSeq=namedtuple('DefSeq','component dimension zerobeyond steps grids')

DefComp=namedtuple('DefComp','date factor before after')

refdate=Time.Time(datetime(2000,1,1))
class TimeEvent:
    def __init__(self,time,f0,time1,f1):
        self.time0=time
        self.time1=time1
        self.days=time.daysAfter(refdate)
        self.f0=f0
        self.f1=f1

    def __str__(self):
        return 'Time:{0} f0:{1} {2} f1:{3}'.format(self.time0,self.f0,self.time1,self.f1)

    def __cmp__( self, other ):
        return cmp( (self.time0,self.time1), (other.time0,other.time1) )

    def __repr__(self):
        return str(self)

sequences=[]

for c in m.components():
    print "Adding component:",c.name
    component=c.submodel
    # print component
    tm=c.timeFunction
    mtype=tm.time_function
    if mtype not in ('velocity','step','ramp'):
        raise RuntimeError('Cannot handle temporal model type '+mtype)
    # print type(c.timeComponent.model())
    # print mtype,tm.factor0,tm.factor1,tm.time0,tm.time1

    grids=[]
    dimension=0
    zerobeyond='yes'
    for m in c.spatialModel.models():
        if m.spatial_model != 'llgrid':
            raise RuntimeError('Cannot handle spatial model type '+spatial_model)
        # print m.model().gridSpec()
        #print '    ',m.model().gridFile()
        #print '    ',m.columns
        grids.append(m.model().gridFile())
        dimension=len(m.columns)
        if not m.spatial_complete:
            zerobeyond='no'

    # Reverse grids so that contained grids occur before containing grids..
    grids.reverse()
    step=TimeStep(mtype,tm.time0,tm.factor0,tm.time1,tm.factor1)

    found = False
    for s in sequences:
        if (s.component == component and
            s.zerobeyond == zerobeyond and
            s.grids == grids ):
            s.steps.append(step)
            found=True
            break

    if not found:
        sequences.append(DefSeq(component,dimension,zerobeyond,[step],grids))

# print sequences

seqdef='''

DEFORMATION_SEQUENCE {name}
DIMENSION {dimension}
START_DATE 1-Jan-1850
END_DATE 1-Jan-2100
ZERO_BEYOND_RANGE {zerobeyondrange}
NESTED_SEQUENCE yes
DESCRIPTION
{description}
END_DESCRIPTION
'''

griddef= '''

DEFORMATION_COMPONENT {gridfile}
MODEL_TYPE grid
REF_DATE {refdate}
TIME_MODEL {tmodel}
DESCRIPTION
{description}
END_DESCRIPTION
'''

SeqComp=namedtuple('SeqComp','time factor before after nested')
small=0.00001
gdf_files={}
seqcomps={}

for sequence in sequences:
    if sequence.component not in seqcomps:
        seqcomps[sequence.component] = 0
    seqcomps[sequence.component] += 1

seqcomps={c:'' for c in seqcomps if seqcomps[c] > 1}

for sequence in sequences:
    compname=sequence.component
    ncomp=0
    while compname in seqcomps:
        ncomp += 1
        compname=sequence.component + '_c{0}'.format(ncomp)
    seqcomps[compname]=sequence.component

    print "Sequence:",compname

    subsequences = []
    events=[]
    for s in sequence.steps:
        if s.mtype == 'velocity':
            subsequences.append([s.t0,'VELOCITY'])
        elif s.mtype == 'step':
            events.append(TimeEvent(s.t0,s.f0,s.t0,s.f1))
        elif s.mtype == 'ramp':
            events.append(TimeEvent(s.t0,s.f0,s.t1,s.f1))
        else:
            raise RuntimeError('Cannot handle time model type '+s.mtype)

    pwm=''
    if len(events) > 0:
        time_model=[]
        events.sort()
        e0=None
        for e in events:
            if e0 and e0.time1.daysAfter(e.time0) > 0.001:
                raise RuntimeError('Cannot handle overlapping time events in series')
            v0 = 0.0 if len(time_model) == 0 else time_model[-1][1]
            for t in time_model:
                t[1] += e.f0
            time_model.append([e.time0,e.f0+v0])
            time_model.append([e.time1,e.f1+v0])
        for i in reversed(range(len(time_model)-1)):
             t0=time_model[i]
             t1=time_model[i+1]
             if abs(t0[0].daysAfter(t1[0])) < 0.001 and abs(t0[1]-t1[1]) < 0.00001:
                 time_model[i:i+1]=[]
        pwm='PIECEWISE_LINEAR {0}'.format(time_model[0][1])
        i0=1 if time_model[1][0].daysAfter(time_model[0][0]) < 0.001 else 0
        for t in time_model[i0:]:
            pwm = pwm+" {0} {1}".format( t[0].strftime('%d-%b-%Y'),t[1] )

        subsequences.append([refdate,pwm])

    # Now construct the sequences in the definition file...

    factors=[]
    zgrids=[]

    for nsseq,sseq in enumerate(subsequences):
        # print "Subsequence {0}, {1}".format(nsseq,sseq)
        sseqname='_s{0}'.format(nsseq+1) if len(subsequences) > 1 else ''
        deffile.write(seqdef.format(
            name=sequence.component+sseqname,
            dimension=sequence.dimension,
            zerobeyondrange=sequence.zerobeyond,
            description=sequence.component
        ))

        for ngf,g in enumerate(sequence.grids):
            gf=joinpath(md,'model',g)
            if gf not in gdf_files:
                gfname='_g{0}'.format(ngf+1) if len(sequence.grids) > 1 else ''
                sgfbasename=defname+'_'+compname+sseqname+gfname
                sgfname=sgfbasename+'.gdf'
                nsgf=0
                while sgfname in gdf_files.values():
                    nsgf += 1
                    sgfname=sgfbasename+'_f{0}'.format(nsgf)+'.gdf'
                cmd=['gridtool','read','csv',gf]
                cmd.extend(['write_linzgrid','NZGD2000','Deformation grid',sgfname,''])
                cmd.append(joinpath(bd,sgfname))
                call(cmd)
                gdf_files[gf]=sgfname
            sgfname=gdf_files[gf]
            deffile.write(griddef.format(
                gridfile=sgfname,
                refdate=sseq[0].strftime('%d-%b-%Y'),
                description='Deformation grid',
                tmodel=sseq[1],
            ))

deffile.close()

print "Building binary deformation files"
os.chdir(bd)
call(['makelinzdefmodel.pl','-f','LINZDEF2B',defname+'.def',defname+'b.bin'])
call(['makelinzdefmodel.pl','-f','LINZDEF2L',defname+'.def',defname+'l.bin'])

