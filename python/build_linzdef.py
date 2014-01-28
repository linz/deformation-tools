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
if not isdir(md) or not isdir(joinpath(md,'model')) or not isdir(joinpath(md,'tools')):
    raise RuntimeError("Invalid model directory: "+md)

bd=args.target_dir
if not isdir(bd):
    if os.path.exists(bd):
        raise RuntimeError("Invalid target build directory: "+bd)
    else:
        os.makedirs(bd)

defname=args.linzdef_file
deffile=open(joinpath(bd,defname+'.def'),"w")

sys.path.append(joinpath(md,'tools'))
from LINZ.DeformationModel import Model, Time

m=Model.Model(joinpath(md,'model'))

header= '''
DEFORMATION_MODEL NZGD2000 deformation model
FORMAT LINZDEF1B
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
    def __init__(self,time,f0,f1):
        self.time=time
        self.days=time.daysAfter(refdate)
        self.f0=f0
        self.f1=f1

    def __str__(self):
        return 'Time:{0} f0:{1} f1:{2}'.format(self.time,self.f0,self.f1)

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
    zerobeyond='no'
    for m in c.spatialModel.models():
        if m.spatial_model != 'llgrid':
            raise RuntimeError('Cannot handle spatial model type '+spatial_model)
        # print m.model().gridSpec()
        #print '    ',m.model().gridFile()
        #print '    ',m.columns
        grids.append(m.model().gridFile())
        dimension=len(m.columns)
        zerobeyond='yes' if m.spatial_complete else 'no'

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
DATA_TYPE {dtype}
DIMENSION {dimension}
START_DATE 1-Jan-1850
END_DATE 1-Jan-2100
ZERO_BEYOND_RANGE {zerobeyondrange}
DESCRIPTION
{description}
END_DESCRIPTION
'''

griddef= '''

DEFORMATION_COMPONENT {gridfile}
MODEL_TYPE grid
REF_DATE {refdate}
BEFORE_REF_DATE {before}
AFTER_REF_DATE {after}
DESCRIPTION
{description}
END_DESCRIPTION
'''

SeqComp=namedtuple('SeqComp','mtype time factor before after')
small=0.00001

for sequence in sequences:
    compname=sequence.component
    print "Sequence:",compname
    basegrids=[]
    for i,g in enumerate(sequence.grids):
        gfile=joinpath(md,'model',g)
        bfile=joinpath(bd,'basegrid'+str(i)+'.csv')
        cmd=['gridtool','read','csv',gfile]
        for b in basegrids:
            cmd.extend(['subtract','csv',b])
        cmd.extend(['write','csv',bfile])
        call(cmd)
        basegrids.append(bfile)

    # print "Grids\n  "+"\n  ".join(basegrids)

    # Construct the time sequence... Really messy stuff as new deformation
    # model doesn't really fit format used in Landonline and SNAP...
    # ... for the moment!

    # Build a series of time events with before and after values...
    # Velocities are handled differently, so put in their own list

    subsequences = []

    velocities=[]
    events=[]
    for s in sequence.steps:
        if s.mtype == 'velocity':
            subsequences.append([SeqComp('velocity',s.t0,1.0,'interpolate','interpolate')])
            continue
        if s.mtype == 'step':
            if len(events) == 0:
                events.append(TimeEvent(s.t0,s.f0,s.f1))
                continue
            da=s.t0.daysAfter(refdate)
            offset=s.f0
            found = False
            for i,e in enumerate(list(events)):
                # Event before step
                if found or e.days < da-0.5:
                    e.f0 += offset
                    e.f1 += offset
                # Event coincident with step
                elif e.days < da+0.5:
                    e.f0 += offset
                    offset=s.f1
                    e.f1 += offset
                    found = True
                # Step before next event
                else:
                    if i == 0:
                        offset=e.f0
                    else:
                        e0=events[i-1]
                        offset=(e0.f1*(e.days-da)+d1.f0*(da-e0.days))/(e.days-e0.days)
                    events.insert(i,TimeEvent(s.t0,s.f0+offset,s.f1+offset))
                    offset=s.f1
                    e.f0 += offset
                    e.f1 += offset
                    found = True
            continue
        if s.mtype == 'ramp':
            if len(events) == 0:
                events.append(TimeEvent(s.t0,s.f0,s.f0))
                events.append(TimeEvent(s.t1,s.f1,s.f1))
                continue
            d0=s.t0.daysAfter(refdate)
            d1=s.t1.daysAfter(refdate)
            found0=False
            found1=False
            insert0=-1
            insert1=-1
            offset0=None
            offset1=None
            de0=events[0].days-3000.0
            fe0=events[0].f0

            for i, e in enumerate(events):
                de=e.days
                if de > d0-0.5 and not found0:
                    if de > d0+0.5:
                        insert0=i
                        offset0=(fe0*(de-d0)+e.f0*(d0-de0))/(de-de0)
                    found0=True
                if de > d1-0.5 and not found1:
                    if de > d1+0.5:
                        insert1=i
                        offset1=(fe0*(de-d1)+e.f0*(d1-de0))/(de-de0)
                    found1=True
                fe0=e.f1
                de0=de
                
                if de < d0:
                    offset=s.f0
                elif de > d1:
                    offset=s.f1
                else:
                    offset=(s.f0*(d1-de)+s.f1*(de-d0))/(d1-d0)
                e.f0 += offset
                e.f1 += offset
            if insert1 >= 0:
                events.insert(insert1,TimeEvent(s.t1,s.f1+offset1,s.f1+offset1))
            if insert0 >= 0:
                events.insert(insert0,TimeEvent(s.t0,s.f0+offset0,s.f0+offset0))
            if not found0:
                events.append(TimeEvent(s.t0,s.f0+fe0,s.f0+fe0))
            if not found1:
                events.append(TimeEvent(s.t1,s.f1+fe0,s.f1+fe0))

    # Define the deformation sequences for each event
    if events:
        lastevent=''
        subseq=[]
        for i,e in enumerate(events):
            # print "Event ",e
            split=abs(e.f0 - e.f1) > small
            before='interpolate' if i>0 else 'zero' if abs(e.f0) < small else 'fixed'
            after='interpolate' if i<len(events)-1 else 'zero' if abs(e.f1) < small else 'fixed'
            if abs(e.f0 - e.f1) > small and before != 'zero' and after != 'zero':
                subseq.append(SeqComp('deformation',e.time,e.f0,before,'zero'))
                subsequences.append(subseq)
                subseq=[]
                before='zero'
            subseq.append(SeqComp('deformation',
                                  e.time,
                                  e.f0 if after == 'zero' else e.f1,
                                  before,after))
        subsequences.append(subseq)

    # Now construct the sequences in the definition file...

    factors=[]
    zgrids=[]

    for nsseq,sseq in enumerate(subsequences):
        # print "Subsequence {0}, {1}".format(nsseq,sseq)
        sseqname='_s{0}'.format(nsseq+1) if len(subsequences) > 1 else ''
        for ngf,gf in enumerate(basegrids):
            gfname='_g{0}'.format(ngf+1) if len(basegrids) > 1 else ''
            deffile.write(seqdef.format(
                name=sequence.component+sseqname+gfname,
                dtype=sseq[0].mtype,
                dimension=sequence.dimension,
                zerobeyondrange=sequence.zerobeyond,
                description=sequence.component
            ))
            for iss,seq in enumerate(sseq):
                seqname='_ss{0}'.format(iss) if len(sseq) > 1 else ''
                sgfname=defname+'_'+compname+sseqname+gfname+seqname+'.gdf'
                bgrid=gf
                if abs(seq.factor) < small:
                    zgrid=joinpath(bd,'zerogrid'+str(ngf)+'.csv')
                    if zgrid not in zgrids:
                        with open(bgrid,"r") as f:
                            hl=f.readline()
                            l0=f.readline()
                            l1=l0
                            for l in f:
                                if ',' in l:
                                    l1 = l
                        with open(zgrid,"w") as f:
                            f.write(hl)
                            parts=l0.strip().split(',')
                            lon0,lat0=parts[:2]
                            parts=l1.strip().split(',')
                            lon1,lat1=parts[:2]
                            suffix=','.join(['0.0']*(len(parts)-2))
                            f.write('{0},{1},{2}\n'.format(lon0,lat0,suffix))
                            f.write('{0},{1},{2}\n'.format(lon1,lat0,suffix))
                            f.write('{0},{1},{2}\n'.format(lon0,lat1,suffix))
                            f.write('{0},{1},{2}\n'.format(lon1,lat1,suffix))
                        zgrids.append(zgrid)
                    bgrid=zgrid
                
                cmd=['gridtool','read','csv',bgrid]
                cmd.extend(['multiply',str(seq.factor)])
                cmd.extend(['write_linzgrid','NZGD2000','Deformation grid',sgfname,''])
                cmd.append(joinpath(bd,sgfname))
                call(cmd)
                deffile.write(griddef.format(
                    gridfile=sgfname,
                    refdate=seq.time.strftime('%d-%b-%Y'),
                    before=seq.before,
                    after=seq.after,
                    description='Deformation grid',
                ))

    for bg in basegrids:
        os.remove(bg)
    for zg in zgrids:
        os.remove(zg)

deffile.close()

print "Building binary deformation files"
os.chdir(bd)
call(['makelinzdefmodel.pl','-f','LINZDEF1B',defname+'.def',defname+'b.bin'])
call(['makelinzdefmodel.pl','-f','LINZDEF1L',defname+'.def',defname+'l.bin'])

