#!/usr/bin/python

# Script to create the master file GeoTIFF file creation scripts to encode 
# the LINZ Deformation model

from __future__ import print_function

from collections import namedtuple, OrderedDict
from datetime import datetime
from os.path import isdir
from os.path import join as joinpath
from subprocess import call
import argparse
import hashlib
import os
import os.path
import re
import json
import sys
from LINZ.DeformationModel import Model, Time

argparser=argparse.ArgumentParser('Build a LINZDEF deformation file from a published deformation model')
argparser.add_argument('model_dir',help='Base directory of model, containing model and tools directories')
argparser.add_argument('target_dir',help='Target directory in which to build model')
argparser.add_argument('def_model_base',help='Base name of final model')
argparser.add_argument('-a','--all-versions',action='store_true',help='Create master files for all versions of model')

args = argparser.parse_args()
allversions=args.all_versions

gridfields=[
    {'none':[],'horizontal':['horizontal'],'vertical':['vertical'],'3d':['horizontal','vertical']},
    {'none':[],'horizontal':['horizontal_uncertainty'],'vertical_uncertainty':['vertical'],'3d':['horizontal_uncertainty','vertical_uncertainty']},
]

component_uncertainty=[0.01,0.01]

md=args.model_dir
if not isdir(md):
    raise RuntimeError("Invalid model directory: "+md)

bd=args.target_dir
if not isdir(bd):
    if os.path.exists(bd):
        raise RuntimeError("Invalid target build directory: "+bd)
    else:
        os.makedirs(bd)
defname=args.def_model_base

m=Model.Model(joinpath(md,'model'))

clean_description=lambda d: re.sub(r'[\r\n]+','\n',d.strip())

TimeStep=namedtuple('TimeStep','mtype t0 f0 t1 f1')
DefSeq=namedtuple('DefSeq','component description gridtype versionstart versionend zerobeyond steps grids extent')
DefComp=namedtuple('DefComp','date factor before after')

timeformat='%Y-%m-%dT00:00:00Z'
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

class Extent:
    def __init__(self,minlon,maxlon,minlat,maxlat):
        self._minlon=minlon
        self._maxlon=maxlon
        self._minlat=minlat
        self._maxlat=maxlat

    def copy( self ):
        return Extent(self._minlon,self._maxlon,self._minlat,self._maxlat)
    
    def unionWith( self, other ):
        if other._minlon < self._minlon:
            self._minlon=other._minlon
        if other._maxlon > self._maxlon:
            self._maxlon=other._maxlon
        if other._minlat < self._minlat:
            self._minlat=other._minlat      
        if other._maxlat > self._maxlat:
            self._maxlat=other._maxlat  

    def spec( self ):
        return [self._minlon,self._minlat,self._maxlon,self._maxlat]

sequences=[]

for c in m.components(allversions=allversions):
    print("Adding source component:",c.name)
    component=c.submodel
    # print component
    tm=c.timeFunction
    mtype=tm.time_function
    if mtype not in ('velocity','step','ramp'):
        raise RuntimeError('Cannot handle temporal model type '+mtype)
    # print type(c.timeComponent.model())
    # print mtype,tm.factor0,tm.factor1,tm.time0,tm.time1

    grids=[]
    gridtype='none'
    zerobeyond='yes'
    if c.versionRevoked != '0' and not allversions:
        continue
    versionstart=c.versionAdded
    versionend=c.versionRevoked
    extent=None    
    descriptions=OrderedDict()
    for sm in c.spatialModel.models():
        if sm.spatial_model != 'llgrid':
            raise RuntimeError('Cannot handle spatial model type '+sm.spatial_model)
        descriptions[sm.description]=1
        # print sm.model().gridSpec()
        #print '    ',sm.model().gridFile()
        #print '    ',sm.columns
        gridtype=sm.displacement_type+':'+sm.error_type
        grids.append(sm.model().gridFile())
        if not sm.spatial_complete:
            zerobeyond='no'
        mextent=Extent(sm.min_lon,sm.max_lon,sm.min_lat,sm.max_lat)
        if extent is None:
            extent=mextent
        else:
            extent.unionWith(mextent)

    description="\n".join([d for d in descriptions])
    description=description.replace("\r","")

    step=TimeStep(mtype,tm.time0,tm.factor0,tm.time1,tm.factor1)

    found = False
    for s in sequences:
        if (s.component == component and
            s.zerobeyond == zerobeyond and
            s.versionstart == versionstart and
            s.versionend == versionend and
            s.grids == grids and 
            s.gridtype == gridtype  ):
            s.steps.append(step)
            s.extent.unionWith(extent)
            found=True
            break

    if not found:
        sequences.append(DefSeq(component,description,gridtype,versionstart,versionend,zerobeyond,[step],grids,extent))

# Compile components and versions

SeqComp=namedtuple('SeqComp','time factor before after nested')
small=0.00001
gdf_files={}
seqcomps={}

for sequence in sequences:
    if sequence.component not in seqcomps:
        seqcomps[sequence.component] = 0
    seqcomps[sequence.component] += 1

seqcomps={c:'' for c in seqcomps if seqcomps[c] > 1}

components=[]
grids={}

for sequence in sequences:
    compname=sequence.component
    ncomp=0
    while compname in seqcomps:
        ncomp += 1
        compname=sequence.component + '_c{0}'.format(ncomp)
    seqcomps[compname]=sequence.component

    print("Building proj component:",compname)

    subsequences = []
    events = []
    timefuncs=[]

    
    for s in sequence.steps:
        if s.mtype == 'velocity':
            timefuncs.append((OrderedDict([
                ('type', 'velocity'),
                ('parameters',{'reference_epoch': s.t0.strftime(timeformat)})
            ]),s.t0))
        elif s.mtype == 'step':
            events.append(TimeEvent(s.t0,s.f0,s.t0,s.f1))
        elif s.mtype == 'ramp':
            events.append(TimeEvent(s.t0,s.f0,s.t1,s.f1))
        else:
            raise RuntimeError('Cannot handle time model type '+s.mtype)

    time_model=[]
    eventtime=None
    if len(events) > 0:
        time_model=[]
        events.sort()
        eventtime=events[0].time0
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

    if len(time_model) == 2 and time_model[0][0] == time_model[1][0] and time_model[0][1] == 0.0 and time_model[1][1] == 1.0 and time_model[0][0]:
        timefuncs.append((OrderedDict([
            ('type','step'),
            ('parameters',{'step_epoch':time_model[0][0].strftime(timeformat)})
        ]),eventtime))
    elif len(time_model) == 2 and time_model[0][0] == time_model[1][0] and time_model[0][1] == -1.0 and time_model[1][1] == 0.0 and time_model[0][0]:
        timefuncs.append((OrderedDict([
            ('type','reverse_step'),
            ('parameters',{'step_epoch':time_model[0][0].strftime(timeformat)})
        ]),eventtime))    
    elif len(time_model) > 1:
        timefuncs.append((OrderedDict([
            ('type','piecewise'),
            ('parameters', OrderedDict([
            ('before_first', 'zero' if time_model[0][1] == 0.0 else 'constant'),
            ('after_last', 'zero' if time_model[-1][1] == 0.0 else 'constant'),
            ('model', [OrderedDict([('epoch',t[0].strftime(timeformat)),('scale_factor',round(t[1],4))]) for t in time_model])
            ]))
        ]),eventtime))

    # Now construct the sequences in the definition file...

    fields=[]
    for f,t in zip(gridfields,sequence.gridtype.split(':')):
        fields.extend(f[t])

    gridspec=None
    gkey=':'.join(sorted(sequence.grids))
    if gkey in grids:
        gridspec=grids[gkey]
    else:
        gname=sequence.grids[0]
        gname=os.path.dirname(gname)
        gname=os.path.basename(gname)
        gname=gname.replace('patch_','')
        gname=gname.replace('_','')
        gname=defname+'-'+gname
        ngname=0
        while True:
            ngname += 1
            gridname=gname+'-grid{0:02d}'.format(ngname)+'.tif'
            if gridname not in grids:
                break
        gridfile=os.path.join(bd,gridname)
        with open(gridfile,'w') as gf:
            gf.write("Grid type: {0}\n".format(','.join(fields)))
            gf.write("Source files:\n")
            for g in sequence.grids:
                gf.write("    {0}\n".format(g))
        gfdata=open(gridfile).read()
        md5=hashlib.md5()
        md5.update(gfdata)
        gridspec=OrderedDict([
            ('attributes',fields),
            ('type','GeoTIFF'),
            ('filename', gridname),
            ('md5_checksum', md5.hexdigest())
        ])
        grids[gridname]=gridspec

    # Add a component for each time function (should only be one, but could be multiple velocities in theory)

    for nfunc,functime in enumerate(timefuncs):
        components.append((sequence,functime[1],OrderedDict([
            ('description',sequence.description),
            ('bbox',sequence.extent.spec()),
            ('default_uncertainty',component_uncertainty),
            ('spatial_model',gridspec),
            ('time_function',functime[0])
        ])))

basename=defname+'-defmod'
# mdfile=os.path.join(md,'model','metadata.xml')
# mdname=basename+'-metadata.xml'
# with open(mdfile,'rb') as mdf:
#     metadata=mdf.read()
# with open(os.path.join(bd,mdname),'wb') as mdf:
#     mdf.write(metadata)
# md5=hashlib.md5()
# md5.update(metadata)
# mdspec=OrderedDict([('filename',mdname),('md5_checksum',md5.hexdigest())])

# Links to information about the deformation model. 

about=OrderedDict([
    ("href","https://www.linz.govt.nz/nzgd2000"),
    ("rel","about"),
    ("type","text/html"),
    ("title","About the NZGD2000 deformation model")
])
source=OrderedDict([
    ("href","https://www.geodesy.linz.govt.nz/download/nzgd2000_deformation_model"),
    ("rel","source"),
    ("type","application/zip"),
    ("title","Authoritative source of the NZGD2000 deformation model")
])
metadatafunc=lambda v: OrderedDict([
    ("href","https://www.geodesy.linz.govt.nz/download/nzgd2000/metadata/nzgd2000_deformation_{version}_metadata.xml"
        .replace('{version}',v)),
    ("rel","metadata"),
    ("type","application/xml"),
    ("title"," ISO 19115 XML encoded metadata regarding the deformation model")
])
versions=[v for v in m.versions()] if allversions else [m.currentVersion()]
for v in versions:
    vseq=[c for c in components if c[0].versionstart <= v and (c[0].versionend=='0' or c[0].versionend > v)]
    vseq.sort(key=lambda c: c[1])
    if len(vseq) == 0:
        continue
    vcomps=[s[2] for s in vseq]
    vextent=vseq[0][0].extent.copy()
    for s in vseq[1:]:
        vextent.unionWith(s[0].extent)
    modelspec=OrderedDict([
        ('name',m.metadata('model_name')),
        ('version',v),
        ('publication_date',m.versionInfo(v).release_date.strftime(timeformat)),
        ('description',m.metadata('description').replace("\r","")),
        ('authority',OrderedDict([
            ('name',m.metadata('authority')),
            ('url',m.metadata('authority_website')),
            ('address',m.metadata('authority_address')),
            ('email',m.metadata('authority_email'))
        ])),
        ('links',[about,source,metadatafunc(v)]),
        ('model_crs', 'EPSG:4167'),
        ('reference_epoch',refdate.strftime(timeformat)),
        ('bbox', vextent.spec()),
        ('components', vcomps)
    ])
    modeljson=json.dumps(modelspec,indent=2)
    mergefunc=lambda m: m.group(1)+re.sub(r'\s+','',m.group(2))
    modeljson=re.sub(r'(\"extent\"\:\s)((?:[^\]]*\]){3})',mergefunc,modeljson)
    modeljson=re.sub(r'(\"fields\"\:\s)([^\]]*\])',mergefunc,modeljson)
    modeljson=re.sub(r'(\"default_uncertainty\"\:\s)([^\]]*\])',mergefunc,modeljson)
    modeljson=re.sub(r'(\{)(\s*\"epoch\"[^\}]+\"scale_factor\"[^\}]+)',mergefunc,modeljson)
    deffile=basename+'-'+v+'.json'
    with open(os.path.join(bd,deffile),'w') as dfh:
        dfh.write(modeljson)
    print("Created proj deformation master file {0}".format(deffile))
