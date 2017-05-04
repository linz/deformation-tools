#!/usr/bin/python
# Script to build a set of nested grids to meet tolerances etc.

from collections import namedtuple
from shapely.geometry import MultiPolygon,Polygon,LineString,Point,shape
from shapely.prepared import prep
from shapely.ops import unary_union
from shapely import affinity
from shapely.wkt import loads as wkt_load
from subprocess import call
import accuracy_standards
import argparse
import csv
import datetime
import defgrid
import math
import numpy as np
import os
import os.path
import re
import shutil
import numbers
import sys

calc_okada=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'okada','calc_okada')
if not os.path.exists(calc_okada):
    raise RuntimeError('Cannot find calc_okada program at '+calc_okada)

gridtool=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'gridtool','gridtool')
if not os.path.exists(gridtool):
    raise RuntimeError('Cannot find gridtool program at '+gridtool)

# Logging parameters

runtime=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
runtime_version=datetime.datetime.now().strftime("%Y%m%d")
logfile=None
wktfile=None
wktgridfile=None
wktgridlines=False

Unspecified=object()

class Config( object ):
    '''
    Class used to contain configuration information for the process.  This 
    also provides a few static methods for generating information related to
    the configuration, such as filenames.
    '''

    verbose=True
    patchpath=''
    patchroot='patch'
    patchname='patch'
    model=None
        
    # Grids will be calculated using integer representation based origin.
    # Subgrids will be obtained by dividing basesize by 4

    # Origin for calculating grids (integer offsets from this)
    origin=(166.0,-48)

    # Base size for grids, grids are built to this size and smaller
    base_size=(0.15,0.125)

    # Extents over which test grid is built (ie initial search for affected area,
    # and to which grids are trimmed). 

    # Extents based on EEZ
    base_extents=((158,-58),(194,-25))

    # Conversion metres to degrees
    scale_factor=(1.0e-5/math.cos(math.radians(40)),1.0e-5)

    # Factor by which accuracy specs are multiplied to provide grid tolerances
    tolerance_factor=0.4

    # Values for deformation patch generation

    # Test column for determining extents of patch
    base_limit_test_column='err'
    base_limit_bounds_column='ds'

    # Minium value for test column within patch (note - column is in ppm)
    base_limit_tolerance=accuracy_standards.NRF.lap*tolerance_factor*1.0e6
    base_limit_sea_tolerance=accuracy_standards.BGN.lap*tolerance_factor*1.0e6

    # Precision for scaling down from edge of patch
    base_ramp_tolerance=accuracy_standards.NRF.lap*tolerance_factor
    base_ramp_sea_tolerance=accuracy_standards.BGN.lap*tolerance_factor

    # Absolute limit on extents of patch - not required where magnitude (ds)
    # less than this level.
    base_limit_bounds_tolerance=accuracy_standards.NRF.lac*tolerance_factor
    base_limit_sea_bounds_tolerance=accuracy_standards.BGN.na*tolerance_factor

    # Values controlling splitting cells into subgrids

    cell_split_factor=4
    subcell_resolution_tolerance=accuracy_standards.NRF.lac*tolerance_factor
    subcell_resolution_sea_tolerance=accuracy_standards.BGN.na*tolerance_factor

    # Replace entire grid with subcells if more than this percentage 
    # of area requires splitting

    max_subcell_ratio=0.5

    # Maximum split level

    max_split_level=5

    # Hybrid patch grid maximum ppm distortion in forward patch

    forward_patch_test_column='err'
    default_forward_patch_max_distortion=10
    forward_patch_max_level=2

    land_areas=None

# Outputs required ...

    # Shift model component (for updating Landonline)
    build_shift_model=True

    shift_model_header='''
SHIFT_MODEL_MODEL Shift model {name}
FORMAT LINZSHIFT1B
VERSION_NUMBER 1.0
COORDSYS NZGD2000
DESCRIPTION
Model built on {runtime}
Source files: {modelname}
END_DESCRIPTION
'''

    shift_model_component='''
SHIFT_MODEL_COMPONENT {gridfile}
MODEL_TYPE grid
NEGATIVE no
DESCRIPTION
{component}
END_DESCRIPTION
'''

    # LINZ published deformation model

    build_published=True

    published_component_columns='''
        version_added
        version_revoked
        reverse_patch
        component
        priority
        min_lon
        max_lon
        min_lat
        max_lat
        spatial_complete
        min_date
        max_date
        time_complete
        npoints1
        npoints2
        displacement_type
        error_type
        max_displacement
        spatial_model
        time_function
        time0
        factor0
        time1
        factor1
        decay
        file1
        description
        '''.split()

    # Alternative values for creating a Landonline parcel shifting patch
    @staticmethod
    def configureForParcelShift():
        raise RuntimeError('Currently not implemented for hybrid patches')
        Config.base_limit_test_column='ds'
        Config.base_limit_tolerance= 0.05   
        Config.tolerance=accuracy_standards.BGN.lap*tolerance_factor

    # Reduce the accuracies to produce smaller data sets for testing...
    @staticmethod
    def configureForTesting():
        Config.subcell_resolution_tolerance *= 100
        Config.max_split_level -= 2

    # Set patch file name
    @staticmethod
    def setPatchFilenameRoot( patchroot ):
        Config.patchroot=patchroot
        Config.patchpath=os.path.dirname(patchroot)
        Config.patchname=os.path.basename(patchroot)

    @staticmethod
    def setFaultModel( model ):
        Config.model=model

    @staticmethod
    def filename( name=None, extension=None ):
        name = Config.patchname if name is None else name
        extension = extension or ''
        return os.path.join(Config.patchpath,name)+extension

    @staticmethod
    def loadLandAreas( polygon_file ):
        from shapely.wkt import loads
        try:
            with open(polygon_file) as laf:
                wkt=laf.read()
                Config.land_areas=loads(wkt)
            Logger.write( "Using land area definition from "+polygon_file)
            return
        except:
            raise RuntimeError("Cannot load land area definition from "+polygon_file)

    @staticmethod
    def writeTo( log ):
        if callable(log):
            writelog=log
        else:
            writelog=lambda x: log.write(x+"\n")
        writelog("Patch grid calculation parameters:")
        writelog("    Grid origin: {0} {1}".format(*Config.origin))
        writelog("    Base grid cell size: {0} {1}".format(*Config.base_size))
        writelog("    Base extents: {0} {1} to {2} {3}".format(
            Config.base_extents[0][0], Config.base_extents[0][1], 
            Config.base_extents[1][0], Config.base_extents[1][1]))
        writelog("    Approx metres to degrees factor: {0} {1}".format(*Config.scale_factor))

        writelog("    Tolerated distortion outside patch (ppm): land {0} sea {1}".format( 
            Config.base_limit_tolerance, Config.base_limit_sea_tolerance))
        writelog("    Tolerated dislocation outside patch (mm): land {0} sea {1}".format( 
            Config.base_limit_bounds_tolerance*1000, Config.base_limit_sea_bounds_tolerance*1000))
        writelog("    Patch ramp tolerance (ppm): land {0} sea {1}".format(
            Config.base_ramp_tolerance,Config.base_ramp_sea_tolerance))
        writelog("    Tolerated gridding error before splitting (mm)".format(
            Config.subcell_resolution_tolerance*1000, Config.subcell_resolution_sea_tolerance*1000))
        writelog("    Grid cell split factor: {0}".format(Config.cell_split_factor))
        writelog("    Maximum split level: {0}".format(Config.max_split_level))
        writelog("    Maximum percent coverage of level with nested grid: {0}".format(Config.max_subcell_ratio*100)) 
        writelog("")

class Logger( object ):

    logfile=None
    wktfile=None
    wktgridfile=None
    wktgridlines=False
    
    verbose=False

    @staticmethod
    def createLog( logfile=None, wktfile=None, gridwkt=None, showgridlines=None ):
        if logfile:
            Logger.logfile=open(logfile,'w')
        if wktfile:
            Logger.wktfile=open(wktfile,'w')
        if gridwkt:
            Logger.wktgridfile=open(gridwkt,'w')
        if showgridlines is not None:
            Logger.wktgridlines=showgridlines

    @staticmethod
    def setVerbose( verbose=True ):
        Logger.verbose=verbose

    @staticmethod
    def write( message, level=0, prefix="" ):
        prefix='  '*level+prefix
        message="\n".join([prefix+m for m in message.split("\n")])
        if Logger.logfile:
            Logger.logfile.write(message)
            Logger.logfile.write("\n")
            Logger.logfile.flush()
        if Logger.verbose or level < 0:
            sys.stdout.flush()

    @staticmethod
    def dumpGridList( message, gridlist ):
        gridtxt=message+"\n"+"\n".join(["   "+g.shortStr() for g in gridlist])
        Logger.write(gridtxt)

    @staticmethod
    def writeGridDefWkt( name, level, griddef ):
        if not Logger.wktgridfile:
            return
        gridspec=griddef.spec
        name=griddef.name
        level=griddef.level
        ncol=gridspec.ncol
        nrow=gridspec.nrow
        if Logger.wktgridlines:
            for l in spec.lines():
                Logger.wktgridfile.write("{0}|{1}|{2}|{3}|{4}\n".format(
                    name,level,ncol,nrow,l.to_wkt()))
        else:
            Logger.wktgridfile.write("{0}|{1}|{2}|{3}|{4}\n".format(
                name,level,ncol,nrow,spec.boundingPolygonWkt()))


    @staticmethod
    def writeWkt( item, level, wkt ):
        if not Logger.wktfile:
            return
        Logger.wktfile.write("{0}|{1}|{2}\n".format(item,level,wkt))
        Logger.wktfile.flush()


class Util( object ):

    @staticmethod
    def asMultiPolygon( areas ):
        if type( areas) != MultiPolygon:
            if type(areas) == Polygon:
                areas=[areas]
            if type(areas) == list:
                areas=MultiPolygon(areas)
        if not areas.is_valid:
            Logger.write("Attempting to fix invalid geometry in asMultiPolygon")
            areas=areas.buffer(0.0)
            if type(areas) == Polygon:
                areas=MultiPolygon([areas])
        return areas

    @staticmethod
    def bufferedPolygon( polygon, buffer ):
        '''
        Buffer a polygon by an extents in metres.  Approximately scales to metres, applies buffer,
        and scales back. Returns a multipolygon
        '''
        if buffer > 0:
            polygon=affinity.scale(polygon,xfact=1.0/Config.scale_factor[0],yfact=1.0/Config.scale_factor[1],origin=Config.origin)
            polygon=polygon.buffer(buffer)
            polygon=affinity.scale(polygon,xfact=Config.scale_factor[0],yfact=Config.scale_factor[1],origin=Config.origin)
        return Util.asMultiPolygon(polygon)

class GridSpec( object ):
    '''
    Class representing a basic grid definition - min and max coordinates and number of 
    grid cells in each direction.  Note: the number of grid values in each direction is
    one greater than the number of cells.

    Also acts as simple bounding box, with just one grid cell each way
    '''

    @staticmethod
    def fromBounds( bounds ):
        if type(bounds[0]) == float:
            return GridSpec(*bounds)
        return GridSpec(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1] )

    @staticmethod
    def _calc_min_max( base, minv, maxv, csize, multiple ):
        csize *= multiple
        n0 = int(math.floor((minv-base)/csize+0.0001))
        n1 = int(math.ceil((maxv-base)/csize-0.0001))
        return base+n0*csize,base+n1*csize,(n1-n0)*multiple

    def __init__( self, xmin, ymin, xmax, ymax, ngx=1, ngy=1, cellsize=None ):
        if ngx > 100000 or ngy > 100000:
            raise RuntimeError('Odd!')
        if cellsize is not None:
            if ngx == 1:
                ngx=max(int((xmax-xmin)/cellsize[0]+0.5),1)
            if ngy == 1:
                ngy=max(int((ymax-ymin)/cellsize[1]+0.5),1)
        self.xmin=xmin
        self.ymin=ymin
        self.xmax=xmax
        self.ymax=ymax
        self.ngx=ngx
        self.ngy=ngy

    def __str__( self ):
        return 'grid:{0}:{1}:{2}:{3}:{4}:{5}'.format(
            self.xmin, self.ymin,
            self.xmax, self.ymax,
            self.ngx, self.ngy,
            )

    def copy( self ):
        return GridSpec( self.xmin, self.ymin, self.xmax, self.ymax, self.ngx, self.ngy )

    def alignedTo( self, cellsize, origin=None, multiple=1 ):
        '''
        Calculate grid definition (ie min/max lat/lon, number of cells)
        based on grid size (lon,lat cell size), required extents, and 
        an optional integer.  Grids are aligned using integer multiples
        of size from the origin.  
        '''
        if origin is None:
            origin=Config.origin

        xmin,xmax,nln=GridSpec._calc_min_max(origin[0],self.xmin,self.xmax,cellsize[0],multiple)
        ymin,ymax,nlt=GridSpec._calc_min_max(origin[1],self.ymin,self.ymax,cellsize[1],multiple)
        return GridSpec( xmin, ymin, xmax, ymax, nln, nlt )

    def split( self, splitfactor ):
        return GridSpec( self.xmin, self.ymin, self.xmax, self.ymax, 
                        self.ngx*splitfactor, self.ngy*splitfactor )

    def cellSize( self ):
        '''
        Cell size in the XY coordinate system
        '''
        return ((self.xmax-self.xmin)/self.ngx,
                (self.ymax-self.ymin)/self.ngy)

    def cellDimension( self ):
        '''
        Maximum dimension of the cell in metres, converted using
        Config.scale_factor
        '''
        cellsize=self.cellSize()
        return max(cellsize[0]/Config.scale_factor[0],cellsize[1]/Config.scale_factor[1])

    def mergeWith( self, other ):
        pt0 = (min(self.xmin,other.xmin),min(self.ymin,other.ymin))
        pt1 = (max(self.xmax,other.xmax),max(self.ymax,other.ymax))
        dln = ((self.xmax-self.xmin)+(other.xmax-other.xmin))/(self.ngx+other.ngx)
        dlt = ((self.ymax-self.ymin)+(other.ymax-other.ymin))/(self.ngy+other.ngy)
        return GridSpec(pt0[0],pt0[1],pt1[0],pt1[1],cellsize=self.cellSize())

    def expand( self, expandby ):
        '''
        Return extents expanded by an east/north amount.
        expandby is a tuple [de,dn].
        '''
        return GridSpec (
                self.xmin-expandby[0],
                self.ymin-expandby[1],
                self.xmax+expandby[0],
                self.ymax+expandby[1],
                cellsize=self.cellSize())

    def intersectWith( self, other ):
        '''
        Generate bounded extents around a polygon.
        buffer is buffer extents in metres
        cellbuffer is (lon,lat) buffer size 
        trim is region to which buffered extents are trimmed
        '''
        return GridSpec (
                max(self.xmin,other.xmin),
                max(self.ymin,other.ymin),
                min(self.xmax,other.xmax),
                min(self.ymax,other.ymax),
                cellsize=self.cellSize())

    def overlaps( self, other ):
        return not (
            self.xmin > other.xmax or
            other.xmin > self.xmax or
            self.ymin > other.ymax or
            other.ymin > self.ymax)

    def boundingPolygon( self ):
        bl=(self.xmin,self.ymin)
        br=(self.xmin,self.ymax)
        tr=(self.xmax,self.ymax)
        tl=(self.xmax,self.ymin)
        return Polygon((bl,br,tr,tl,bl))

    def boundingPolygonWkt( self ):
        return self.boundingPolygon().to_wkt()

    def _lonLatValues( self ):
        lonv = list(self.xmin+(self.xmax-self.xmin)*float(x)/self.ngx
                for x in range(self.ngx+1))
        latv = list(self.ymin+(self.ymax-self.ymin)*float(x)/self.ngy
                for x in range(self.ngy+1))
        return lonv,latv

    def nodes( self ):
        '''
        Iterator returning the nodes of the grid
        '''
        lonv,latv=self._lonLatValues()
        for lt in latv:
            for ln in lonv:
                yield (ln,lt)

    def lines( self ):
        '''
        Iterator returning lines (as coordinate pairs) forming the internal and
        external edges of the grid.
        '''
        lonv,latv=self._lonLatValues()
        for lat in latv:
            for ln0,ln1 in zip(lonv[:-1],lonv[1:]):
                yield LineString((ln0,lat),(ln1,lat))
        for lon in lonv:
            for lt0,lt1 in zip(latv[:-1],latv[1:]):
                yield LineString((lon,latv[i]),(lon,latv[i+1]))

    def area( self ):
        return (self.xmax-self.xmin)*(self.ymax-self.ymin)

    def mergeIntoList( self, deflist ):
        merged=self.copy()
        while True:
            overlap = None
            for def1 in deflist:
                if self.overlaps(def1):
                    overlap=def1
                    break
            if overlap:
                deflist.remove(overlap)
                merged=self.mergeWith(overlap)
            else:
                break
        deflist.append(merged)

class FaultModel( object ):

    filename_templates=('{0}','{0}.model','fault_models/{0}','fault_models/{0}.model')

    class PatchVersion( object ) :

        TimePoint=namedtuple('TimePoint','date factor')

        def __init__( self, version, event, date, reverse=False, time_model=None, hybrid=False ):
            self.version=version
            self.event=event
            self.date=date
            self.reverse=reverse
            self.hybrid=hybrid
            if time_model is None:
                time_model=[FaultModel.PatchVersion.TimePoint(date,1.0)]
            self.time_model=time_model

        
        def __str__( self ):
            return "\n".join((
                "Version: {0}".format(self.version),
                "Event: {0}".format(self.event),
                "Date: {0}".format(self.date.strftime('%Y-%m-%d')),
                "Reverse: {0}".format(self.reverse),
                "Hybrid: {0}".format(self.hybrid),
                "TimeModel: {0}\n  ".format(
                    "\n  ".join(("{0} {1:.2f}".format(x[0].strftime('%Y-%m-%d'),x[1]) 
                                 for x in self.time_model)))
                ))

    def __init__( self, modelfile ):
        self.patchversions=[]
        self.modelpath=None
        self.modelname=None
        self.modeldate=None
        self.description=None
        self.hybridTolerance=None
        for template in FaultModel.filename_templates:
            mf=template.format(modelfile)
            if os.path.exists(mf):
                self.modelpath=mf
                break
        if self.modelpath is None:
            raise RuntimeError('Cannot find fault model file {0}'.format(modelfile))

        self.loadModelFile()
        self.loadPatchSpecs()


    def loadModelFile( self ):
        '''
        Loads a patch model file, and returns a list of patch versions associated with it.
        '''
        
        self.modelname=os.path.splitext(os.path.basename(self.modelpath))[0]
        with open(self.modelpath) as f:
            event = f.readline()
            model = f.readline()
            version = f.readline()
            description = event+' '+model+' '+version
            match=re.search(r'(\d{1,2})\s+(\w{3})\w*\s+(\d{4})\s*$',event)
            modeldate = None
            if match:
                try:
                    datestr=match.group(1)+' '+match.group(2)+' '+match.group(3)
                    modeldate = datetime.datetime.strptime(datestr,'%d %b %Y')
                except:
                    pass
            if not modeldate:
                raise RuntimeError("Event record at start of {0} does not end with a valid date".format(modelfile))
            self.description=description
            self.modeldate=modeldate.strftime('%Y-%m-%d')

    def _parseDate( self, datestr ):
        datestr=datestr.strip()
        for format in ('%d %b %Y','%d %B %Y','%Y-%m-%d'):
            try:
                return datetime.datetime.strptime(datestr,format)
            except:
                pass
        raise ValueError("Cannot parse date {0}".format(datestr))

    def _buildTimeModel(self,time_model,modeldate):
        if time_model is None:
            return [FaultModel.PatchVersion.TimePoint(self._parseDate(modeldate),1.0)]
        model=[]
        ntime=len(time_model)
        if ntime < 2:
            raise RuntimeError("Invalid or empty TimeModel")
        if ntime % 2 != 0:
            raise RuntimeError("TimeModel must consist of paired date and value")
        for i in range(0,ntime,2):
            fdate=self._parseDate(time_model[i])
            try:
                factor=float(time_model[i+1])
            except:
                raise RuntimeError("Invalid scale factor {0} in TimeModel"
                                   .format(time_model[i+1]))
            model.append(FaultModel.PatchVersion.TimePoint(fdate,factor))

        for i in range(len(model)-1):
            if model[i].date > model[i+1].date:
                raise RuntimeError("Dates out of order in TimeModel")
        return model

    def loadPatchSpecs( self ):
        patchfile=os.path.splitext(self.modelpath)[0]+'.patch'
        if not os.path.exists(patchfile):
            raise RuntimeError('Cannot find fault model patch version file {0}'.format(patchfile))

        patchversions=[]
        version=None
        reverse=False
        hybrid=False
        hybrid_tol=None
        nested=None
        time_model=None
        in_model=False

        with open(patchfile) as f:
            for l in f:
                # Ignore comments and blank lines
                if re.match(r'^\s*(\#|$)',l):
                    continue
                # split out command and value...
                cmdmatch=re.match(r'^\s*(\w+)\s*\:\s*(.*?)\s*$',l)
                if cmdmatch:
                    in_model=False
                    command=cmdmatch.group(1).lower()
                    value=cmdmatch.group(2)
                    if command == 'version':
                        if version is not None:
                            if value <= version:
                                raise RuntimeError("Patch versions out of order in {0}: {1} <= {2}"
                                                   .format(patchfile,value,version))
                            try:
                                patchversions.append(FaultModel.PatchVersion(
                                    version,
                                    event,
                                    self._parseDate(modeldate),
                                    reverse=reverse,
                                    hybrid=hybrid,
                                    time_model=self._buildTimeModel(time_model,modeldate)))
                            except Exception as ex:
                                raise RuntimeError("Error in {0}: {1}"
                                                   .format(patchfile,ex.message))
                        version=value
                        if not re.match(r'^[12]\d\d\d[01]\d[0123]\d$',version):
                            raise RuntimeError("Invalid deformation model version {0} in {1}"
                                               .format(version,patchfile))
                        continue
                    if not version:
                        raise RuntimeError("Version must be the first item in patch file {0}"
                                           .format(patchfile))
                    if command == 'event':
                        event=value
                    elif command == 'date':
                        modeldate=value
                    elif command == 'type':
                        if value.lower() == 'forward':
                            reverse=False
                        elif value.lower() == 'reverse':
                            reverse=True
                        elif value.lower() == 'hybrid':
                            reverse=True
                            hybrid=True
                        else:
                            raise RuntimeError("Invalid Type - must be forward or reverse in {0}"
                                               .format(patchfile))
                    elif command == 'forwardpatchmaxdistortionppm':
                        fvalue=float(value)
                        if hybrid_tol is None:
                            hybrid_tol=fvalue
                        elif hybrid_tol != falue:
                            raise RuntimeError("Inconsistent forwardPatchMaxDistortionPpm {0} in {1}"
                                               .format(value,patchfile))
                    elif command == 'subgridmethod':
                        if value.lower() == 'nested':
                            vnested=True
                        elif value.lower() == 'indpendent':
                            vnested=False
                        else:
                            raise RuntimeError("Invalid SubgridMethod - must be nested or independent in {0}"
                                               .format(patchfile))
                        if nested is None:
                            nested=vnested
                        elif nested != vnested:
                            raise RuntimeError("Inconsistent SubgridMethod - all patch versions must use the same method in {0}"
                                               .format(patchfile))

                    elif command == 'timemodel':
                        time_model=value.split()
                        in_model=True
                    else:
                        raise RuntimeError("Unrecognized command {0} in {1}".format(strip(),patchfile))

                elif in_model: # Not a command line
                    time_model.extend(l.split())
                else:
                    raise RuntimeError("Invalid data \"{0}\" in {1}".format(l,patchfile))
                    
            if version is None:
                raise RuntimeError("No version specified in {0}".format(patchfile))

            try:
                patchversions.append(FaultModel.PatchVersion(
                    version,
                    event,
                    self._parseDate(modeldate),
                    reverse=reverse,
                    hybrid=hybrid,
                    time_model=self._buildTimeModel(time_model,modeldate)))
            except Exception as ex:
                raise RuntimeError("Error in {0}: {1}"
                                   .format(patchfile,ex.message))

            for i in range(len(patchversions)-1):
                if patchversions[i].version >= patchversions[i+1].version:
                    raise RuntimeError("Deformation model versions out of order in {0}".format(patchfile))

            self.patchversions=patchversions
            self.hybrid=hybrid
            self.hybrid_tol=hybrid_tol

    def _cacheKey( self, xy ):
        key="{0:.10f}\t{1:.10f}".format(float(xy[0]),float(xy[1]))
        key=re.sub(r'(\.\d+)0*(\t|$)',r'\1\2',key)
        return key

    def _cacheFile( self ):
        return Config.filename(name=self.modelname,extension='_grid.cache')

    def _loadGridFromCache( self, spec, gridfile ):
        keys={}
        for i,n in enumerate(spec.nodes()):
            keys[self._cacheKey(n)]=i
        values=[None]*len(keys)

        keyre=re.compile('^[^\t]+\t[^\t]+')
        header=None
        if os.path.exists(self._cacheFile()):
            with open(self._cacheFile()) as cf:
                header=cf.next()
                for l in cf:
                    match=keyre.match(l)
                    if match:
                        k=match.group()
                        if k in keys:
                            values[keys[k]]=l
        
        if os.path.exists(gridfile):
            os.remove(gridfile)

        missing=[k for k in keys if values[keys[k]] is None]
        if missing:
            return missing

        with open(gridfile,'w') as gf:
            gf.write(header)
            for v in values:
                gf.write(v)

    def _addToCache( self, missing ):
        global calc_okada
        nmissing=len(missing)
        cachefile=self._cacheFile()
        missingfile=Config.filename(name=self.modelname,extension='_grid.cache.missing.tmp')
        calcfile=Config.filename(name=self.modelname,extension='_grid.cache.calcs.tmp')

        blocksize=1000
        ncalc=0
        for imissing in range(0,nmissing,blocksize):
            sys.stdout.write("          ... calculated {0}/{1} grid values\r".format(ncalc,nmissing))
            sys.stdout.flush()
            with open(missingfile,'w') as mf:
                for l in missing[imissing:imissing+blocksize]:
                    mf.write(l)
                    mf.write('\n')

            modelpath=self.modelpath
            params=[calc_okada,'-f','-x','-l','-s',modelpath,missingfile,calcfile]
            call(params)

            header=None
            if os.path.exists(cachefile):
                with open(cachefile) as cf:
                    header=cf.next()

            with open(calcfile) as calcf:
                calcheader=calcf.next()
                mode='w' if calcheader != header else 'a'
                with open(cachefile,mode) as cf:
                    if mode == 'w':
                        cf.write(calcheader)
                    for l in calcf:
                        parts=l.split('\t')
                        key=self._cacheKey(parts[:2])
                        cf.write(key)
                        cf.write('\t')
                        cf.write('\t'.join(parts[2:]))
                        ncalc += 1
            os.remove(missingfile)
            os.remove(calcfile)
        if nmissing:
            sys.stdout.write("\n")

        
    def calcGrid( self, gridspec, gridfile ):
        global calc_okada
        gridspecstr=str(gridspec)
        modelpath=self.modelpath
        # if verbose:
        #     print "Calculating deformation on {0}".format(gridfile)
        #     print "Model {0}".format(modeldef)
        #     print "Grid spec {0}".format(gridspec)

        params=[calc_okada,'-f','-x','-l','-s',modelpath,gridspecstr,gridfile]
        meta='\n'.join(params)
        metafile=gridfile+'.metadata'
        built = False
        # Check if grid file is already built - mainly convenience for
        # script development as this takes much longer than anything else!
        if os.path.exists(gridfile) and os.path.exists(metafile):
            with open(metafile) as f:
                oldmeta=f.read()
                if oldmeta.strip() == meta.strip():
                    built=True
            if os.path.getmtime(gridfile) > os.path.getmtime(metafile):
                built=False
            if built:
                if os.path.getmtime(modelpath) > os.path.getmtime(gridfile):
                    built=False
        if not built:
            if os.path.exists(gridfile):
                os.remove(gridfile)
            print "          Building grid {0} ({1} x {2} points) ...".format(
                gridfile,gridspec.ngx,gridspec.ngy)
            missing=self._loadGridFromCache(gridspec,gridfile)
            if missing:
                self._addToCache(missing)
                missing=self._loadGridFromCache(gridspec,gridfile)
                if missing:
                    raise RuntimeError("Failed to calculate all required grid values {0}\n{1}".
                                   format(gridfile," ".join(params)))

            # call(params)
            if not os.path.exists(gridfile):
                raise RuntimeError("Failed to calculate grid file {0}\n{1}".
                                   format(gridfile," ".join(params)))
            with open(metafile,'w') as f:
                f.write(meta)
        else:
            print "          Reloading existing grid {0} ...".format(gridfile)
        return gridfile

class PatchGridDef:
                                 
    def __init__( self, name=None, subpatch=None, extents=None, buffered_extents=None, level=1, spec=None, parent=None, source=None ):
        self.level=level
        self.parent=parent
        self.name=name
        self.subpatch=subpatch
        self._extents=extents
        self._extentsWkt=None
        self._bufferedExtents=buffered_extents or extents
        self._bufferedExtentsWkt=None
        self.spec=spec
        self._grid=None
        self._strainGrid=None
        self._csv=None
        self.gdf=None
        self.base=None
        self.basename=None
        self.isforward=None  # Component is forward/reverse in hybrid grid else None (= either)
        self._buffered_extents=None 
        self.source=source
        if spec is not None:
            Logger.writeWkt('PatchGridDef {0} spec {1}'.format(name,str(spec)),
                            level,spec.boundingPolygonWkt())
        if extents is not None:
            Logger.writeWkt('PatchGridDef {0} extents'.format(name,subpatch),
                            level,extents.to_wkt())
            Logger.writeWkt('PatchGridDef {0} buffered extents'.format(name,subpatch),
                            level,self._bufferedExtents.to_wkt())
            if not self._extents.is_valid:
                raise RuntimeError('Invalid extents defined for grid {0}'.format(name))

    def updatedWith( self, 
                    name=Unspecified, 
                    subpatch=Unspecified, 
                    extents=Unspecified, 
                    buffered_extents=Unspecified, 
                    level=Unspecified, 
                    parent=Unspecified, 
                    isforward=Unspecified, 
                    source=Unspecified ):
        pgd=PatchGridDef(
            name=name if name is not Unspecified else self.name,
            subpatch=subpatch if subpatch is not Unspecified else self.subpatch,
            extents=extents if extents is not Unspecified else self.extents,
            buffered_extents=buffered_extents if buffered_extents is not Unspecified else self.bufferedExtents,
            level=level if level is not Unspecified else self.level,
            source=source if source is not Unspecified else self.source,
            parent=parent if parent is not Unspecified else None
            )
        pgd.isforward=isforward if isforward is not Unspecified else self.isforward
        pgd.base=self.base
        pgd.basename=self.basename
        return pgd

    def setSource( self, source ):
        self.source=source
        self.grid=None

    @property
    def filename( self ):
        return self._fileWithExtension('.grid')

    def _loadGrid( self, build_only=False ):
        if self.source is None:
            raise RuntimeError('PatchGridDef._loadGrid: Grid source not defined')
        gridfile=self.filename
        source=self.source
        if callable( source):
            source(self.spec, gridfile )
        elif os.path.exists( source ):
            if( source != gridfile ):
                shutil.copyfile(source,gridfile)
                try:
                    shutil.copystat(source,gridfile)
                    self.source=gridfile
                    self.grid=None
                except:
                    pass
        else:
            raise RuntimeError('PatchGridDef._loadGrid: Invalid grid source {0}'
                               .format(source))
        if not os.path.exists(gridfile):
            raise RuntimeError('PatchGridDef._loadGrid: Grid file {0} is not built'
                               .format(gridfile))

        if build_only:
            return None
        grid=defgrid.DeformationGrid(gridfile)
        # Should we reload the grid here
        return grid

    @property
    def builtFilename( self ):
        filename=self.filename
        if not os.path.exists(filename):
            self._loadGrid(build_only=True)
        return filename

    @property
    def grid( self ):
        if self.spec is None:
            raise RuntimeError('{0} grid spec not calculated before calculating grid values'
                               .format(self.name))
        if self._grid is None:
            self._grid=self._loadGrid()
        return self._grid

    @property
    def csvfile( self ):
        global gridtool
        if self._csv is None:
            gridfile=self.builtFilename
            csvfile = re.sub(r'(\.grid)?$','.csv',gridfile,re.I)
            self.writeCsvFile( csvfile, -1 )
            self._csv=csvfile
        return self._csv

    def writeCsvFile( self, filename, ndp=4 ):
        global gridtool
        if self._csv is None:
            gridfile=self.builtFilename
            csvfile = re.sub(r'(\.grid)?$','.csv',gridfile,re.I)
            commands=[gridtool,
                      'read','maxcols','3',gridfile,
                      'write','ndp',str(ndp),'csv','dos',filename]
            call(commands)
        return filename

    @property
    def gdffile( self ):
        global gridtool
        if self._gdf is None:
            gridfile=self.builtFilename
            gdffile = re.sub(r'(\.grid)?$','.gdf',gridfile,re.I)
            header2='Model '+Config.model.description
            header1='Patch component '+self.name
            header3=runtime
            commands=[gridtool,
                      'read','maxcols','3',gridfile,
                      'write_linzgrid','NZGD2000',
                      header1,
                      header2,
                      header3,
                      gdffile]
            call(commands)
            self._gdf=gdffile
        return self._gdf

    def _fileWithExtension( self, extension ):
        return Config.filename(name=self.name,extension=extension)

    @property
    def extents( self ):
        '''
        Extents over which the grid is required to ensure that the gridded model
        matches the physical model to the required accuracy
        '''
        return self._extents

    @property
    def extentsWktFile( self ):
        if self._extentsWkt is None and self.extents is not None:
            extentfile=self._fileWithExtension('extent.wkt')
            with open(extentfile,'w') as wktf:
                for pgn in self.extents:
                    wktf.write(pgn.to_wkt())
                    wktf.write('\n')
            self._extentsWkt=extentfile
        return self._extentsWkt

    @property
    def bufferedExtents( self ):
        '''
        Return the buffered extents over which the grid may differ from its
        parent.

        Returned as a single MultiPolygon
        '''
        if self._bufferedExtents is None and self.extents is not None:
            self._bufferedExtents=self.extents
        return self._bufferedExtents

    @property
    def bufferedExtentsWktFile( self ):
        '''
        Return the buffered extents over which the grid may differ from its
        parent.

        Returned as a single MultiPolygon
        '''
        if self._bufferedExtentsWkt is None and self.bufferedExtents is not None:
            extentfile=self._fileWithExtension('buffer.wkt')
            with open(extentfile,'w') as wktf:
                for pgn in self.extents:
                    wktf.write(pgn.to_wkt())
                    wktf.write('\n')
            self._bufferedExtentsWkt=extentfile
        return self._bufferedExtentsWkt

    def mergeIntoParent( self ):
        '''
        Modify the grid to match transition to the parent grid between
        its extents and its buffered extents.  

        Grid is not modified inside its extents (ie where it is distinctly different).
        Merge is between extents and buffered extents.  Outside buffered extents and 
        on boundary it is always matched to the parent.
        '''
        global gridtool
        if self.parent is None:
            return
        gridfile=self.builtFilename
        extwkt=self.extentsWktFile
        bufwkt=self.bufferedExtentsWktFile
        pgridfile=self.parent.builtFilename
        commands=[gridtool,
                  'read','maxcols','3',gridfile,
                  'trimto','buffer','1','wkt',bufwkt,
                  'trimto',pgridfile,
                  'zero','outside',extwkt,'edge','1',
                  'add','maxcols','3',pgridfile,'outside',extwkt,'edge','1',
                  'smooth','linear','outside',extwkt,'not','outside',bufwkt,'not','edge','1',
                  'write',gridfile
                 ]
        call(commands)
        self.setSource(gridfile)

    def maxShiftOutsideAreas( self, areas ):
        grid=self.grid
        total_area = prep(areas.buffer(0.0001))
        dscol = grid.getIndex('ds')
        return max((n[dscol] for n in grid.nodes() 
                     if not total_area.contains(Point(n[0],n[1]))))

    def regionsExceedingLevel( self, column, tolerance ):
        grid=self.grid
        if column not in grid.columns:
            if self._strainGrid is None:
                self._strainGrid=grid.strainComponents()
            grid=self._strainGrid
        return MultiPolygon(grid.regionsExceedingLevel( column, tolerance ))

    def areasNeedingFinerGrid( self, tolerance, celldimension  ):
        '''
        Determine the areas needing a finer grid than cellsize based
        on a required tolerance
        '''
        grid=self.grid
        gridr = grid.calcResolution(tolerance,precision=0.0000001)
        subcell_areas=gridr.regionsLessThanLevel('reqsize',celldimension )
        if len(subcell_areas) <= 0:
            return None
        return subcell_areas

    # Recursive function (recursing on grid level) that populates the array 
    # gridlist with PatchGridDef objects.

    def splitToSubgrids( self ):
        griddef=self
        subgrids=[self]
        Logger.write("Building level {0} patch {1}".format(self.level,self.name),self.level,self.name)

        subcell_areas=[]
        if self.level < Config.max_split_level:
            # Try calculating subcells if necessary

            # Form a trial grid split the cellsize for the next level
            splitspec=self.spec.split(Config.cell_split_factor)
            name="{0}_SL{1}".format(self.name,self.level+1)
            level1=self.level+1
            grid2=PatchGridDef(
                name=name,
                level=level1,
                spec=splitspec,
                source=self.source
                )

            # Determine required resolution of subcells, and find the extents for which 
            # the subcell resolution is required (ie the required resolution is less than
            # the parent cell size)

            resolution=self.spec.cellDimension()
            subcell_areas=grid2.areasNeedingFinerGrid( 
                Config.subcell_resolution_tolerance, resolution)
            if subcell_areas is not None:
                subcell_areas=Util.asMultiPolygon(subcell_areas)
                Logger.writeWkt("{0} subcells".format(name),level1,subcell_areas.to_wkt())
                
                if Config.land_areas is not None:
                    subcell_areas=subcell_areas.intersection(Config.land_areas)
                    Logger.writeWkt("{0} subcells on land area".format(name),
                                    level1,subcell_areas.to_wkt())

                    # If limiting to land areas, then also select areas where sea resolution
                    # is to be applied.
                
                    sea_areas= grid2.areasNeedingFinerGrid( 
                        Config.subcell_resolution_sea_tolerance,
                        resolution)
                    if sea_areas is not None:
                        sea_areas=Util.asMultiPolygon(sea_areas)
                        Logger.writeWkt("{0} sea subcells".format(name),
                                        level1,sea_areas.to_wkt())
                        subcell_areas=subcell_areas.union(Util.asMultiPolygon(sea_areas))
                        Logger.writeWkt("{0} combined subcells".format(name),
                                        level1,subcell_areas.to_wkt())

                # Note: previous version of code provided for replacing parent 
                # grid where subgrids covered significant portion of parent
                #
                # This will be replaced with an optimize grids function 
                # when generating the published grids

                subcells=PatchGridList.gridsForAreas( 
                    subcell_areas, 
                    resolution,
                    splitspec.cellSize(),
                    keepWithin=self.spec,
                    level=self.level+1,
                    subpatch=self.subpatch,
                    parent=self,
                    source=self.source
                    )

                # Now split each subgrid and add to the total list of grids

                for s in subcells:
                    subgrids.extend(s.splitToSubgrids())
            
        return subgrids

    def isTopLevelGrid( self ):
        return self.parent is None

    def topLevelGrid( self ):
        return self if self.isTopLevelGrid() else self.parent.topLevelGrid()

    def owns( self, other ):
        ''' 
        True if other belongs to this grid (as child or if it is this grid
        '''
        if other == self:
            return True
        elif other is None:
            return False
        else:
            return self.owns(other.parent)

    def family( self, gridlist ):
        family=[g for g in gridlist if self.owns(g)]
        family.sort( key=lambda x: (x.level,x.file) )
        return family

    def shortStr( self ):
        return "Grid: {0} parent={1} level={2} F/R={3}".format(
            self.name,
            self.parent.name if self.parent is not None else '-',
            self.level,
            'F' if self.isforward is True else
                'R' if self.isforward is False else
                '-'
            )

    def __str__( self ):
        return "\n".join((
            "Grid Name: {0}".format(self.name),
            "    Subpatch: {0}".format(self.subpatch),
            "    Grid spec: {0}".format(self.spec),
            "    Parent grid: {0}".format(self.parent.name if self.parent else ""),
            "    Grid file: {0}".format(self.filename),
            "    Extent file: {0}".format(self.extentsWktFile),
            "    Level: {0}".format(self.level),
            "    Csv: {0}".format(self._csv),
            "    Gdf: {0}".format(self.gdf),
            "    Base: {0}".format(self.base),
            "    Base name: {0}".format(self.basename),
            "    Is Forward?: {0}".format(self.isforward),
            ))

class PatchGridList( list ):

    @staticmethod
    def gridsForAreas( areas, buffered_areas, cellsize, keepWithin=None, level=1, 
                      name=None,subpatch='',parent=None,source=None):
        '''
        Calculate a set of grids to cover the polygon areas defined in areas.
        The areas are expanded by cellsize to provide a buffer for smoothing on 
        to parent grid or zero. 
        The areas are then allocated to rectangular grids of the specified 
        cellsize and aligned with the origin.

        When the areas are expanded to grids they can end up overlapping, in
        which case they are merged. 
        
        NOTE that it could be possible to split areas up to provide a more 
        efficient coverage.  For example narrow areas running to NE or NW could 
        possibly be split into multiple offset adjacent grids.  This is better done 
        as a separate postprocessing step to optimise total grid size.
        '''

        # pname used for debugging info
        pname=name
        if pname is None:
            pname=Config.patchname
            if subpatch:
                pname=pname+"_"+subpatch+'??'
            pname=pname+"_L{0}".format(level)

        areas=Util.asMultiPolygon(areas)
        Logger.writeWkt('Required extents for {0}/{1}'.format(pname,subpatch),
                        level,areas.to_wkt())
        if not areas.is_valid:
            raise RuntimeError('Invalid area defined for splitting {0}/{1}'.
                               format(pname,subpatch))

        if isinstance(buffered_areas,numbers.Number):
            buffered_areas=Util.bufferedPolygon( areas, buffered_areas )
        buffered_areas=Util.asMultiPolygon(buffered_areas)
        Logger.writeWkt('Buffered extents for {0}/{1}'.format(pname,subpatch),
                        level,buffered_areas.to_wkt())

        merged_specs=[]
        for a in buffered_areas:
            if not a.intersects(areas):
                continue
            bounds = GridSpec.fromBounds(a.bounds).expand(cellsize)
            if keepWithin is not None:
                bounds=bounds.intersectWith(keepWithin)
            bounds=bounds.alignedTo(cellsize)
            bounds.mergeIntoList( merged_specs )

        merged_specs.sort(key=lambda x: x.area(), reverse=True)

        if name is None:
            if len(merged_specs) > 1:
                subpatch=(subpatch or '')+"P{0}"
            name=Config.patchname
            if subpatch:
                name=name+"_"+subpatch
            name=name+"_L{0}".format(level)

        patchdefs=PatchGridList()
        for i,spec in enumerate(merged_specs):
            extents=MultiPolygon([a for a in areas if a.intersects(spec.boundingPolygon())])
            buffered_extents=MultiPolygon([a for a in buffered_areas if a.intersects(spec.boundingPolygon())])
            patchname=name.format(i)
            gsubpatch=subpatch.format(i)
            patchdefs.append(PatchGridDef(name=patchname,subpatch=gsubpatch,spec=spec,
                                          extents=extents,buffered_extents=buffered_extents,
                                          level=level,parent=parent,source=source))

        return patchdefs


def split_forward_reverse_patches( trialgrid, hybrid_tol, gridlist ):
    Logger.write("Splitting patch grids into forward/reverse patches")
    Logger.write("Forward patch tolerated distortion (ppm): {0}".format(hybrid_tol))
    reverse_extents = trialgrid.regionsExceedingLevel(Config.forward_patch_test_column,hybrid_tol)
    Logger.write("{0} regions found".format(len(reverse_extents)))
    Logger.writeWkt("Extents requiring reverse patch (base tolerance)",0,reverse_extents.to_wkt())
    reverse_grid_extents=[g.bufferedExtents
                          for g in gridlist if g.level > Config.forward_patch_max_level]
    if reverse_grid_extents:
        reverse_grid_extents.append(reverse_extents)
        reverse_extents=unary_union(reverse_grid_extents)
        if type(reverse_extents) == Polygon:
            reverse_extents=MultiPolygon([reverse_extents])

    reversewkt=Config.filename(extension='_reverse_patch_extents.wkt')
    with open(reversewkt,'w') as f:
        f.write(reverse_extents.to_wkt())

    #Logger.write("Grid list:")
    #for g in gridlist:
    #    Logger.write("{0}".format(g),prefix="    ")

    # Find grids are fully forward patches, and calculate forward
    # patches for those which are not.
    # 
    # For patches up the forward_patch_max_level the significant part 
    # of the patch (defined by extents) is either completely outside the 
    # reverse patch area, in which case the grid gets directly copied,
    # or it overlaps it, in which case a forward patch and reverse patch
    # is constructed from it.

    not_reverse=set()

    # Forward grid that applies for each grid.
    forward_grids={}
    reverse_grids={}
    hybrid_patches={}
    
    splitgrids=PatchGridList()
    
    for g in sorted(gridlist, key=lambda x: x.level):

        # Identify grids which do not require a reverse patch

        #print "Splitting",g.name,g.extents.is_valid,g.buffer,g.bufferedExtents.is_valid
        #print reverse_extents.is_valid

        if (g.parent in not_reverse 
                or not g.bufferedExtents.intersects(reverse_extents)):
            if g.level > Config.forward_patch_max_level:
                raise RuntimeError('Forward patch required for level greater than {0}'
                                   .format(Config.forward_patch_max_level))
            Logger.write('Copying {0} as forward patch grid only'.format(g.name))
            not_reverse.add(g)
            fgridname=g.name+'_F'
            gf=g.updatedWith(isforward=True,name=fgridname,parent=forward_grids.get(g.parent,None),
                             source=g.builtFilename)
            gf.mergeIntoParent()
            forward_grids[g]=gf
            splitgrids.append(gf)
            continue

        # Create smoothed grid for forward patch, smoothing over reverse patch
        # extents
        #
        # Have tried different smoothing options (including datumgrid), but
        # none seems to offer particular advantage.  Could be an area for more
        # research

        if g.level <= Config.forward_patch_max_level:
            fgridname=g.name+'_F'
            basefg=forward_grids.get(g.parent,None)

            gf=g.updatedWith(name=fgridname,isforward=True,source=None,
                            parent=basefg)
            fgridfile=gf.filename

            Logger.write('Creating smoothed forward patch grid for {0}'.format(g.name))
            commands=[
                gridtool,
                'read','maxcols','3',g.builtFilename
                #,'write','csv',g._fileWithExtension('_presmooth.csv')
                ]
            if basefg is not None:
                commands.extend(['replace','maxcols','3',basefg.builtFilename,'where','inside',reversewkt])
            commands.extend(['smooth','linear','inside',reversewkt])
            if basefg is not None:
                commands.extend(['not','on_grid',basefg.builtFilename])
            commands.extend(['write',fgridfile])
            #commands.extend(['write','csv',g._fileWithExtension('_postsmooth.csv')])
        
            call(commands)

            if not os.path.exists(fgridfile):
                raise RuntimeError('Cannot build smoothed forward patch '+fgridfile)
            gf.setSource(fgridfile)
            Logger.write('Created forward patch for {0}'.format(fgridfile))

            gf.mergeIntoParent()
            forward_grids[g]=gf
            splitgrids.append(gf)
        elif g.parent is not None:
            # If not creating a forward patch, then use the forward patch
            # associated with parent grid
            forward_grids[g]=forward_grids[g.parent]
        else:
            raise RuntimeError('Split grid {0} has no parent'.format(g.name))


        # Now create the reverse patch files by subtracting the forward patch
        # from them

        rgridname=g.name+'_R'
        rgrid=g.updatedWith(name=rgridname,parent=reverse_grids.get(g.parent,None),isforward=False)
        rgridfile=rgrid.filename
        fgrid=forward_grids[g]
        Logger.write('Creating reverse patch grid for {0}'.format(g.name))
        commands=[
            gridtool,
            'read','maxcols','3',g.builtFilename,
            'subtract',fgrid.builtFilename,
            'trim','1',
            'write',rgridfile]
        call(commands)
        if not os.path.exists(rgridfile):
            raise RuntimeError('Cannot build reverse patch '+rgridfile)
        rgrid.setSource(rgridfile)
        Logger.write('Created reverse patch for {0}'.format(fgridfile))
        rbase=reverse_grids.get(g,None)
        reverse_grids[g]=rgrid
        splitgrids.append(rgrid)

    Logger.write("\nSplit grids")
    for g in splitgrids:
        Logger.write("-------------------\n{0}".format(g))
    
    return splitgrids

def build_deformation_grids( splitbase=True ):
    '''
    Function builds all the grids that are used.  Basically sets up a trial grid to determine
    the extents on which patches are required, then for each patch required within the extents
    calls create_grids to build the actual patches.  create_grids is a recursive function that 
    will break the grids down further if the grid resolution compromises the accuracy requirements

    The returned list of grid definitions includes a field called base.  If splitbase is true
    then this will identify the separate patch areas, otherwise it will be a single patch area.
    '''

    # Calculate a trial grid to determine extents of patches

    gridsource=Config.model.calcGrid

    basespec=GridSpec.fromBounds(Config.base_extents)
    trialgriddef=basespec.alignedTo(Config.base_size)
    trialgridfile=Config.patchname+"_trial"

    Logger.write("Building trial grid using {0} to {1}".format(
        trialgriddef,trialgridfile))
    Logger.writeWkt("Trial extents",0,trialgriddef.boundingPolygonWkt())

    trialgrid=PatchGridDef(trialgridfile,spec=trialgriddef,
                          source=gridsource)


    # Determine the extents on which the base tolerance is exceed - returns
    # a list of polygons

    real_extents = trialgrid.regionsExceedingLevel(Config.base_limit_test_column,Config.base_limit_tolerance)
    Logger.write("{0} patch extents found".format(len(real_extents)))
    Logger.writeWkt("Extents requiring patch (base tolerance)",0,real_extents.to_wkt())

    # Bounds beyond which patch not required because guaranteed by local accuracy bounds
    bounds_extents = trialgrid.regionsExceedingLevel(
        Config.base_limit_bounds_column,Config.base_limit_bounds_tolerance)
    Logger.writeWkt("Absolute bounds on patch",0,bounds_extents.to_wkt())

    real_extents=real_extents.intersection( bounds_extents )
    real_extents=Util.asMultiPolygon(real_extents)
    Logger.writeWkt("Intersected extents",0,real_extents.to_wkt())

    # Now want to find maximum shift outside extents of test.. 
    # Prepare a multipolygon for testing containment.. add a buffer to
    # handle points on edge of region where patch runs against base polygon

    dsmax = trialgrid.maxShiftOutsideAreas( real_extents )
    buffersize = dsmax/Config.base_ramp_tolerance
    Logger.write("Maximum shift outside patch extents {0}".format(dsmax))
    Logger.write("Buffering patches by {0}".format(buffersize))

    buffered_extents=Util.bufferedPolygon( real_extents, buffersize )
    buffered_extents=buffered_extents.intersection( bounds_extents )
    buffered_extents=Util.asMultiPolygon(buffered_extents)
    Logger.writeWkt("Buffered extents",0,buffered_extents.to_wkt())

    # Test areas against land definitions, buffer, and merge overlapping grid definitions

    if Config.land_areas:
        Logger.write("Intersecting areas with land extents".format(len(real_extents)))
        real_extents=real_extents.intersection(Config.land_areas)
        real_extents=Util.asMultiPolygon(real_extents)
        buffered_extents=buffered_extents.intersection( Config.land_areas )
        buffered_extents=Util.asMultiPolygon(buffered_extents)
        Logger.writeWkt("Bounds after land area intersection",0,real_extents.to_wkt())
        Logger.writeWkt("Buffered_bounds after land area intersection",0,buffered_extents.to_wkt())

        # If limiting to land areas then also need to calculate sea areas where
        # lower tolerance applies

        sea_extents = trialgrid.regionsExceedingLevel(
            Config.base_limit_test_column,
            Config.base_limit_sea_tolerance)
        if len(sea_extents) == 0:
            Logger.write("No potential sea extents found")
        else:
            Logger.write("{0} patch sea extents found".format(len(sea_extents)))
            Logger.writeWkt("Sea extents requiring patch (base tolerance)",0,sea_extents.to_wkt())

            # Sea bounds beyond which patch not required because guaranteed by local accuracy bounds
            sea_bounds_extents = trialgrid.regionsExceedingLevel(
                Config.base_limit_bounds_column,Config.base_limit_sea_bounds_tolerance)

            if len(sea_bounds_extents) > 0:
                Logger.writeWkt("Sea absolute bounds on patch",0,sea_bounds_extents.to_wkt())
                sea_extents=sea_extents.intersection( sea_bounds_extents )
                sea_extents=Util.asMultiPolygon(sea_extents)
                Logger.writeWkt("Intersected sea extents",0,sea_extents.to_wkt())

            # Now want to find maximum shift outside extents of test.. 
            # Prepare a multipolygon for testing containment.. add a buffer to
            # handle points on edge of region where patch runs against base polygon

            sea_outside_patch=sea_extents.union( Config.land_areas )
            sea_dsmax = trialgrid.maxShiftOutsideAreas( sea_outside_patch )
            sea_buffersize = sea_dsmax/Config.base_ramp_sea_tolerance
            Logger.write("Maximum shift outside sea patch extents {0}".format(sea_dsmax))
            Logger.write("Buffering sea patches by {0}".format(sea_buffersize))

            sea_buffered_extents=Util.bufferedPolygon( sea_extents, sea_buffersize )
            if len(sea_bounds_extents) > 0:
                sea_buffered_extents=sea_buffered_extents.intersection( sea_bounds_extents )
            sea_buffered_extents=Util.asMultiPolygon(sea_buffered_extents)
            Logger.writeWkt("Sea buffered extents",0,sea_buffered_extents.to_wkt())

            real_extents=Util.asMultiPolygon(real_extents.union(sea_extents))
            buffered_extents=Util.asMultiPolygon(buffered_extents.union(sea_buffered_extents))
            Logger.writeWkt("Bounds after union with sea extents",0,real_extents.to_wkt())
            Logger.writeWkt("Buffered bounds after union with sea extents",0,buffered_extents.to_wkt())

    # Form merged buffered areas

    griddefs=PatchGridList.gridsForAreas(real_extents,
                                        buffered_extents,
                                        cellsize=Config.base_size,
                                        keepWithin=basespec,
                                        level=1,
                                        source=gridsource,
                                       )

    # Now process the extents...

    gridlist=[]
    for griddef in griddefs:
        patchgrids=griddef.splitToSubgrids()
        for p in patchgrids:
            p.base=griddef.name if splitbase else Config.patchname
            p.basename=Config.patchname
            gridlist.append(p)

    # If splitting into a hybrid model then separate out forward and reverse patches.  
    # This also merges grids into parents.
    #
    # Otherwise just need to directly merge grids into parent grids

    
    if Config.model.hybrid:
        Logger.dumpGridList("Grids before hybrid/merge",gridlist)
        hybrid_tol=Config.model.hybridTolerance or Config.default_forward_patch_max_distortion
        gridlist=split_forward_reverse_patches( trialgrid, hybrid_tol, gridlist )
        Logger.dumpGridList("Grids after hybrid/merge",gridlist)
    else:
        for g in gridlist:
            g.mergeIntoParent()
    return gridlist


def build_published_component( gridlist, comppath, cleandir=False, ndp=4 ):
    '''
    Creates the component.csv file and grid components used as a published
    component of a LINZ published deformation model
    '''

    patchdir='patch_'+Config.model.modelname
    modeldesc=Config.model.description
    modeldate=Config.model.modeldate

    # If the model date doesn't contain a date, then append it
    if not re.search(r'[12]\d\d\d[01]\d[0123]\d',patchdir):
        patchdir=patchdir+'_'+modeldate.strftime('%Y%m%d')

    if not gridlist:
        raise RuntimeError("Cannot build shift model - haven't built grid CSV files")
    comppath = os.path.join( comppath, 'model', patchdir )

    Logger.write("Writing published model submodel {0}".format(comppath))
    if not os.path.isdir(comppath):
        os.makedirs(comppath)
    if cleandir:
        for f in os.listdir(comppath):
            fn=os.path.join(comppath,f)
            if not os.path.isdir(fn):
                Logger.write("   Removing existing file {0}".format(f))
                os.remove(fn)


    compcsv=os.path.join(comppath,'component.csv')

    with open(compcsv,"w") as ccsvf:
        ccsv=csv.writer(ccsvf)
        ccsv.writerow(Config.published_component_columns)
        csvdata=dict(
            version_added=0,
            version_revoked=0,
            reverse_patch='N',
            component=0,
            priority=0,
            min_lon=0,
            max_lon=0,
            min_lat=0,
            max_lat=0,
            spatial_complete='Y',
            min_date=0,
            max_date=0,
            time_complete='Y',
            npoints1=0,
            npoints2=0,
            displacement_type='3d',
            error_type='none',
            max_displacement=0,
            spatial_model='llgrid',
            time_function='step',
            time0=modeldate,
            factor0=0,
            time1=modeldate,
            factor1=1,
            decay=0,
            file1='',
            description=modeldesc
        )
        
        ir0=1
        toplevelgrids=[g for g in gridlist if g.isTopLevelGrid()]
        patchversions=Config.model.patchversions
        for nv in range(len(patchversions)):
            versionspec=patchversions[nv]
            reverse_patch=versionspec.reverse
            csvdata['version_added']=versionspec.version
            csvdata['version_revoked']=0
            if nv < len(patchversions)-1:
                csvdata['version_revoked']=patchversions[nv+1].version
            for tg in toplevelgrids:
                ir1 = ir0
                ir0 += len(versionspec.time_model)
                for priority, grid in enumerate(gridlist):
                    if grid.topLevelGrid() != tg:
                        continue
                    rvs=reverse_patch if grid.isforward is None else not grid.isforward
                    csvdata['reverse_patch']='Y' if rvs else 'N'
                    csvname=os.path.split(grid.builtFilename)[1]
                    csvname=os.path.splitext(csvname)[0]+'.csv'
                    csvname = re.sub(r'^(grid_)?','grid_',csvname)
                    compcsv=os.path.join(comppath,csvname)
                    grid.writeCsvFile( compcsv, ndp )
                    
                    gd=defgrid.DeformationGrid(compcsv)
                    min_lon=round(np.min(gd.column('lon')),10)
                    max_lon=round(np.max(gd.column('lon')),10)
                    min_lat=round(np.min(gd.column('lat')),10)
                    max_lat=round(np.max(gd.column('lat')),10)
                    compdata=dict(
                        min_lon=min_lon,
                        max_lon=max_lon,
                        min_lat=min_lat,
                        max_lat=max_lat,
                        npoints1=gd.array.shape[1],
                        npoints2=gd.array.shape[0],
                        max_displacement=round(math.sqrt(np.max(
                            gd.column('de')*gd.column('de') +
                            gd.column('dn')*gd.column('dn') +
                            gd.column('du')*gd.column('du'))),5),
                        file1=csvname,
                        )
                    Logger.writeWkt('Published '+csvname,grid.level,
                            GridSpec(min_lon,min_lat,max_lon,max_lat).boundingPolygonWkt())
                    
                    rdate0=0
                    for ir,r in enumerate(versionspec.time_model):
                        rdate=r.date.strftime("%Y-%m-%d")
                        rvalue=r.factor
                        if ir == 0:
                            if rvalue==0:
                                continue
                            compdata['time_function']='step'
                            compdata['time0']=rdate
                            compdata['time1']=rdate
                        else:
                            compdata['time_function']='ramp'
                            compdata['time0']=rdate0
                            compdata['time1']=rdate
                            rvalue -= versionspec.time_model[ir-1].factor
                        rdate0=rdate
                        compdata['factor0']=-rvalue if rvs else 0
                        compdata['factor1']=0 if rvs else rvalue
                        if not additive:
                            compdata['component']=ir+ir1
                            compdata['priority']=priority
                        csvdata.update(compdata)
                        ccsv.writerow([csvdata[c] for c in Config.published_component_columns])


if __name__ == "__main__":

    # Process arguments

    parser=argparse.ArgumentParser(description='Build set of grids for deformation patch')
    parser.add_argument('patch_file',help='Base name used for output files')
    parser.add_argument('model_file',help='Model file(s) used to calculate deformation, passed to calc_okada')
    # parser.add_argument('--shift-model-path',help="Create a linzshiftmodel in the specified directory")
    parser.add_argument('--submodel-path',help="Create publishable component in the specified directory")
    parser.add_argument('--clean-dir',action='store_true',help="Clean publishable component subdirectory")
    subnest_group=parser.add_mutually_exclusive_group()
    subnest_group.add_argument('--subgrids-nest',action='store_true',help="Grid CSV files subgrids replace parent values to calculate deformation")
    subnest_group.add_argument('--subgrids-additive',action='store_true',help="Grid CSV files subgrids added to calculate total deformation")
#    parser.add_argument('--apply-time-model-scale',action='store_true',help="Scale by the time model final value")
    parser.add_argument('--max-level',type=int,help="Maximum number of split levels to generate (each level increases resolution by 4)")
    parser.add_argument('--base-tolerance',type=float,help="Base level tolerance - depends on base column")
    parser.add_argument('--split-base',action='store_true',help="Base deformation will be split into separate patches if possible")
    parser.add_argument('--no-trim-subgrids',action='store_false',help="Subgrids will not have redundant rows/columns trimmed")
    parser.add_argument('--precision',type=int,default=5,help="Precision (ndp) of output grid displacements")
    parser.add_argument('--land-area',help="WKT file containing area land area over which model must be defined")
    # parser.add_argument('--parcel-shift',action='store_true',help="Configure for calculating parcel_shift rather than rigorous deformation patch")
    # parser.add_argument('--test-settings',action='store_true',help="Configure for testing - generate lower accuracy grids")
    parser.add_argument('--write-grid-lines',action='store_true',help="xx.grid.wkt contains grid cell lines rather than polygon")
    parser.add_argument('-v','--verbose',action='store_true',help="More verbose output (if not already configured)")

    args=parser.parse_args()

    Config.setPatchFilenameRoot(args.patch_file)
        
    split_base=args.split_base
#    shift_path=args.shift_model_path
    comp_path=args.submodel_path
#    if comp_path and shift_path:
#        raise RuntimeError("Cannot create published patch and shift model")
    additive=args.subgrids_additive
    if additive:
        raise RuntimeError("Non-nested grids not supported in current implementation")
    trimgrid=args.no_trim_subgrids
#    if args.parcel_shift: 
#        Logger.write("Configuring model to use parcel shift parameters")
#        args.configureForParcelShift()
#    if args.test_settings: 
#        Logger.write("Configuring model with low accuracy testing settings")
#        args.configureForTesting()
#    if comp_path: 
#        args.apply_time_model_scale=False

    if args.max_level:
        Config.max_split_level=args.max_level

    if args.base_tolerance:
        Config.base_limit_tolerance=args.base_tolerance 
        
    wktgridlines=args.write_grid_lines

    if not os.path.isdir(Config.patchpath):
        os.makedirs(Config.patchpath)

    Logger.createLog(
        Config.filename(extension=".build.log"),
         wktfile=Config.filename(extension=".build.wkt"),
         gridwkt=Config.filename(extension=".grids.wkt"),
         showgridlines=args.write_grid_lines
        )
    if Config.verbose or args.verbose:
        Logger.setVerbose()
    Config.writeTo(Logger.write)

    try:
        if args.land_area:
            Config.loadLandAreas(args.land_area)
                    
        model=FaultModel( args.model_file )
        Config.setFaultModel(model)

        gridlist = build_deformation_grids( splitbase=split_base )

        if comp_path:
            build_published_component( gridlist, comp_path, args.clean_dir, args.precision )

    except Exception as ex:
        Logger.write('\n\nFailed with error: {0}'.format(ex.message),level=-1)
        raise

