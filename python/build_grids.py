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
import sys

calc_okada=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'okada','calc_okada')
if not os.path.exists(calc_okada):
    raise RuntimeError('Cannot find calc_okada program at '+calc_okada)

gridtool=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'gridtool','gridtool')
if not os.path.exists(gridtool):
    raise RuntimeError('Cannot find gridtool program at '+gridtool)

# Log file ...

runtime=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
runtime_version=datetime.datetime.now().strftime("%Y%m%d")
logfile=None
wktfile=None
wktgridfile=None
wktgridlines=False

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

    # Conversion degrees to metres
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
        writelog("    Approx degrees to metres factor: {0} {1}".format(*Config.scale_factor))

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
            print message
            sys.stdout.flush()

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


class Util( object ):

    @staticmethod
    def buffered_polygon( polygon, buffer ):
        '''
        Buffer a polygon by an extents in metres.  Approximately scales to metres, applies buffer,
        and scales back. Returns a multipolygon
        '''
        if buffer > 0:
            polygon=affinity.scale(polygon,xfact=1.0/scale_factor[0],yfact=1.0/scale_factor[1],origin=origin)
            polygon=polygon.buffer(buffer)
            polygon=affinity.scale(polygon,xfact=scale_factor[0],yfact=scale_factor[1],origin=origin)
        if 'geoms' not in dir(polygon):
            polygon=MultiPolygon([polygon])
        return polygon

class GridSpec( object ):
    '''
    Class representing a basic grid definition - min and max coordinates and number of 
    grid cells in each direction.  Note: the number of grid values in each direction is
    one greater than the number of cells.

    Also acts as simple bounding box, with just one grid cell each way
    '''

    @staticmethod
    def fromBounds( bounds ):
        return GridSpec(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1] )

    @staticmethod
    def _calc_min_max( base, minv, maxv, csize, multiple ):
        csize *= multiple
        n0 = int(math.floor((minv-base)/csize))
        n1 = int(math.ceil((maxv-base)/csize))
        return base+n0*csize,base+n1*csize,(n1-n0)*multiple

    def __init__( self, xmin, ymin, xmax, ymax, ngx=1, ngy=1, cellsize=None ):
        if cellsize is not None:
            if ngx == 1:
                ngx=max(int((xmax-xmin)/cellsize[0]+0.5)+1,1)
            if ngy == 1:
                ngy=max(int((ymax-ymin)/cellsize[1]+0.5)+1,1)
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
        return GridSpec(pt0[0],pt0[1],pt1[0],pt1[1],cellsize=self.cellsize())

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

    def lines( self ):
        '''
        Iterator returning lines (as coordinate pairs) forming the internal and
        external edges of the grid.
        '''
        lonv = list(self.xmin+(self.xmax-self.xmin)*float(x)/self.ngx
                for x in range(self.ngx+1))
        latv = list(self.ymin+(self.ymax-self.ymin)*float(x)/self.ngy
                for x in range(self.ngy+1))
        for lat in latv:
            for ln0,ln1 in zip(lonv[:-1],lonv[1:]):
                yield LineString((ln0,lat),(ln1,lat))
        for lon in lonv:
            for lt0,lt1 in zip(latv[:-1],latv[1:]):
                yield LineString((lon,latv[i]),(lon,latv[i+1]))

    def area( self ):
        return (self.xmax-self.xmin)*(self.ymax-self.ymin)

    def mergeIntoList( self, deflist ):
        while True:
            overlap = None
            for def1 in deflist:
                if self.overlaps(def1):
                    overlap=def1
                    break
            if overlap:
                deflist.remove(overlap)
                self=self.mergeWith(overlap)
            else:
                break
        deflist.append(self)

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
        
        self.modelname=os.path.basename(self.modelpath)
        with open(self.modelpath) as f:
            event = f.readline()
            model = f.readline()
            version = f.readline()
            description = event+model+version
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
            call(params)
            if not os.path.exists(gridfile):
                raise RuntimeError("Failed to calculate grid file {0}\n{1}".
                                   format(gridfile," ".join(params)))
            with open(metafile,'w') as f:
                f.write(meta)
        else:
            print "          Reloading existing grid {0} ...".format(gridfile)
        grid=defgrid.DeformationGrid(gridfile)
        return grid

class PatchGridDef:
                                 
    def __init__( self, name=None, subpatch=None, extents=None, buffer=0.0, level=1, spec=None, parent=None ):
        self.buffer=0.0
        self.level=level
        self.parent=parent
        self.name=name
        self.subpatch=subpatch
        self._extents=extents
        self._extentsWkt=None
        self._bufferedExtents=None
        self._bufferedExtentsWkt=None
        self.spec=spec
        self._grid=None
        self._strainGrid=None
        self.csv=None
        self.gdf=None
        self.base=None
        self.basename=None
        self.isforward=None  # Component is forward/reverse in hybrid grid else None (= either)
        self._buffered_extents=None 

    def update( self, level=None, file=None, parent=None, extentfile=None, isforward=None ):
        pgd=PatchGridDef(
            level=level if level is not None else self.level,
            file=file if file is not None else self.file,
            extentfile=extentfile if extentfile is not None else self.extentfile,
            parent=parent)
        pgd.isforward=isforward
        pgd.base=self.base
        pgd.basename=self.basename
        return pgd

    # This is not quite correct - model should be a procedure or object that
    # is called with the grid spec ... may not alway be using Okada.

    def calcGrid( self, modeldef=None ):
        gridfile=self._fileWithExtension('.grid')
        return Config.model.calcGrid( self.spec, gridfile )

    @property
    def grid( self ):
        if self.spec is None:
            raise RuntimeError('{0} grid spec not calculated before calculating grid values'
                               .format(self.name))
        if self._grid is None:
            self._grid=self.calcGrid()
        return self._grid

    def _fileWithExtension( self, extension ):
        return os.path.join(Config.patchpath,self.name+extension)

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
    def bufferedExtents( self, buffer=None ):
        '''
        Return the buffered extents over which the grid may differ from its
        parent.

        Returned as a single MultiPolygon
        '''
        if self.buffer != buffer:
            self._bufferedExtents=None
            self._bufferedExtentsWkt=None
            self.buffer=buffer
        if self._bufferedExtents is None and self.extents is not none:
            self._bufferedExtents=Util.bufferedPolygon(self.extents,self.buffer)
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

    def fitToParent( self ):
        '''
        Modify the grid to match transition to the parent grid between
        its extents and its buffered extents
        '''
        if self.parent is None:
            return
        raise NotImplementedError('... still needed!')

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
        grid=self.grid()
        gridr = grid.calcResolution(tolerance,precision=0.0000001)
        subcell_areas=gridr.regionsLessThanLevel('reqsize',celldimension )
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
            grid2=PatchGridDef(
                name="{0}_SL{1}".format(self.name,self.level+1),
                level=self.level+1,
                spec=self.spec.split(Config.cell_split_factor))

            # Determine required resolution of subcells, and find the extents for which 
            # the subcell resolution is required (ie the required resolution is less than
            # the parent cell size)

            resolution=self.spec.cellDimension()
            subcell_areas=grid2.areasNeedingFinerGrid( 
                Config.subcell_resolution_tolerance, resolution)

            if land_areas:
                valid_areas=[]
                for sc in subcell_areas:
                    p=sc.intersection(land_areas)
                    if 'geoms' not in dir(p):
                        p=[p]
                    for g in p:
                        if type(g) == Polygon:
                            valid_areas.append(g)
                subcell_areas=valid_areas

                # If limiting to land areas, then also select areas where sea resolution
                # is to be applied.
            
                subcell_areas.extend(grid2.areasNeedingFinerGrid( 
                    Config.subcell_resolution_sea_tolerance,
                    resolution))

            # Buffer to expand subcell areas by parent grid size (resolution) 
            # This is to provide a sufficient buffer for smoothing on to
            # parent grid without 

            # Note: previous version of code provided for replacing parent 
            # grid 

            subcells=PatchGridDef.gridsForAreas( 
                subcell_areas, 
                subcellsize,
                buffer=resolution,
                keepWithin=self.spec,
                level=self.level+1,
                subpatch=self.subpatch
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

    def __str__( self ):
        return "\n".join((
            "Grid Name: {0}".format(self.name),
            "    Parent grid: {0}".format(self.parent.name if self.parent else ""),
            "    Grid file: {0}".format(self.file),
            "    Extent file: {0}".format(self.extentfile),
            "    Level: {0}".format(self.level),
            "    Csv: {0}".format(self.csv),
            "    Gdf: {0}".format(self.gdf),
            "    Base: {0}".format(self.base),
            "    Base name: {0}".format(self.basename),
            "    Is Forward?: {0}".format(self.isforward),
            ))

class PatchGridList( list ):

    @staticmethod
    def gridsForAreas( areas, cellsize, buffer=0.0, keepWithin=None, level=1, 
                      name=None,subpatch=''):
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
        merged_specs=[]
        for a in areas:
            bounds = GridSpec.fromBounds(a.bounds).expand(cellsize)
            if keepWithin is not None:
                bounds=bounds.intersectWith(keepWithin)
            bounds=bounds.alignTo(cellsize)
            bounds.mergeIntoList( merged_specs )
        merged_specs.sort(key=lambda x: x.griddef[0][1] )

        if name is None:
            if len(merged_specs > 1):
                subpatch=(subpatch or '')+"P{0}"
            name=Config.patchname
            if subpatch:
                name=name+"_"+subpatch
            name=name+"_L{0}".format(level)

        boundsbuffer=(buffer/Config.scale_factor[0],buffer/Config.scale_factor[1])

        patchdefs=PatchGridList()
        for i,b in enumerate(merged_specs):
            extents=MultiPolygon([a for a in areas if a.intersects(b.boundingPolygon())])
            patchname=name.format(i)
            spec=GridSpec.fromBounds(bounds).expand(boundsbuffer).alignedTo(cellsize)
            patchdefs.append(PatchGridDef(name=patchname,subpatch=subpatch,spec=spec,extents=extents,buffer=buffer,level=level))

        return patchdefs


def grid_file_version( gridfile, version ):
    path,ext=os.path.splitext(gridfile)
    return path+version+ext

def split_forward_reverse_patches( patchpath, patchname, trialgrid, hybrid_tol, gridlist ):
    Logger.write("Splitting patch grids into forward/reverse patches")
    Logger.write("Forward patch tolerated distortion (ppm): {0}".format(hybrid_tol))
    reverse_extents = regionsExceedingLevel(trialgrid, forward_patch_test_column,hybrid_tol)
    Logger.write("{0} regions found".format(len(reverse_extents)))
    Logger.writeWkt("Extents requiring reverse patch (base tolerance)",0,reverse_extents.to_wkt())
    reverse_grid_extents=[g.bufferedExtents() 
                          for g in gridlist if g.level >= forward_patch_max_level]
    if reverse_grid_extents:
        reverse_grid_extents.append(reverse_extents)
        reverse_extents=unary_union(reverse_grid_extents)
        if type(reverse_extents) == Polygon:
            reverse_extents=MultiPolygon([reverse_extents])

    reversewkt=os.path.join(patchpath,patchname+'_reverse_patch_extents.wkt')
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

        if (g.parent in not_reverse 
                or not g.bufferedExtents().intersects(reverse_extents)):
            if g.level > forward_patch_max_level:
                raise RuntimeError('Forward patch required for level greater than {0}'
                                   .format(forward_patch_max_level))
            not_reverse.add(g)
            Logger.write('Copying {0} as forward patch grid only'.format(g.name))
            fgridfile=grid_file_version(g.name,'_F')
            gf=g.update(isforward=True,file=fgridfile,parent=forward_grids.get(gc.parent,None))
            gf.mergeIntoParent()
            forward_grids[g]=gf
            splitgrids.append(gcf)

        # Create smoothed grid for forward patch, smoothing over reverse patch
        # extents
        #
        # Have tried different smoothing options (including datumgrid), but
        # none seems to offer particular advantage.  Could be an area for more
        # research

        if g.level <= forward_patch_max_level:
            fgridfile=grid_file_version(g.file,'_F')
            basefg=forward_grids.get(g.parent,None)

            Logger.write('Creating smoothed forward patch grid for {0}'.format(g.name))
            commands=[
                gridtool,
                'read','maxcols','3',g.file
                ]
            if basefg is not None:
                commands.extend(['replace','maxcols','3',basefg.file,'where','inside',reversewkt])
            commands.extend(['smooth','linear','inside',reversewkt])
            if basefg is not None:
                commands.extend(['not','on_grid',basefg.file])
            commands.extend(['write',fgridfile])
        
            call(commands)

            if not os.path.exists(fgridfile):
                raise RuntimeError('Cannot build smoothed forward patch '+fgridfile)
            Logger.write('Created forward patch for {0}'.format(fgridfile))

            gf=g.update(name=fgridfile,isforward=True)
            gf.mergeIntoParent()
            forward_grids[g]=gf
            splitgrids.append(gf)
        else:
            # If not creating a forward patch, then use the forward patch
            # associated with parent grid
            forward_grids[g]=forward_grids.get(g.parent,None)


        # Now create the reverse patch files by subtracting the forward patch
        # from them

        rgridfile=grid_file_version(g.file,'_R')
        fgrid=forward_grids[g]
        Logger.write('Creating reverse patch grid for {0}'.format(g.name))
        commands=[
            gridtool,
            'read','maxcols','3',g.file,
            'subtract',fgrid.file,
            'trim','1',
            'write',rgridfile]
        call(commands)
        if not os.path.exists(rgridfile):
            raise RuntimeError('Cannot build reverse patch '+rgridfile)
        rbase=reverse_grids.get(g,None)
        rgrid=g.update(file=rgridfile,parent=reverse_grids.get(g.parent,None),isforward=False)
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

    basespec=GridSpec.fromBounds(Config.base_extents)
    trialgriddef=basespec.alignedTo(Config.base_size)
    trialgridfile=Config.filename(extension="_trial.grid")

    Logger.write("Building trial grid using {0} to {1}".format(
        trialgriddef,trialgridfile))
    Logger.writeWkt("Trial extents",0,trialgriddef.boundingPolygonWkt())

    trialgrid=PatchGridDef(trialgridfile,spec=trialgriddef)

    # Determine the extents on which the base tolerance is exceed - returns
    # a list of polygons

    real_extents = trialgrid.regionsExceedingLevel(Config.base_limit_test_column,Config.base_limit_tolerance)
    Logger.write("{0} patch extents found".format(len(real_extents)))
    Logger.writeWkt("Extents requiring patch (base tolerance)",0,real_extents.to_wkt())

    # Bounds beyond which patch not required because guaranteed by local accuracy bounds
    bounds_extents = trialgrid.regionsExceedingLevel(base_limit_bounds_column,base_limit_bounds_tolerance)

    # Now want to find maximum shift outside extents of test.. 
    # Prepare a multipolygon for testing containment.. add a buffer to
    # handle points on edge of region where patch runs against base polygon

    dsmax = trialgrid.maxShiftOutsideAreas( real_extents )
    buffersize = dsmax/base_ramp_tolerance
    Logger.write("Maximum shift outside patch extents {0}".format(dsmax))
    Logger.write("Buffering patches by {0}".format(buffersize))

    buffered_extent=Util.bufferedPolygon( real_extents, buffersize )
    Logger.writeWkt("Buffered extents",0,buffered_extent.to_wkt())
    Logger.writeWkt("Absolute bounds on patch",0,bounds_extents.to_wkt())

    real_extents=buffered_extent.intersection( bounds_extents )
    if 'geoms' not in dir(real_extents):
        real_extents = MultiPolygon([real_extents])
    Logger.writeWkt("Bounds before potential land intersection",0,real_extents.to_wkt())


    # Test areas against land definitions, buffer, and merge overlapping grid definitions

    if Config.land_areas:
        Logger.write("Intersecting areas with land extents".format(len(real_extents)))
        valid_areas=[]
        for sc in real_extents:
            p=sc.intersection(Config.land_areas)
            if 'geoms' not in dir(p):
                p=[p]
            for g in p:
                if type(g) == Polygon:
                    valid_areas.append(g)
        real_extents=MultiPolygon(valid_areas)
        Logger.writeWkt("Bounds after land area intersection",0,real_extents.to_wkt())

        # If limiting to land areas then also need to calculate sea areas where
        # lower tolerance applies

        sea_extents = trialgrid.regionsExceedingLevel(base_limit_test_column,base_limit_sea_tolerance)
        if len(sea_extents) == 0:
            Logger.write("No potential sea extents found")
        else:
            Logger.write("{0} patch sea extents found".format(len(sea_extents)))
            Logger.writeWkt("Sea extents requiring patch (base tolerance)",0,sea_extents.to_wkt())

            # Sea bounds beyond which patch not required because guaranteed by local accuracy bounds
            sea_bounds_extents = trialgrid.regionsExceedingLevel(
                base_limit_bounds_column,base_limit_sea_bounds_tolerance)

            # Now want to find maximum shift outside extents of test.. 
            # Prepare a multipolygon for testing containment.. add a buffer to
            # handle points on edge of region where patch runs against base polygon

            sea_outside_patch=sea_extents.union( land_area )
            sea_dsmax = trialgrid.maxShiftOutsideAreas( sea_outside_patch )
            sea_buffersize = sea_dsmax/base_ramp_sea_tolerance
            Logger.write("Maximum shift outside sea patch extents {0}".format(sea_dsmax))
            Logger.write("Buffering sea patches by {0}".format(sea_buffersize))

            sea_buffered_extent=Util.bufferedPolygon( sea_extents, sea_buffersize )
            Logger.writeWkt("Sea buffered extents",0,sea_buffered_extent.to_wkt())

            if len(sea_bounds_extents) > 0:
                Logger.writeWkt("Sea absolute bounds on patch",0,sea_bounds_extents.to_wkt())
                sea_extents=sea_buffered_extent.intersection( sea_bounds_extents) 
                if 'geoms' not in dir(sea_extents):
                    sea_extents = MultiPolygon([sea_extents])

            Logger.writeWkt("Sea bounds before union with land extents",0,sea_extents.to_wkt())

            real_extents=real_extents.union(sea_extents)
            if 'geoms' not in dir(real_extents):
                real_extents=MultiPolygon([real_extents])
            Logger.writeWkt("Bounds after union with sea extents",0,real_extents.to_wkt())

    # Form merged buffered areas

    griddefs=PatchGridDef.gridsForAreas(real_extents,
                                        cellsize=Config.base_size,
                                        keepWithin=basespec,
                                        level=1
                                       )

    # Now process the extents...

    gridlist=[]
    for griddef in enumerate(griddefs):
        patchgrids=griddef.splitToSubgrids()
        for p in patchgrids:
            p.base=griddef.name if splitbase else Config.patchname
            p.basename=Config.patchname
            gridlist.append(p)
    if Config.model.hybrid:
        hybrid_tol=Config.model.hybridTolerance
        gridlist=split_forward_reverse_patches( patchpath, patchname, trialgrid, hybrid_tol, gridlist )
    return gridlist

def create_patch_csv( patchlist, modelname, additive=True, trimgrid=True, linzgrid=False, precision=5 ):
    '''
    Function takes a set of grids of deformation data and compiles CSV and optionally LINZ 
    ascii grid format files defining the patch.  Each grid file has associated extent and 
    buffer WKT files, which define the extents over which they are actually required (ie
    within which the parent grid resolution is not adequate).  The buffer wkt defines the
    extents over which the grid deformation is smoothed into that of the parent grid.
    So the main work of this routine is to smooth each grid into its parent grid.

    By default assumes that the patch is "additive", that is each component is
    added together to build the total deformation.  If additive is false, then the 
    CSV files are constructed such that each subgrid overrides its parent grid over the
    area in which it applies.

    By default trims grids to remove redundant zero rows/columns at edges.  This 
    can be turned off to keep grid boundaries aligned with parent grid (makes no
    difference to calculated values - just looks nicer!?)
    '''
    # Note: Assumes patchlist is sorted such that parents are before children
    Logger.write("Creating patch CSV files")
    for patch in patchlist:
        # Get the files used for this component of the patch

        gridfile = patch.file
        component = patch.name
        csvfile = re.sub(r'(\.grid)?$','.csv',gridfile,re.I)
        gdffile = re.sub(r'(\.grid)?$','.gdf',gridfile,re.I)
        Logger.write("Creating patch file file {0}".format(csvfile))
        parent = patch.parent
        extfile = patch.extentfile+'.extent.wkt'
        buffile = patch.extentfile+'.buffer.wkt'
        gridparams = dict(
            gridfile=gridfile,
            csvfile=csvfile,
            gdffile=gdffile,
            parentgrid=parent.file if parent else 'undefined',
            parentcsv=parent.csv if parent else 'undefined',
            extents=extfile,
            buffer=buffile,
            header1='Patch component '+component.replace('"',"'"),
            header2='Model source '+modelname.replace('"',"'"),
            header3=runtime
        )

        # Build a grid tool command to process the file
        
        # Read in the existing file (only need de, dn, du components)
        merge_commands=gridtool+'\nread maxcols 3 gridfile '
 
        # If have a parent, then subtract it as want to merge into it
        if parent: merge_commands = merge_commands + '\nsubtract maxcols 3 parentgrid'

        # Set to zero outside area extents for which patch subcell resolution is required
        merge_commands=merge_commands + '\nzero outside extents'

        # Ensure zero at edge
        merge_commands=merge_commands + '\nzero edge 1'

        # Smooth out to buffer to avoid sharp transition
        merge_commands=merge_commands + '\nsmooth linear outside extents not outside buffer not edge 1'

        # Trim redundant zeros around the edge of the patch
        merge_commands=merge_commands + '\ntrim 1'

        # If LINZ grid required then write it out.  LINZ grid is always additive
        # so this comes before adding the parent grid back on
        if linzgrid:
            merge_commands=merge_commands + '\nwrite_linzgrid NZGD2000 header1 header2 header3 gdffile'

        # If we have a parent and are not making additive patches then add the parent back
        # on.
        if parent and not additive: merge_commands = merge_commands + '\nadd csv parentcsv'

        # Set the output coordinate precision
        merge_commands = merge_commands + '\nprecision {0}'.format(precision)

        # Finally write out as a CSV file
        merge_commands=merge_commands + '\nwrite csv csvfile'

        # Run the commands
        gridtoolcmd=[gridparams.get(x,x) for x in merge_commands.split() ]

        def expand_param(m):
            x=m.group(1)
            x=gridparams.get(x,x)
            x='"'+x+'"' if re.search(r'\s',x) else x
            return x
    
        Logger.write("Running "+
                  re.sub(r'(\w+)',expand_param,merge_commands)
                  .replace("\n","\n    ")
                   )

        call(gridtoolcmd)

        if not os.path.exists( csvfile ):
            raise RuntimeError("Could not create CSV file "+csvfile)
        patch.csv=csvfile
        if linzgrid:
            if not os.path.exists( gdffile ):
                raise RuntimeError("Could not create CSV file "+gdffile)
            patch.gdf=gdffile


def build_linzshift_model( gridlist, modeldef, additive, shiftpath, splitbase=False ):
    '''
    Creates shift definition files using list of grids.  Assumes that grids have
    been generated using linzgrid=True
    '''

    if not gridlist or not gridlist[0].gdf:
        raise RuntimeError("Cannot build shift model - haven't built LINZ grid format files")
    if not additive:
        raise RuntimeError("Cannot build shift model - grids must be built as additive")

    shiftname=os.path.basename(shiftpath)
    shiftpath = os.path.dirname(shiftpath)
    if not os.path.isdir(shiftpath):
        os.makedirs(shiftpath)
    patchname=gridlist[0].basename
    plen=len(patchname)
    baselist=set([patch.base for patch in gridlist])
    for base in baselist:
        shiftdef=shiftname+base[plen:]+'_shift.def'
        shiftbin=shiftname+base[plen:]+'_shift.bin'
        shiftfile=os.path.join(shiftpath,shiftdef)
        Logger.write("Creating shift definition file {0}".format(shiftfile))
        with open(shiftfile,'w') as f:
            replace=dict(
                name=base,
                runtime=runtime,
                modelname=modelname
            )
            f.write(re.sub(
                r'\{(\w+)\}',
                lambda m: replace[m.group(1)],
                shift_model_header
            ))
            for patch in gridlist:
                if patch.base != base and splitbase:
                    continue
                gdf=os.path.basename(patch.gdf)
                gdf=shiftname+gdf[plen:]
                gdffile=os.path.join(shiftpath,gdf)
                if gdffile != patch.gdf:
                    shutil.copyfile(patch.gdf,gdffile)
                replace['component']=patch.name
                replace['gridfile']=gdf
                f.write(re.sub(
                    r'\{(\w+)\}',
                    lambda m: replace[m.group(1)],
                    shift_model_component
                ))
        call(['makelinzshiftmodel.pl',shiftdef,shiftbin],cwd=shiftpath)


def build_published_component( gridlist, modeldef, modelname, additive, comppath, cleandir=False ):
    '''
    Creates the component.csv file and grid components used as a published
    component of a LINZ published deformation model
    '''

    patchversions=get_patch_spec( modeldef )
    if not patchversions:
        raise RuntimeError("Patch definition file missing for {0}".format(modeldef))
    patchdir='patch_'+modelname
    modeldesc=patchversions[0].event
    modeldate=patchversions[0].date

    # If the model date doesn't contain a date, then append it
    if not re.search(r'[12]\d\d\d[01]\d[0123]\d',patchdir):
        patchdir=patchdir+'_'+modeldate.strftime('%Y%m%d')

    if not gridlist or not gridlist[0].csv:
        raise RuntimeError("Cannot build shift model - haven't built grid CSV files")
    comppath = os.path.join( comppath, 'model', patchdir )
    if not os.path.isdir(comppath):
        os.makedirs(comppath)
    if cleandir:
        for f in os.listdir(comppath):
            fn=os.path.join(comppath,f)
            if not os.path.isdir(fn):
                os.remove(fn)

    Logger.write("Writing published model submodel {0}".format(comppath))

    compcsv=os.path.join(comppath,'component.csv')

    with open(compcsv,"w") as ccsvf:
        ccsv=csv.writer(ccsvf)
        ccsv.writerow(published_component_columns)
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
                    gd=defgrid.DeformationGrid(grid.csv)
                    (gridpath,csvname)=os.path.split(grid.csv)
                    rvs=reverse_patch if grid.isforward is None else not grid.isforward
                    csvdata['reverse_patch']='Y' if rvs else 'N'
                    csvname = re.sub(r'^(grid_)?','grid_',csvname)
                    compcsv=os.path.join(comppath,csvname)
                    if compcsv != grid.csv:
                        shutil.copyfile(grid.csv,compcsv)
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
                              bounds_poly(((min_lon,min_lat),(max_lon,max_lat))).to_wkt())
                    
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
                        ccsv.writerow([csvdata[c] for c in published_component_columns])


if __name__ == "__main__":

    # Process arguments

    parser=argparse.ArgumentParser(description='Build set of grids for deformation patch')
    parser.add_argument('patch_file',help='Base name used for output files')
    parser.add_argument('model_file',help='Model file(s) used to calculate deformation, passed to calc_okada')
    # parser.add_argument('--shift-model-path',help="Create a linzshiftmodel in the specified directory")
    parser.add_argument('--submodel-path',help="Create publishable component in the specified directory")
    parser.add_argument('--clean-dir',action='store_true',help="Clean publishable component subdirectory")
    parser.add_argument('--subgrids-nest',action='store_false',help="Grid CSV files calculated to replace each other rather than total to deformation")
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
    additive=args.subgrids_nest
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

#        build_linzgrid=bool(shift_path)

        create_patch_csv( gridlist, modelname, linzgrid=build_linzgrid, additive=additive, trimgrid=trimgrid, precision=args.precision )

#        if shift_path:
#            build_linzshift_model( gridlist, modeldef, additive, shiftpath=shift_path, splitbase=split_base )

        if comp_path:
            build_published_component( gridlist, modeldef, modelname, additive, comp_path, args.clean_dir )
    except Exception as ex:
        Logger.write('\n\nFailed with error: {0}'.format(ex.message),level=-1)
        raise

