#!/usr/bin/python
# Script to build a set of nested grids to meet tolerances etc.

from collections import namedtuple
from shapely.geometry import MultiPolygon,Polygon,LineString,Point,shape
from shapely.prepared import prep
from shapely.ops import unary_union
from shapely import affinity
from shapely.wkt import loads as wkt_load
from subprocess import call, check_output
import accuracy_standards
import argparse
import copy
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

Unspecified=object()

class PatchGridCriteria( object ):

    '''
    Accuracy criteria governing grid creation

    Base limit are criteria defining the total extent of the patch based
    on the area in which both the distortion and the size of deformation 
    are considered significant.
    Subcell ramp criteria define the ramp used to merge the total
    patch down to zero at the boundary
    
    Grid level split criteria defines the allowable gridding error before 
    it is required to be split into a subgrid 
   
    Forward patch criteria define the partitioning of deformation into
    forward and reverse patch components
    '''

    fields= (
    'base_limit_test_column',           # Column defining distortion
    'base_limit_tolerance',             # Land based/default maximum ignorable distortion
    'base_limit_sea_tolerance',         # Sea based maximum ignorable distortion

    'base_limit_bounds_column',         # Column defining size of deformation
    'base_limit_bounds_tolerance',      # Land based/default maximum ignorable deformation
    'base_limit_sea_bounds_tolerance',  # Se abased maximum ignorable deformation

    'ramp_tolerance',           # Maximum permitted deformation at edge of patch
                                        # over land.
    'ramp_sea_tolerance',       #   "    " over sea
    'grid_level_split_criteria',        # Acceptable gridding error for a grid over land
                                        # before it is split into a subgrid
    'grid_level_split_sea_criteria',    #   "    " over sea
    'max_split_level',                  # Maximum depth to which cells are split

    'forward_patch_test_column',        # Criteria column for partitioning forward/reverse
    'forward_patch_max_distortion',     # Maximum allowable value for test value
    'forward_patch_max_level',          # Maximum allowable level of forward patch
    'subsplit_using_vertical',          # True if splitting is based on vertical resolution
    )

    # Factor by which accuracy standards are multiplied to derive 
    # criteria used in forming grids.  Based rather perversely on the 
    # assumption that some of the tolerance can be assigned to the 
    # modelling/gridding process and the rest is left to the source 
    # uncertainty of coordinates calculated using the patch
    #
    # Basically just a way of trying to include some rationale for 
    # numeric criteria.

    tolerance_factor=0.4 

    def __init__( self, **values ):
        for f in self.fields:
            setattr(self,f,values.pop(f))
        for k in values:
            raise RuntimeError('Invalid field {0}'.format(k))

    def copy( self ):
        return copy.copy(self)

    def update( self, **values ):
        for f in self.fields:
            if f in values:
                setattr(self,f,values.pop(f))
        for k in values:
            raise RuntimeError('Invalid field {0}'.format(k))
        return self


class Config( object ):
    '''
    Class used to contain configuration information for the process.  This 
    also provides a few static methods for generating information related to
    the configuration, such as filenames.
    '''

    patch_filename_templates=('{0}','{0}.patch','fault_models/{0}','fault_models/{0}.patch')
    model_filename_templates=('{0}','{0}.model','fault_models/{0}','fault_models/{0}.model')

    verbose=True
    patchpath=''
    patchroot='patch'
    patchname='patch'
    patchversions=None
        
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

    # Division of grids into subgrids 
    cell_split_factor=4

    # Conversion metres to degrees
    scale_factor=(1.0e-5/math.cos(math.radians(40)),1.0e-5)

    # Areas outside 
    land_areas=None

    # Criteria defining the extents and resolution of grids
    tolerance_factor=PatchGridCriteria.tolerance_factor

    # Grids will be discarded if expanding the child grid to cover 
    # it does not increase the child by more than this amount.

    redundantParentRatio=1.1

    # Splits published grids into approximate size below 
    # with tolerance to allow slightly bigger to be allowed,
    # eg split_size=100 with tolerance 0.2 allows up to 120 
    # before splitting.

    published_grid_split_size=100
    published_grid_split_tol=None

    horizontal_grid_criteria=PatchGridCriteria(
        base_limit_test_column='err',
        base_limit_tolerance=accuracy_standards.NRF.lap*tolerance_factor*1.0e6,
        base_limit_sea_tolerance=accuracy_standards.BGN.lap*tolerance_factor*1.0e6,

        base_limit_bounds_column='ds',
        base_limit_bounds_tolerance=accuracy_standards.NRF.lac*tolerance_factor,
        base_limit_sea_bounds_tolerance=accuracy_standards.BGN.na*tolerance_factor,

        ramp_tolerance=accuracy_standards.NRF.lap*tolerance_factor,
        ramp_sea_tolerance=accuracy_standards.BGN.lap*tolerance_factor,

        grid_level_split_criteria=accuracy_standards.NRF.lac*tolerance_factor,
        grid_level_split_sea_criteria=accuracy_standards.BGN.na*tolerance_factor,
        max_split_level=5,

        forward_patch_test_column='err',
        forward_patch_max_distortion=10,
        forward_patch_max_level=2,
        subsplit_using_vertical=False,
        )

    vertical_grid_criteria=horizontal_grid_criteria.copy().update(
        base_limit_test_column='tiltmax',
        base_limit_bounds_column='du',
        forward_patch_test_column='tiltmax',
        subsplit_using_vertical=True,
        )

    # Replace entire grid with subcells if more than this percentage 
    # of area requires splitting
    #
    # max_subcell_ratio=0.5
    #
    # This grid optimization is not currently enabled..

    # Hybrid patch grid maximum ppm distortion in forward patch


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
        Config.horizontal_grid_criteria.update( 
            base_limit_test_column='ds',
            base_limit_tolerance=0.05   
            )
        Config.vertical_grid_criteria.update( 
            base_limit_test_column='du',
            base_limit_tolerance=0.05   
            )

    # Reduce the accuracies to produce smaller data sets for testing...
    @staticmethod
    def configureForTesting():
        for gc in (Config.horizontal_grid_criteria,Config.vertical_grid_criteria):
            gc.grid_level_split_criteria *= 100
            gc.max_split_level -= 2

    # Set patch file name
    @staticmethod
    def setPatchFilenameRoot( patchroot ):
        Config.patchroot=patchroot
        Config.patchpath=os.path.dirname(patchroot)
        Config.patchname=os.path.basename(patchroot)

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
        writelog("    Grid cell split factor: {0}".format(Config.cell_split_factor))
        writelog("    Approx metres to degrees factor: {0} {1}".format(*Config.scale_factor))
        for name,grid_criteria in (
            ('Horizontal',Config.horizontal_grid_criteria),
            ('Vertical',Config.vertical_grid_criteria)):
            writelog("    {0} grid criteria".format(name))
            writelog("        Tolerated distortion outside patch (ppm): land {0} sea {1}".format( 
                grid_criteria.base_limit_tolerance, grid_criteria.base_limit_sea_tolerance))
            writelog("        Tolerated dislocation outside patch (mm): land {0} sea {1}".format( 
                grid_criteria.base_limit_bounds_tolerance*1000, grid_criteria.base_limit_sea_bounds_tolerance*1000))
            writelog("        Patch ramp tolerance (ppm): land {0} sea {1}".format(
                grid_criteria.ramp_tolerance,grid_criteria.ramp_sea_tolerance))
            writelog("        Tolerated gridding error before splitting (mm)".format(
                grid_criteria.grid_level_split_criteria*1000, grid_criteria.grid_level_split_sea_criteria*1000))
            writelog("        Maximum split level: {0}".format(grid_criteria.max_split_level))
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
            Logger.wktfile.write("Item|level|wkt\n")
        if gridwkt:
            Logger.wktgridfile=open(gridwkt,'w')
            Logger.wktgridfile.write("Name|level|nrow|ncol|wkt\n")
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
        self.xmin=float(xmin)
        self.ymin=float(ymin)
        self.xmax=float(xmax)
        self.ymax=float(ymax)
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

    def nodeCount( self ):
        return (self.ngx+1)*(self.ngy+1)

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

    def xy( self, gx, gy ):
        x=self.xmin+float((self.xmax-self.xmin)*gx)/self.ngx;
        y=self.ymin+float((self.ymax-self.ymin)*gy)/self.ngy;
        return x,y

    def gxy( self, x, y ):
        gx=int(math.ceil((x-self.xmin)*self.ngx/(self.xmax-self.xmin)-0.5))
        gy=int(math.ceil((y-self.ymin)*self.ngy/(self.ymax-self.ymin)-0.5))
        return gx,gy

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

    def contains( self, other ):
        return (
            self.xmin <= other.xmin and
            self.xmax >= other.xmax and
            self.ymin <= other.ymin and
            self.ymax >= other.ymax)

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

    @staticmethod
    def _splitn( nx, ns, tolerance ):
        ng=max(int(float(nx-1)/ns-tolerance),0)+1
        nr=int((nx-1)/ng)+1
        result=[(i,min(i+nr,nx)) for i in range(0,nx,nr)]
        return result

    def splitGridOnSize( self, splitrows=None, tolerance=None ):
        if tolerance is None:
            tolerance=0.5 
        if splitrows is None:
            return [self]
        splitx=self._splitn(self.ngx,splitrows,tolerance)
        splity=self._splitn(self.ngy,splitrows,tolerance)
        if len(splitx) <= 1 and len(splity) <= 1:
            return [self]
        specs=[]
        for nxmin,nxmax in splitx:
            for nymin,nymax in splity:
                xmin,ymin=self.xy(nxmin,nymin)
                xmax,ymax=self.xy(nxmax,nymax)
                specs.append(GridSpec(xmin,ymin,xmax,ymax,nxmax-nxmin,nymax-nymin))
        return specs

class PatchVersion( object ) :

    TimePoint=namedtuple('TimePoint','date factor')

    @staticmethod
    def loadPatchDefinition( patchfile, version=None ):
        '''
        Loads a patch specification, which may include multiple versions 
        of the patch file.
        '''
        for template in Config.patch_filename_templates:
            filename=template.format(patchfile)
            if os.path.exists(filename):
                patchfile=filename
                break
        if not os.path.exists(patchfile):
            raise RuntimeError('Cannot find patch definition file {0}'.format(patchfile))

        patchversions=[]
        version=None
        event=None
        eventdate=None
        patchdir=os.path.dirname(patchfile)
        patchname=os.path.splitext(os.path.basename(patchfile))[0]
        default_model_filename=patchname
        model_filename=None
        loaded_model_filename=None
        model=None
        patch_type=None
        hybrid_tol=None
        nested=None
        time_model=None
        in_time_model=False
        overrides={}

        with open(patchfile) as f:
            lineno=0
            for l in f:
                lineno += 1
                try:
                    # Ignore comments and blank lines
                    if re.match(r'^\s*(\#|$)',l):
                        continue
                    # split out command and value...
                    cmdmatch=re.match(r'^\s*(\w+)\s*\:\s*(.*?)\s*$',l)
                    if cmdmatch:
                        in_time_model=False
                        command=cmdmatch.group(1).lower()
                        value=cmdmatch.group(2)
                        if command == 'event':
                            if version is not None:
                                raise RuntimeError('Event must be defined before versions')
                            elif event is None:
                                event=value
                            else:
                                raise RuntimeError('Patch cannot define multiple events')
                            continue
                        elif command == 'date':
                            if version is not None:
                                raise RuntimeError('Date must be defined before versions')
                            elif eventdate is None:
                                eventdate=PatchVersion._parseDate(value)
                            else:
                                raise RuntimeError('Patch cannot define multiple dates')
                            continue

                        if command == 'model':
                            model_filename=value
                            continue

                        if command == 'version':
                            if not re.match(r'^[12]\d\d\d[01]\d[0123]\d$',value):
                                raise RuntimeError("Invalid deformation model version {0} - must by YYYYDDMM"
                                                   .format(version))
                            if version is not None:
                                if value <= version:
                                    raise RuntimeError("Patch versions out of order: {1} <= {2}"
                                                       .format(value,version))
                                model_filename=model_filename or default_model_filename
                                model_filename=os.path.join(patchdir,model_filename)
                                if model_filename != loaded_model_filename:
                                    loaded_model_filename=model_filename
                                    model=FaultModel(model_filename)
                                patchversions.append(PatchVersion(
                                    patchname=patchname,
                                    event=event,
                                    eventdate=eventdate,
                                    model=model,
                                    version=version,
                                    patch_type=patch_type,
                                    nested=nested,
                                    time_model=PatchVersion._buildTimeModel(time_model),
                                    criteria_overrides=copy.copy(overrides)))
                            version=value
                            continue
                        elif version is None:
                            raise RuntimeError("{0} cannot appear before Version".format(command))

                        if command == 'patchtype':
                            tvalue=value.lower().split()
                            if len(tvalue) not in (1,2):
                                raise RuntimeError("Invalid PatchType {0} - must be just one or two values"
                                               .format(value))
                            for v in tvalue:
                                if v not in ('forward','reverse','hybrid','none'):
                                    raise RuntimeError("Invalid PatchType - must be forward, reverse, hybrid, or none")
                            patch_type=tvalue

                        elif command == 'subgridmethod':
                            if value.lower() == 'nested':
                                vnested=True
                            elif value.lower() == 'additive':
                                vnested=False
                            else:
                                raise RuntimeError("Invalid SubgridMethod - must be nested or additive")
                            if nested is None:
                                nested=vnested
                            elif nested != vnested:
                                raise RuntimeError("Inconsistent SubgridMethod - all patch versions must use the same method")

                        elif command == 'timemodel':
                            time_model=value.split()
                            in_time_model=True
                        else:
                            found=False
                            for f in PatchGridCriteria.fields:
                                fcmd=f.replace('_','').lower()
                                if command == fcmd:
                                    try:
                                        fvalue=[float(x) for x in value.split()]
                                    except:
                                        raise RuntimeError("Invalid value {0} for {1}".
                                                           format(value,command))

                                    if len(fvalue) == 1:
                                        fvalue=fvalue*2
                                    if len(fvalue) != 2:
                                        raise RuntimeError("Invalid value {0} for {1}".
                                                           format(value,command))
                                    overrides[f]=fvalue
                                    found=True
                                    break
                            if not found:
                                raise RuntimeError("Unrecognized command {0} in {1}".format(command,patchfile))

                    elif in_time_model: # Not a command line
                        time_model.extend(l.split())
                    else:
                        raise RuntimeError("Invalid data \"{0}\" in {1}".format(l,patchfile))
                except Exception as ex:
                    msg1="Error in patch definition file {0} at line {1}".format(patchfile,lineno)
                    msg2=ex.message
                    raise RuntimeError(msg1+"\n"+msg2)
                        
            if version is None:
                raise RuntimeError("No version specified in patch definition file {0}".format(patchfile))

            try:
                model_filename=model_filename or default_model_filename
                model_filename=os.path.join(patchdir,model_filename)
                if model_filename != loaded_model_filename:
                    loaded_model_filename=model_filename
                    model=FaultModel(model_filename)
                patchversions.append(PatchVersion(
                    patchname=patchname,
                    event=event,
                    eventdate=eventdate,
                    model=model,
                    version=version,
                    patch_type=patch_type,
                    nested=nested,
                    time_model=PatchVersion._buildTimeModel(time_model),
                    criteria_overrides=copy.copy(overrides)))
            except Exception as ex:
                raise RuntimeError("Error in {0}: {1}"
                                   .format(patchfile,ex.message))

        return patchversions

    @staticmethod
    def _parseDate( datestr ):
        datestr=datestr.strip()
        for format in ('%d %b %Y','%d %B %Y','%Y-%m-%d'):
            try:
                return datetime.datetime.strptime(datestr,format)
            except:
                pass
        raise ValueError("Cannot parse date {0}".format(datestr))

    @staticmethod
    def _buildTimeModel(time_model):
        if time_model is None:
            return None
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
            model.append(PatchVersion.TimePoint(fdate,factor))

        for i in range(len(model)-1):
            if model[i].date > model[i+1].date:
                raise RuntimeError("Dates out of order in TimeModel")
        return model

    class GridSet(namedtuple('GridSet','version ordinates patch_type key')):

        def grid_criteria( self ):
            ordinates=self.ordinates
            if ordinates == 'vertical':
                criteria=Config.vertical_grid_criteria
                nord=1
            else:
                criteria=Config.horizontal_grid_criteria
                nord=0
            overrides=self.version.criteria_overrides
            if overrides:
                override={k:overrides[k][nord] for k in overrides}
                criteria=criteria.copy().update(**override)
            return criteria

        def is_hybrid( self ):
            return self.patch_type=='hybrid'

    def __init__( self, 
                    patchname=None,
                    model=None,
                    version=None,
                    event=None,
                    eventdate=None,
                    patch_type=None,
                    nested=None,
                    time_model=None,
                    criteria_overrides={} ):
        self.patchname=patchname
        self.model=model
        self.event=event or model.eventname
        self.eventdate=eventdate or model.eventdate
        self.version=version
        self.patch_type=patch_type if patch_type is not None else ['forward']
        self.nested=nested
        self.time_model=(time_model if time_model is not None 
                         else [PatchVersion.TimePoint(self.eventdate,1.0)])
        self.criteria_overrides=criteria_overrides

    def grid_sets( self ):
        '''
        String defining the grids that are required for the patch version.
        Intent is to uniquely identify the grid set in a way that allows 
        common grids to be used for multiple versions.

        May be either one or two sets depending on whether horizontal and
        vertical are processed separately.
        '''
        hv=['3d'] if len(self.patch_type)==1 else ['horizontal','vertical']
        sets=[]
        overrides=self.criteria_overrides
        override_key=':'.join([x+':'+str(overrides[x]) for x in sorted(overrides)])
        modelpath=self.model.modelpath
        for ordinates,type in zip(hv,self.patch_type):
            if type == 'none':
                continue
            ishybrid=(type == 'hybrid')
            isvertical=(ordinates == 'vertical')
            key=':'.join([str(x) for x in (modelpath,isvertical,ishybrid,override_key)])
            sets.append(PatchVersion.GridSet(self,ordinates,type,key))
        return sets
    
    def __str__( self ):
        return "\n".join((
            "Version: {0}".format(self.version),
            "Type: {0}".format(self.patch_type),
            "Event: {0}".format(self.event),
            "Date: {0}".format(self.date.strftime('%Y-%m-%d')),
            "TimeModel: {0}\n  ".format(
                "\n  ".join(("{0} {1:.2f}".format(x[0].strftime('%Y-%m-%d'),x[1]) 
                             for x in self.time_model)))
            ))

class FaultModel( object ):

    def __init__( self, modelfile ):
        self.modelpath=None
        self.modelname=None
        self.eventname=None
        self.eventdate=None
        self.description=None
        self.hybridTolerance=None
        for template in Config.model_filename_templates:
            mf=template.format(modelfile)
            if os.path.exists(mf):
                self.modelpath=mf
                break
        if self.modelpath is None:
            raise RuntimeError('Cannot find fault model file {0}'.format(modelfile))
        self.loadModelFile()


    def loadModelFile( self ):
        '''
        Loads a fault model definition
        '''
        
        self.modelname=os.path.splitext(os.path.basename(self.modelpath))[0]
        with open(self.modelpath) as f:
            event = f.readline()
            model = f.readline()
            version = f.readline()
            description = event+' '+model+' '+version
            match=re.search(r'^\s*(\S.*?)\s*(\d{1,2})\s+(\w{3})\w*\s+(\d{4})\s*$',event)
            eventdate = None
            if match:
                eventname=match.group(1)
                try:
                    datestr=match.group(2)+' '+match.group(3)+' '+match.group(4)
                    eventdate = datetime.datetime.strptime(datestr,'%d %b %Y')
                except:
                    pass
            if not eventdate:
                raise RuntimeError("Event record at start of {0} does not end with a valid date"
                                   .format(modelfile))
            self.description=description
            self.eventname=eventname
            self.eventdate=eventdate

    def _cacheKey( self, xy ):
        key="{0:.10f}\t{1:.10f}".format(float(xy[0]),float(xy[1]))
        key=re.sub(r'(\.\d+?)0*(\t|$)',r'\1\2',key)
        return key

    def _cacheFile( self ):
        return Config.filename(name=self.modelname,extension='_grid.cache')

    def _loadValuesFromCache( self, keys ):
        values=[None]*len(keys)
        keyi={k:i for i,k in enumerate(keys)}
        keyre=re.compile('^[^\t]+\t[^\t]+')
        header=None
        if os.path.exists(self._cacheFile()):
            with open(self._cacheFile()) as cf:
                header=cf.next()
                for l in cf:
                    match=keyre.match(l)
                    if match:
                        k=match.group()
                        if k in keyi:
                            values[keyi[k]]=l
        return values

    def _cacheHeader( self ):
        header=None
        if os.path.exists(self._cacheFile()):
            with open(self._cacheFile()) as cf:
                header=cf.next()
        return header

    def _loadGridFromCache( self, spec, gridfile ):
        Logger.write("_loadGridFromCache {0}".format(spec))
        keys=[]
        for n in spec.nodes():
            keys.append(self._cacheKey(n))

        values=self._loadValuesFromCache( keys )
        
        if os.path.exists(gridfile):
            os.remove(gridfile)

        missing=[k for k,v in zip(keys,values) if v is None]
        if missing:
            Logger.write("{0} missing values".format(len(missing)))
            return missing

        header = self._cacheHeader()
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
            if os.path.exists(missingfile):
                os.path.remove(missingfile)
            if os.path.exists(calcfile):
                os.path.remove(calcfile)
            with open(missingfile,'w') as mf:
                for l in missing[imissing:imissing+blocksize]:
                    mf.write(l)
                    mf.write('\n')

            modelpath=self.modelpath
            params=[calc_okada,'-f','-x','-l','-s','-t',modelpath,missingfile,calcfile]
            Logger.write(" ".join(params),1)
            okada_output=check_output(params)
            Logger.write(okada_output,1)

            header=self._cacheHeader()
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
        Logger.write("Running calcgrid ..{0}".format(gridspecstr),0)
        modelpath=self.modelpath
        # if verbose:
        #     print "Calculating deformation on {0}".format(gridfile)
        #     print "Model {0}".format(modeldef)
        #     print "Grid spec {0}".format(gridspec)

        params=[calc_okada,'-f','-x','-l','-s','-t',modelpath,gridspecstr,gridfile]
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
        if built:
            print "          Using cached file {0}".format(gridfile)
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
    '''
    PatchGridDef defines a specific grid within a patch, or used during the calculation of a grid patch.
    May either be directly based on a fault model, or sourced from a derived grid file.  When it is 
    evaluated from a fault model the source is replaced with the calculated file.

    Parameters are:

        name         The base name of the file
        subpatch     The subpatch name used to derive the file name (passed on to grid splitting routine)
        extents      The extents over which the grid is required in order to meet the
                     grid criteria resolution requirements
        buffered_extents  The extents over which the grid may differ from its parent grid - buffered to
                     allow a smooth transition from the parent grid
        level        The level of the grid from 1 (coarsest) to n
        spec         The GridSpec object defining the grid extents and resolution.  Required for fault model
                     based grids.  May be overridden for file based grids
        source       Either a function for evaluating the grid based upon the grid specification, or
                     the name of a source datafile.
        refgrid      A grid used to define a reference used for merging a grid into its parent grid
    '''
                                 
    def __init__( self, name=None, subpatch=None, extents=None, buffered_extents=None, level=1, 
                 spec=None, parent=None, source=None, refgrid=None ):
        self.level=level
        self.parent=parent
        self.name=name
        self.subpatch=subpatch
        self._extents=extents
        self._extentsWkt=None
        self._bufferedExtents=buffered_extents or extents
        self._bufferedExtentsWkt=None
        self._spec=spec
        self._grid=None
        self._strainGrid=None
        self._csv=None
        self.gdf=None
        self.base=None
        self.basename=None
        self.refgrid=refgrid
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
            parent=parent if parent is not Unspecified else None,
            spec=self._spec.copy() if self._spec is not None else None
            )
        pgd.isforward=isforward if isforward is not Unspecified else self.isforward
        pgd.base=self.base
        pgd.basename=self.basename
        pgd.refgrid=self.refgrid
        return pgd

    def _isFileSource( self ):
        return isinstance(self.source,basestring)

    def setSource( self, source ):
        self.source=source
        self.grid=None
        if self._isFileSource():
            self._spec=None

    def setRefGrid( self, refgrid ):
        self.refgrid=refgrid

    @property
    def filename( self ):
        return self.fileWithExtension('.grid')

    @property
    def spec( self ):
        if self._spec is None:
            if self._isFileSource():
                grid=self._loadGrid()
                grid_params=grid.grid_params()
                self._spec=GridSpec(*grid_params)
        return self._spec

    def _loadGrid( self, build_only=False ):
        if self.source is None:
            raise RuntimeError('PatchGridDef._loadGrid: Grid source not defined')
        gridfile=self.filename
        source=self.source
        if not self._isFileSource():
            if self._spec is None:
                raise RuntimeError('{0} grid spec not defined before calculating grid values'
                                   .format(self.name))
            source(self._spec, gridfile )
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
        if self._grid is None:
            self._grid=self._loadGrid()
        return self._grid

    @property
    def csvfile( self ):
        global gridtool
        if self._csv is None:
            gridfile=self.builtFilename
            csvfile = re.sub(r'(\.grid)?$','.csv',gridfile,re.I)
            self.writeCsvFile( csvfile, ndp=-1 )
            self._csv=csvfile
        return self._csv

    def writeCsvFile( self, filename, columns=['de','du','dn'], ndp=4, 
                     subgrid=None, trim=False ):
        global gridtool
        if self._csv is None or filename != self._csv:
            gridfile=self.builtFilename
            csvfile = re.sub(r'(\.grid)?$','.csv',gridfile,re.I)
            colnames='+'.join(columns)
            commands=[gridtool,
                      'read','maxcols','3',gridfile,
                     ]
            if subgrid:
                commands.extend([
                      'trimto','extents',
                       str(subgrid.xmin),str(subgrid.ymin),
                       str(subgrid.xmax),str(subgrid.ymax)
                    ])
            if trim:
                commands.extend([
                    'trim','tolerance','0.00001','leave','1'
                    ])
            commands.extend([
                      'write','csv','columns',colnames,'ndp',str(ndp),'dos',filename,
                      'read','csv',filename
                ])
            Logger.write(" ".join(commands),1)
            gt_output=check_output(commands)
            Logger.write(gt_output,1)

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
            Logger.write(" ".join(commands),1)
            gt_output=check_output(commands)
            Logger.write(gt_output,1)
            self._gdf=gdffile
        return self._gdf

    def fileWithExtension( self, extension ):
        return Config.filename(name=self.name,extension=extension)

    def setExtents( self, extents, bufferedExtents ):
        self._extents=extents
        self._bufferedExtents=bufferedExtents
        self._extentsWkt=None
        self._bufferedExtentsWkt=None

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
            extentfile=self.fileWithExtension('_extent.wkt')
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
            extentfile=self.fileWithExtension('_buffer.wkt')
            with open(extentfile,'w') as wktf:
                for pgn in self.bufferedExtents:
                    wktf.write(pgn.to_wkt())
                    wktf.write('\n')
            self._bufferedExtentsWkt=extentfile
        return self._bufferedExtentsWkt

    def mergeIntoParent( self, saveext=None ):
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
        savefile=None
        if saveext:
            savefile=self.fileWithExtension(saveext+'.grid')
        
        commands=[gridtool,
                  'read','maxcols','3',gridfile,
                 ]
        if savefile:
            commands.extend([
                  'write',savefile
                  ])
        commands.extend([
                  'trimto','buffer','1','wkt',bufwkt,
                  'trimto',pgridfile,
                  'zero','outside',extwkt,'and','edge','1',
                  'add','maxcols','3',pgridfile,'where','outside',extwkt,'and','edge','1',
                 ])
        if self.refgrid is not None:
            commands.extend([
                  'subtract','maxcols','3',self.refgrid
                  ])
        commands.extend([
                  'smooth','linear','outside',extwkt,'not','outside',bufwkt,'not','edge','1',
                  ])
        if self.refgrid is not None:
            commands.extend([
                  'add','maxcols','3',self.refgrid
                  ])
        commands.extend([
                  'write',gridfile,
                  'read',gridfile,
                 ])
        Logger.write(" ".join(commands),1)
        gt_output=check_output(commands)
        Logger.write(gt_output,1)
        self.setSource(gridfile)

    def maxShiftOutsideAreas( self, areas ):
        grid=self.grid
        total_area = prep(areas.buffer(0.0001))
        dscol = grid.getIndex('ds')
        return max((n[dscol] for n in grid.nodes() 
                     if not total_area.contains(Point(n[0],n[1]))))

    def regionsExceedingLevel( self, column, tolerance, absolute=False ):
        grid=self.grid
        if column not in grid.columns:
            if self._strainGrid is None:
                self._strainGrid=grid.strainComponents()
            grid=self._strainGrid
        return MultiPolygon(grid.regionsExceedingLevel( column, tolerance, absolute=absolute ))

    def areasNeedingFinerGrid( self, tolerance, celldimension, vertical ):
        '''
        Determine the areas needing a finer grid than cellsize based
        on a required tolerance
        '''
        grid=self.grid
        gridr = grid.calcResolution(tolerance,precision=0.0000001,vertical=vertical)
        subcell_areas=gridr.regionsLessThanLevel('reqsize',celldimension )
        if len(subcell_areas) <= 0:
            return None
        return subcell_areas

    # Recursive function (recursing on grid level) that populates the array 
    # gridlist with PatchGridDef objects.

    def splitToSubgrids( self, grid_criteria ):
        griddef=self
        subgrids=PatchGridList(self)
        Logger.write("Building level {0} patch {1}".format(self.level,self.name),self.level)

        subcell_areas=[]
        if self.level < grid_criteria.max_split_level:
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
                grid_criteria.grid_level_split_criteria, resolution, grid_criteria.subsplit_using_vertical)
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
                        grid_criteria.grid_level_split_sea_criteria,
                        resolution, grid_criteria.subsplit_using_vertical)
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
                    s.setRefGrid( grid2.builtFilename)
                    subgrids.extend(s.splitToSubgrids(grid_criteria))
            
        else:
            resolution=self.spec.cellDimension()
            subcell_areas=griddef.areasNeedingFinerGrid( 
                grid_criteria.grid_level_split_criteria, resolution, grid_criteria.subsplit_using_vertical)
            if subcell_areas is not None:
                subcell_areas=Util.asMultiPolygon(subcell_areas)
                Logger.writeWkt("{0} not meeting accuracy requirements".format(self.name),self.level,subcell_areas.to_wkt())

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

    def splitGridSpecOnSize(self,nsplit,tolerance=0.5):
        specs=self.spec.splitGridOnSize(nsplit,tolerance)
        if len(specs) > 1:
            buffer=self.bufferedExtents
            specs=[s for s in specs if buffer.intersects(s.boundingPolygon())]
        return specs

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
            "    Grid spec: {0}".format(self._spec),
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

    def __init__( self, *grids ):
        list.__init__(self)
        self.extend(grids)

    def remove( self, grid ):
        for g in self:
            if g.parent == grid:
                g.parent=grid.parent
        grid.parent=None
        list.remove(self,grid)

    def removeRedundantGrids( self ):
        # Remove grids which are completely overlapped by child grid.  Note:
        # should this set the child grid extents to the union with the parent
        # grid extents
        global gridtool
        gridlist=sorted(self,key=lambda x: x.level)
        removed=PatchGridList()
        for g in gridlist:
            p=g.parent
            if p is None:
                continue
            mspec=g.spec.mergeWith(p.spec)
            ratio=float(mspec.nodeCount())/float(g.spec.nodeCount())
            skip=False
            if ratio > Config.redundantParentRatio:
                continue
            Logger.write("Grid {0} candidate for replacement by child grid"
                         .format(p.name))
            for other in gridlist:
                if other != g and other.parent == p:
                    skip=True
                    break
            if skip:
                Logger.write("    Grid retained as has other children")
                continue
            savename=g.fileWithExtension("_unexpanded.grid")
            commands=[
                gridtool,
                'read',g.builtFilename,
                'write',savename,
                'expandto',p.builtFilename,
                'replace',p.builtFilename,
                'replace',g.builtFilename,
                'write',g.builtFilename
                ]
            extents=g.extents.union(p.extents)
            buffered=g.bufferedExtents.union(p.bufferedExtents)
            g.setSource(g.builtFilename)
            g.setExtents(extents,buffered)
            self.remove(p)
            removed.append(p)
            Logger.write("    grid {0} displaced by {1} (size increased by {2:.2f})"
                         .format(p.name,g.name,ratio))
        return removed

def split_forward_reverse_patches( trialgrid, grid_criteria, gridlist ):
    hybrid_tol=grid_criteria.forward_patch_max_distortion
    Logger.write("Splitting patch grids into forward/reverse patches")
    Logger.write("Forward patch tolerated distortion (ppm): {0}".format(hybrid_tol))
    reverse_extents = trialgrid.regionsExceedingLevel(grid_criteria.forward_patch_test_column,hybrid_tol)
    Logger.write("{0} regions found".format(len(reverse_extents)))
    Logger.writeWkt("Extents requiring reverse patch (base tolerance)",0,reverse_extents.to_wkt())
    reverse_grid_extents=[g.bufferedExtents
                          for g in gridlist if g.level > grid_criteria.forward_patch_max_level]
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

        Logger.write("Processing {0}".format(g.name))

        if (g.parent in not_reverse 
                or not g.bufferedExtents.intersects(reverse_extents)):
            if g.level > grid_criteria.forward_patch_max_level:
                raise RuntimeError('Forward patch required for level greater than {0}'
                                   .format(grid_criteria.forward_patch_max_level))
            Logger.write('Copying {0} as forward patch grid only'.format(g.name))
            not_reverse.add(g)
            fgridname=g.name+'_F'
            Logger.writeWkt("Outside forward grid {0} before update".format(g.name),0,g.spec.boundingPolygonWkt())
            gf=g.updatedWith(isforward=True,name=fgridname,parent=forward_grids.get(g.parent,None),
                             source=g.builtFilename)
            Logger.writeWkt("Outside forward grid {0}_F before merge".format(gf.name),0,gf.spec.boundingPolygonWkt())
            #print "Split: forward gf: {0}".format(gf)
            #gf.mergeIntoParent()
            #Logger.writeWkt("Outside forward grid {0}_F after merge".format(gf.name),0,gf.spec.boundingPolygonWkt())
            #print "Split: forward gf merged: {0}".format(gf)
            forward_grids[g]=gf
            splitgrids.append(gf)
            continue

        # Create smoothed grid for forward patch, smoothing over reverse patch
        # extents
        #
        # Have tried different smoothing options (including datumgrid), but
        # none seems to offer particular advantage.  Could be an area for more
        # research

        if g.level <= grid_criteria.forward_patch_max_level:
            Logger.write('Creating smoothed forward patch grid for {0}'.format(g.name))
            fgridname=g.name+'_F'
            basefg=forward_grids.get(g.parent,None)

            Logger.writeWkt("Overlapping forward grid {0} before update".format(g.name),0,g.spec.boundingPolygonWkt())
            gf=g.updatedWith(name=fgridname,isforward=True,source=None,
                            parent=basefg)

            fgridfile=gf.filename

            commands=[
                gridtool,
                'read','maxcols','3',g.builtFilename
                #,'write','csv',g.fileWithExtension('_presmooth.csv')
                ]
            if basefg is not None:
                commands.extend(['replace','maxcols','3',basefg.builtFilename,'where','inside',reversewkt])
            commands.extend(['smooth','linear','inside',reversewkt])
            if basefg is not None:
                commands.extend(['not','on_grid',basefg.builtFilename])
            commands.extend(['write',fgridfile,'read',fgridfile])
            #commands.extend(['write','csv',g.fileWithExtension('_postsmooth.csv')])
        
            Logger.write(" ".join(commands),1)
            gt_output=check_output(commands)
            Logger.write(gt_output,1)

            if not os.path.exists(fgridfile):
                raise RuntimeError('Cannot build smoothed forward patch '+fgridfile)
            gf.setSource(fgridfile)

            #gf.mergeIntoParent()
            Logger.writeWkt("Overlapping forward grid {0}_F".format(gf.name),0,gf.spec.boundingPolygonWkt())
            Logger.write('Created forward patch for {0}'.format(fgridfile))

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

        if not g.extents.intersects(reverse_extents):
            continue
        Logger.write('Creating reverse patch grid for {0}'.format(g.name))
        buffersize=g.spec.cellDimension()*Config.cell_split_factor
        r_extents=Util.asMultiPolygon(g.extents.intersection(reverse_extents));
        r_buffered_extents=Util.asMultiPolygon(
            g.bufferedExtents.intersection(reverse_extents).union(
                Util.bufferedPolygon(r_extents,buffersize)
            ));

        Logger.writeWkt("Constructing reverse grid for {0}".format(g.name),0,g.spec.boundingPolygonWkt())
        rgridname=g.name+'_R'
        rgrid=g.updatedWith(
            extents=r_extents,
            buffered_extents=r_buffered_extents,
            name=rgridname,
            parent=reverse_grids.get(g.parent,None),
            isforward=False,
             )
        fgrid=forward_grids[g]
        rgridfile=rgrid.filename
        refgridfile=rgrid.fileWithExtension('_ref.grid')
        commands=[
            gridtool,
            'read','maxcols','3',g.builtFilename,
            'subtract',fgrid.builtFilename,
            'trimto',fgrid.builtFilename,
            'write',refgridfile,
            'zero','edge','1',
            'trim','tolerance','0.00001','1',
            'write',rgridfile,
            'read',rgridfile]
        Logger.write(" ".join(commands),1)
        gt_output=check_output(commands)
        Logger.write(gt_output,1)
        if not os.path.exists(rgridfile):
            raise RuntimeError('Cannot build reverse patch '+rgridfile)
        rgrid.setSource(rgridfile)
        rgrid.setRefGrid(refgridfile)
        Logger.write('Created reverse patch for {0}'.format(g.name))
        Logger.writeWkt("Reverse grid {0}".format(rgrid.name),0,g.spec.boundingPolygonWkt())
        rbase=reverse_grids.get(g,None)
        reverse_grids[g]=rgrid
        splitgrids.append(rgrid)

    Logger.write("\nSplit grids")
    for g in splitgrids:
        Logger.write("-------------------\n{0}".format(g))
    
    return splitgrids

def build_deformation_grid_set( grid_set, subpatch='', splitbase=True ):
    '''
    Function builds all the grids that are used.  Basically sets up a trial grid to determine
    the extents on which patches are required, then for each patch required within the extents
    calls create_grids to build the actual patches.  create_grids is a recursive function that 
    will break the grids down further if the grid resolution compromises the accuracy requirements

    The returned list of grid definitions includes a field called base.  If splitbase is true
    then this will identify the separate patch areas, otherwise it will be a single patch area.
    '''

    # Calculate a trial grid to determine extents of patches

    
    grid_criteria=grid_set.grid_criteria()
    gridsource=grid_set.version.model.calcGrid

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

    real_extents = trialgrid.regionsExceedingLevel(grid_criteria.base_limit_test_column,
                                                   grid_criteria.base_limit_tolerance)
    Logger.write("{0} patch extents found".format(len(real_extents)))
    Logger.writeWkt("Extents requiring patch (base tolerance)",0,real_extents.to_wkt())

    # Bounds beyond which patch not required because guaranteed by local accuracy bounds
    bounds_extents = trialgrid.regionsExceedingLevel(
        grid_criteria.base_limit_bounds_column,
        grid_criteria.base_limit_bounds_tolerance,absolute=True)
    Logger.writeWkt("Absolute bounds on patch",0,bounds_extents.to_wkt())

    real_extents=real_extents.intersection( bounds_extents )
    real_extents=Util.asMultiPolygon(real_extents)
    Logger.writeWkt("Intersected extents",0,real_extents.to_wkt())

    # Now want to find maximum shift outside extents of test.. 
    # Prepare a multipolygon for testing containment.. add a buffer to
    # handle points on edge of region where patch runs against base polygon

    dsmax = trialgrid.maxShiftOutsideAreas( real_extents )
    buffersize = dsmax/grid_criteria.ramp_tolerance
    Logger.write("Maximum shift outside patch extents {0}".format(dsmax))
    Logger.write("Buffering patches by {0}".format(buffersize))

    buffered_extents=Util.bufferedPolygon( real_extents, buffersize )
    buffered_extents=buffered_extents.intersection( bounds_extents )
    buffered_extents=Util.asMultiPolygon(buffered_extents)
    Logger.writeWkt("Buffered extents",0,buffered_extents.to_wkt())

    # Test areas against land definitions, buffer, and merge overlapping grid definitions

    if Config.land_areas:
        Logger.write("Intersecting areas with land extents".format(len(real_extents)))
        Logger.writeWkt("Land areas",0,Config.land_areas.to_wkt())
        real_extents=real_extents.intersection(Config.land_areas)
        real_extents=Util.asMultiPolygon(real_extents)
        buffered_extents=buffered_extents.intersection( Config.land_areas )
        buffered_extents=Util.asMultiPolygon(buffered_extents)
        Logger.writeWkt("Bounds after land area intersection",0,real_extents.to_wkt())
        Logger.writeWkt("Buffered_bounds after land area intersection",0,buffered_extents.to_wkt())

        # If limiting to land areas then also need to calculate sea areas where
        # lower tolerance applies

        sea_extents = trialgrid.regionsExceedingLevel(
            grid_criteria.base_limit_test_column,
            grid_criteria.base_limit_sea_tolerance)
        if len(sea_extents) == 0:
            Logger.write("No potential sea extents found")
        else:
            Logger.write("{0} patch sea extents found".format(len(sea_extents)))
            Logger.writeWkt("Sea extents requiring patch (base tolerance)",0,sea_extents.to_wkt())

            # Sea bounds beyond which patch not required because guaranteed by local accuracy bounds
            sea_bounds_extents = trialgrid.regionsExceedingLevel(
                grid_criteria.base_limit_bounds_column,
                grid_criteria.base_limit_sea_bounds_tolerance,absolute=True)

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
            sea_buffersize = sea_dsmax/grid_criteria.ramp_sea_tolerance
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
                                        subpatch=subpatch,
                                        source=gridsource,
                                       )

    # Now split the grids into subgrids to meet grid resolution accuracy criteria

    gridlist=PatchGridList()
    for griddef in griddefs:
        patchgrids=griddef.splitToSubgrids(grid_criteria)
        for p in patchgrids:
            p.base=griddef.name if splitbase else Config.patchname
            p.basename=Config.patchname
            gridlist.append(p)

    # Merge grids into parent grids

    for g in gridlist:
        g.mergeIntoParent(saveext='_unmerged')

    #for g in gridlist:
    #    Logger.writeWkt("{0} bounds".format(g.name),0,g.spec.boundingPolygonWkt())
    #sys.exit()

    # If splitting into a hybrid model then separate out forward and reverse patches.  
    # This also merges grids into parents.
    
    if grid_set.is_hybrid():
        Logger.dumpGridList("Grids before hybrid/merge",gridlist)
        gridlist=split_forward_reverse_patches( trialgrid, grid_criteria, gridlist )
        Logger.dumpGridList("Grids after hybrid/merge",gridlist)

    # Optimize grids.  Remove parent grids completely overlapped by child grid.
    # This could be expanded to include previous implementation which removed levels which
    # are not significantly more extensive than child levels.
    #
    # Also could consider splitting grid level into multiple grids so that redundant areas of
    # the grid can be discarded (when much of the grid is outside the buffered extents over 
    # which it differs from its parent.

    for g in gridlist.removeRedundantGrids():
        Logger.write("Removing redundant grid {0}".format(g.name))
    return gridlist

def build_published_component( patchversions, built_gridsets, 
                              comppath, cleandir=False, ndp=4,
                              splitsize=None, splittol=0.5 ):
    '''
    Creates the component.csv file and grid components used as a published
    component of a LINZ published deformation model
    '''

    lastversion=patchversions[-1]
    modeldesc=lastversion.model.description
    eventdate=lastversion.eventdate
    patchdir='patch_'+lastversion.patchname

    # If the model date doesn't contain a date, then append it
    if not re.search(r'[12]\d\d\d[01]\d[0123]\d',patchdir):
        patchdir=patchdir+'_'+eventdate.strftime('%Y%m%d')

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

    Logger.write("Splitting large grids to approx size {0} rows/cols tolerance {1}" 
                 .format(splitsize,splittol))

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
            time0=eventdate,
            factor0=0,
            time1=eventdate,
            factor1=1,
            decay=0,
            file1='',
            description=modeldesc
        )
        
        ir0=1
        for nv in range(len(patchversions)):
            versionspec=patchversions[nv]
            csvdata['version_added']=versionspec.version
            csvdata['version_revoked']=0
            if nv < len(patchversions)-1:
                csvdata['version_revoked']=patchversions[nv+1].version
            for gs in versionspec.grid_sets():
                gridlist=built_gridsets[gs.key]
                toplevelgrids=[g for g in gridlist if g.isTopLevelGrid()]
                ordinates=gs.ordinates
                columns=( ['de','dn'] if ordinates == 'horizontal'
                         else ['du'] if ordinates == 'vertical'
                         else ['de','dn','du'] )
                patch_type=gs.patch_type
                for tg in toplevelgrids:
                    ir1 = ir0
                    ir0 += len(versionspec.time_model)
                    priority=0
                    Logger.write("Publishing grids for top level grid {0}".format(tg.name))
                    for grid in gridlist:
                        if grid.topLevelGrid() != tg:
                            continue
                        rvs=not grid.isforward if patch_type == 'hybrid' else (patch_type != 'forward')
                        csvdata['reverse_patch']='Y' if rvs else 'N'
                        csvname=os.path.split(grid.builtFilename)[1]
                        csvname=os.path.splitext(csvname)[0]
                        csvname = re.sub(r'^(grid_)?','grid_',csvname)
                        subgrids=grid.splitGridSpecOnSize(splitsize,splittol)
                        if len(subgrids) > 1:
                            Logger.write("Splitting grid into {0} smaller grids"
                                         .format(len(subgrids)))
                            csvname=csvname+"_{0:02d}"
                        Logger.writeWkt('Published '+csvname+' bounds',grid.level,
                            grid.bufferedExtents.wkt)
                        for isubgrid, subgrid in enumerate(subgrids):
                            subcsvname=csvname.format(isubgrid)+".csv"
                            subcsvpath=os.path.join(comppath,subcsvname)
                            Logger.write("Building published grid {0} for extents {1}"
                                            .format(subcsvname,subgrid))
                            grid.writeCsvFile( subcsvpath, columns, ndp, subgrid=subgrid, trim=True )
                            priority += 1
                            gd=defgrid.DeformationGrid(subcsvpath)
                            min_lon=round(np.min(gd.column('lon')),10)
                            max_lon=round(np.max(gd.column('lon')),10)
                            min_lat=round(np.min(gd.column('lat')),10)
                            max_lat=round(np.max(gd.column('lat')),10)
                            disp2=None
                            for c in columns:
                                col=gd.column(c)
                                disp2=col*col if disp2 is None else disp2+col*col
                            compdata=dict(
                                min_lon=min_lon,
                                max_lon=max_lon,
                                min_lat=min_lat,
                                max_lat=max_lat,
                                npoints1=gd.array.shape[1],
                                npoints2=gd.array.shape[0],
                                max_displacement=round(math.sqrt(np.max(disp2)),5),
                                file1=subcsvname,
                                displacement_type=gs.ordinates,
                                )
                            Logger.writeWkt('Published '+subcsvname,grid.level,
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
    parser.add_argument('patch_file',help='Model file(s) used to calculate deformation, passed to calc_okada')
    parser.add_argument('calc_file_root',help='Base name used for calculation files')
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
    parser.add_argument('--published-grid-size',type=int,default=100,help="Published grids will be split into subgrids of about this size")
    # parser.add_argument('--no-trim-subgrids',action='store_false',help="Subgrids will not have redundant rows/columns trimmed")
    parser.add_argument('--precision',type=int,default=5,help="Precision (ndp) of output grid displacements")
    parser.add_argument('--land-area',help="WKT file containing area land area over which model must be defined")
    # parser.add_argument('--parcel-shift',action='store_true',help="Configure for calculating parcel_shift rather than rigorous deformation patch")
    # parser.add_argument('--test-settings',action='store_true',help="Configure for testing - generate lower accuracy grids")
    parser.add_argument('--write-grid-lines',action='store_true',help="xx.grid.wkt contains grid cell lines rather than polygon")
    parser.add_argument('-v','--verbose',action='store_true',help="More verbose output (if not already configured)")

    args=parser.parse_args()

    Config.setPatchFilenameRoot(args.calc_file_root)
        
    split_base=args.split_base
#    shift_path=args.shift_model_path
    comp_path=args.submodel_path
#    if comp_path and shift_path:
#        raise RuntimeError("Cannot create published patch and shift model")
    additive=args.subgrids_additive
    if additive:
        raise RuntimeError("Non-nested grids not supported in current implementation")
#    trimgrid=args.no_trim_subgrids
#    if args.parcel_shift: 
#        Logger.write("Configuring model to use parcel shift parameters")
#        args.configureForParcelShift()
#    if args.test_settings: 
#        Logger.write("Configuring model with low accuracy testing settings")
#        args.configureForTesting()
#    if comp_path: 
#        args.apply_time_model_scale=False

    if args.max_level:
        Config.horizontal_grid_criteria.update(max_split_level=args.max_level)
        Config.vertical_grid_criteria.update(max_split_level=args.max_level)

    if args.base_tolerance:
        Config.horizontal_grid_criteria.update(base_limit_tolerance=args.base_tolerance)
        Config.vertical_grid_criteria.update(base_limit_tolerance=args.base_tolerance)

    Config.published_grid_split_size=args.published_grid_size
        
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
                    
        patchversions=PatchVersion.loadPatchDefinition(args.patch_file)

        gridsets={}
        subpatches=[]
        for v in patchversions:
            for gs in v.grid_sets():
                if gs.key not in gridsets:
                    subpatch=''
                    if gs.ordinates != '3d':
                        subpatch=gs.ordinates[0].upper()
                    sp=subpatch
                    while subpatch in subpatches:
                        nsp += 1
                        subpatch="{0}S{1}".format(sp,nsp)
                    subpatches.append(sp)
                    gridsets[gs.key]=build_deformation_grid_set( gs, subpatch=subpatch, splitbase=split_base )

        if comp_path:
            build_published_component( patchversions, gridsets, comp_path, 
                                      args.clean_dir, args.precision,
                                      Config.published_grid_split_size,
                                      Config.published_grid_split_tol
                                     )

    except Exception as ex:
        Logger.write('\n\nFailed with error: {0}'.format(ex.message),level=-1)
        raise

