#!/usr/bin/python
# Script to build a set of nested grids to meet tolerances etc.
#
# Not a pretty script, if this was being used a lot it could be refactored to use a lot more
# object type stuff (ie not a good example of python code!).  And broken up
# into smaller bits!
#

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
import weakref

calc_okada=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'okada','calc_okada')
if not os.path.exists(calc_okada):
    raise RuntimeError('Cannot find calc_okada program at '+calc_okada)

gridtool=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'gridtool','gridtool')
if not os.path.exists(gridtool):
    raise RuntimeError('Cannot find gridtool program at '+gridtool)

class PatchGridDef:

    def __init__( self, level=0, file=None, parent=None, extentfile=None ):
        self.level=level
        self.file=file
        self.name=re.sub(r'\.grid$','',os.path.basename(file),re.I)
        self.parent=parent
        self.extentfile=extentfile
        self.csv=None
        self.gdf=None
        self.base=None
        self.basename=None
        self.isforward=None  # Component is forward/reverse in hybrid grid else None (= either)
        self._extents=None 
        self._buffered_extents=None 

    def update( self, level=None, file=None, parent=None, extentfile=None, isforward=None ):
        pgd=PatchGridDef(
            level=level if level is not None else self.level,
            file=file if file is not None else self.file,
            extentfile=extentfile if extentfile is not None else self.extentfile,
            parent=parent)
        pgd.isforward=isforward
        return pgd

    def grid( self ):
        return defgrid.DeformationGrid(self.file)

    def extents( self ):
        if self.extentfile is None:
            return None
        if self._extents is None:
            areas=[]
            with open(self.extentfile+'.extent.wkt') as wktf:
                for l in wktf:
                    if l:
                        areas.append(wkt_load(l))
            self._extents= MultiPolygon(areas)
        return self._extents

    def bufferedExtents( self ):
        if self.extentfile is None:
            return None
        if self._buffered_extents is None:
            with open(self.extentfile+'.buffer.wkt') as wktf:
                self._buffered_extents=wkt_load(wktf.read())
        return self._buffered_extents

    def isTopLevelGrid( self ):
        return self.parent is None

    def topLevelGrid( self ):
        return self if self.isTopLevelGrid() else self.parent.topLevelGrid()

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

class PatchVersionDef:

    TimePoint=namedtuple('TimePoint','date factor')

    def __init__( self, version, event, date, reverse=False, nested=False, time_model=None, hybrid=False, hybrid_tol=10.0 ):
        self.version=version
        self.event=event
        self.date=date
        self.reverse=reverse
        self.nested=nested
        self.hybrid=hybrid
        self.hybrid_tol=hybrid_tol
        if time_model is None:
            time_model=[PatchVersionDef.TimePoint(date,1.0)]
        self.time_model=time_model

    
    def __str__( self ):
        return "\n".join((
            "Version: {0}".format(self.version),
            "Event: {0}".format(self.event),
            "Date: {0}".format(self.date.strftime('%Y-%m-%d')),
            "Reverse: {0}".format(self.reverse),
            "Nested: {0}".format(self.nested),
            "Hybrid: {0}".format(self.hybrid),
            "Forward patch max distortion: {0} ppm".format(self.hybrid_tol),
            "TimeModel: {0}\n  ".format(
                "\n  ".join(("{0} {1:.2f}".format(x[0].strftime('%Y-%m-%d'),x[1]) 
                             for x in self.time_model)))
            ))

verbose=True
    
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

# Alternative values for creating a Landonline parcel shifting patch
def configure_for_parcel_shift():
    global base_limit_test_column, base_limit_tolerance, base_ramp_tolerance
    base_limit_test_column='ds'
    base_limit_tolerance= 0.05   
    base_ramp_tolerance=accuracy_standards.BGN.lap*tolerance_factor

# Reduce the accuracies to produce smaller data sets for testing...
def configure_for_testing():
    global base_limit_test_column, base_limit_tolerance, base_ramp_tolerance
    global subcell_resolution_tolerance
    global max_split_level
    subcell_resolution_tolerance *= 100
    max_split_level -= 2


def log_configuration():
    write_log("Patch grid calculation parameters:")
    write_log("    Grid origin: {0} {1}".format(*origin))
    write_log("    Base grid cell size: {0} {1}".format(*base_size))
    write_log("    Base extents: {0} {1} to {2} {3}".format(
        base_extents[0][0], base_extents[0][1], base_extents[1][0], base_extents[1][1]))
    write_log("    Approx degrees to metres factor: {0} {1}".format(*scale_factor))

    write_log("    Tolerated distortion outside patch (ppm): land {0} sea {1}".format( 
        base_limit_tolerance, base_limit_sea_tolerance))
    write_log("    Tolerated dislocation outside patch (mm): land {0} sea {1}".format( 
        base_limit_bounds_tolerance*1000, base_limit_sea_bounds_tolerance*1000))
    write_log("    Patch ramp tolerance (ppm): land {0} sea {1}".format(
        base_ramp_tolerance,base_ramp_sea_tolerance))
    write_log("    Tolerated gridding error before splitting (mm)".format(
        subcell_resolution_tolerance*1000, subcell_resolution_sea_tolerance*1000))
    write_log("    Grid cell split factor: {0}".format(cell_split_factor))
    write_log("    Maximum split level: {0}".format(max_split_level))
    write_log("    Maximum percent coverage of level with nested grid: {0}".format(max_subcell_ratio*100)) 
    write_log("")

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

# Log file ...

runtime=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
runtime_version=datetime.datetime.now().strftime("%Y%m%d")
logfile=None
wktfile=None
wktgridfile=None
wktgridlines=False

def grid_spec( size, extents, multiple=1 ):
    '''
    Calculate grid definition (ie min/max lat/lon, number of cells)
    based on grid size (lon,lat cell size), required extents, and 
    an optional integer.  Grids are aligned using integer multiples
    of size from the origin.  

    Returns grid string definition used by calc_okada
    '''
    def calc_min_max( base, minv, maxv, csize ):
        csize *= multiple
        n0 = int(math.floor((minv-base)/csize))
        n1 = int(math.ceil((maxv-base)/csize))
        return base+n0*csize,base+n1*csize,(n1-n0)*multiple
    lnmin,lnmax,nln=calc_min_max(origin[0],extents[0][0],extents[1][0],size[0])
    ltmin,ltmax,nlt=calc_min_max(origin[1],extents[0][1],extents[1][1],size[1])
    return 'grid:{0}:{1}:{2}:{3}:{4}:{5}'.format(lnmin,ltmin,lnmax,ltmax,nln,nlt)

def calc_grid( modeldef, griddef, gridfile, subtract_grids=[] ):
    global calc_okada
    gridspec=grid_def_spec( griddef )
    # if verbose:
    #     print "Calculating deformation on {0}".format(gridfile)
    #     print "Model {0}".format(modeldef)
    #     print "Grid spec {0}".format(gridspec)

    params=[calc_okada,'-f','-x','-l','-s',modeldef,gridspec,gridfile]
    meta='\n'.join(params)
    metafile=gridfile+'.metadata'
    built = False
    # Check if grid file is already built - mainly convenience for
    # script development
    if os.path.exists(gridfile) and os.path.exists(metafile):
        with open(metafile) as f:
            oldmeta=f.read()
            if oldmeta.strip() == meta.strip():
                built=True
        if os.path.getmtime(gridfile) > os.path.getmtime(metafile):
            built=False
        if built:
            for mf in modeldef.split('+'):
                if '*' in mf:
                    mf=mf.split('*',1)[1]
                if os.path.getmtime(mf) > os.path.getmtime(gridfile):
                    built=False
                    break
    if not built:
        npts=griddef[2][0]*griddef[2][1]
        print "          Building {0} ({1} points) ...".format(gridfile,npts)
        call(params)
        if not os.path.exists(gridfile):
            raise RuntimeError("Failed to calculate grid file {0}".format(gridfile))
        with open(metafile,'w') as f:
            f.write(meta)
    print "          Loading {0} ...".format(gridfile)
    grid=defgrid.DeformationGrid(gridfile)
    if subtract_grids is not None and len(subtract_grids) > 0:
        raise NotImplementedError("subtract_grids not implemented in calc_grid!")
    return grid

def bounds_grid_def( bounds, cellsize, multiple=1 ):
    '''
    Calculate grid definition (ie min/max lat/lon, number of cells)
    based on grid size (lon,lat cell size), required extents, and 
    an optional integer.  Grids are aligned using integer multiples
    of size from the origin.  

    Returns grid string definition used by calc_okada
    '''
    def calc_min_max( base, minv, maxv, csize ):
        csize *= multiple
        n0 = int(math.floor((minv-base)/csize))
        n1 = int(math.ceil((maxv-base)/csize))
        return base+n0*csize,base+n1*csize,(n1-n0)*multiple
    lnmin,lnmax,nln=calc_min_max(origin[0],bounds[0][0],bounds[1][0],cellsize[0])
    ltmin,ltmax,nlt=calc_min_max(origin[1],bounds[0][1],bounds[1][1],cellsize[1])
    return ((lnmin,ltmin),(lnmax,ltmax),(nln,nlt))

def grid_def_spec( griddef ):
    return 'grid:{0}:{1}:{2}:{3}:{4}:{5}'.format(
        griddef[0][0], griddef[0][1],
        griddef[1][0], griddef[1][1],
        griddef[2][0], griddef[2][1],
        )

def grid_def_lines( griddef ):
    lonv = list(griddef[0][0]+(griddef[1][0]-griddef[0][0])*float(x)/griddef[2][0]
            for x in range(griddef[2][0]+1))
    latv = list(griddef[0][1]+(griddef[1][1]-griddef[0][1])*float(x)/griddef[2][1]
            for x in range(griddef[2][1]+1))
    for lat in latv:
        for i in range(len(lonv)-1):
            yield ((lonv[i],lat),(lonv[i+1],lat))
    for lon in lonv:
        for i in range(len(latv)-1):
            yield ((lon,latv[i]),(lon,latv[i+1]))

def grid_def_merge( def1, def2 ):
    pt0 = (min(def1[0][0],def2[0][0]),min(def1[0][1],def2[0][1]))
    pt1 = (max(def1[1][0],def2[1][0]),max(def1[1][1],def2[1][1]))
    dln = ((def1[1][0]-def1[0][0])+(def2[1][0]-def2[0][0]))/(def1[2][0]+def2[2][0])
    dlt = ((def1[1][1]-def1[0][1])+(def2[1][1]-def2[0][1]))/(def1[2][1]+def2[2][1])
    nln= int(math.floor(0.5+(pt1[0]-pt0[0])/dln))
    nlt= int(math.floor(0.5+(pt1[1]-pt0[1])/dlt))
    return (pt0,pt1,(nln,nlt))

def grid_def_list_merge( deflist, griddef ):
    while True:
        overlap = None
        for def1 in deflist:
            if bounds_overlap(def1,griddef):
                overlap=def1
                break
        if overlap:
            deflist.remove(overlap)
            griddef=grid_def_merge(griddef,overlap)
        else:
            break
    deflist.append(griddef)

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

def expanded_bounds( polygon, buffer=0, cellbuffer=(0,0), trim=base_extents ):
    '''
    Generate bounded extents around a polygon.
    buffer is buffer extents in metres
    cellbuffer is (lon,lat) buffer size 
    trim is region to which buffered extents are trimmed
    '''
    bounds=buffered_polygon(polygon,buffer).bounds
    return (
        (
            max(bounds[0]-cellbuffer[0],trim[0][0]),
            max(bounds[1]-cellbuffer[1],trim[0][1])
        ),
        (
            min(bounds[2]+cellbuffer[0],trim[1][0]),
            min(bounds[3]+cellbuffer[1],trim[1][1])
        )
    )

def bounds_poly( bounds ):
    bl=bounds[0]
    br=(bounds[0][0],bounds[1][1])
    tr=bounds[1]
    tl=(bounds[1][0],bounds[0][1])
    return Polygon((bl,br,tr,tl,bl))

def bounds_wkt( bounds ):
    return bounds_poly( bounds ).to_wkt()

def bounds_area( bounds ):
    return (bounds[1][0]-bounds[0][0])*(bounds[1][1]-bounds[0][1])

def bounds_overlap( bounds1, bounds2 ):
    return not (
        bounds1[0][0] > bounds2[1][0] or
        bounds2[0][0] > bounds1[1][0] or
        bounds1[0][1] > bounds2[1][1] or
        bounds2[0][1] > bounds1[1][1])


def maxShiftOutsideAreas( grid, areas ):
    total_area = prep(MultiPolygon(areas).buffer(0.0001))
    dscol = grid.getIndex('ds')
    return max((n[dscol] for n in grid.nodes() 
                 if not total_area.contains(Point(n[0],n[1]))))

# Calculate regions exceeding level.  Assume this is based on a strain component
# column which may be present in source grid, if from calc_grid directly, otherwise
# calculated from it.  Use a weak reference to allow caching the calculated grid.

# Cached strain grid components
source_grid_ref=None
source_grid_strain_grid=None

def regionsExceedingLevel( grid, column, tolerance ):
    global source_grid_ref, source_grid_strain
    if column not in grid.columns:
        if source_grid_ref is None or source_grid_ref == grid:
            source_grid_ref=weakref.ref(grid)
            source_grid_strain=grid.strainComponents()
        grid=source_grid_strain
    return grid.regionsExceedingLevel( column, tolerance )

# Recursive function (recursing on grid level) that populates the array 
# gridlist with PatchGridDef objects.

def create_grids( gridlist, modeldef, patchpath, name, level, grid_def, cellsize, grid1=None, parent=None, extentfile=None, subtract_grids=[] ):
    write_log("Building level {0} patch {1}".format(level,name),level)
    def1 = grid_def
    if not grid1:
        gridfile = os.path.join(patchpath,"{0}_L{1}.grid".format(name,level))
        grid1 = calc_grid(modeldef,def1,gridfile, subtract_grids)

    subcell_areas=[]
    if level < max_split_level:
        # Try calculating subcells if necessary

        subcellsize = list(x/cell_split_factor for x in cellsize)
        def2 = bounds_grid_def( def1, subcellsize, cell_split_factor )
        gridfile2 = os.path.join(patchpath,"{0}_SL{1}.grid".format(name,level+1))
        grid2 = calc_grid(modeldef,def2,gridfile2,subtract_grids)

        # Determine required resolution of subcells, and find the extents for which 
        # the subcell resolution is required (ie the required resolution is less than
        # the parent cell size)

        grid2r = grid2.calcResolution(subcell_resolution_tolerance,precision=0.0000001)
        parentsize=cellsize[1]/scale_factor[1]
        subcell_areas=grid2r.regionsLessThanLevel('reqsize',parentsize)

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
        
            grid2r = grid2.calcResolution(subcell_resolution_sea_tolerance,precision=0.0000001)
            subcell_areas.extend(grid2r.regionsLessThanLevel('reqsize',parentsize))

        grid2r=None
        # Expand subcell areas by parent grid size to avoid rounding issues
        if subcell_areas:
            expanded_area=buffered_polygon(unary_union(subcell_areas),parentsize)
            if expanded_area.type=='Polygon':
                subcell_areas = [expanded_area]
            else:
                subcell_areas = list(expanded_area.geoms)

    subcell_grid_defs=[]
    total_area = 0.0

    if subcell_areas:
        write_log("{0} areas identified requiring subcells".format(len(subcell_areas)),level)

        # Expand the areas up to grid sizes that will contain them

        areabuffersize=max(x/scale_factor[i]*2 for i,x in enumerate(subcellsize))
        write_log("Buffering areas by {0}".format(areabuffersize),level)
        
        for i, area in enumerate(subcell_areas):
            write_wkt("{0} subcell area {1}".format(name,i+1),level+1,area.to_wkt()) 
            area_bounds = expanded_bounds( area, buffer=areabuffersize, trim=def1 )
            write_wkt("{0} subcell expanded_bounds {1}".format(name,i+1),level+1,bounds_wkt(area_bounds)) 
            area_grid_def = bounds_grid_def( area_bounds, subcellsize, cell_split_factor )
            write_wkt("{0} subcell grid area {1}".format(name,i+1),level+1,bounds_wkt(area_grid_def)) 
            grid_def_list_merge( subcell_grid_defs, area_grid_def )

        if len(subcell_grid_defs) < len(subcell_areas):
            write_log("Merged to {0} subcell grids".format(len(subcell_grid_defs)),level)

        for i,area_grid_def in enumerate(subcell_grid_defs):
            total_area += bounds_area( area_grid_def )
            write_wkt("{0} final grid area {1}".format(name,i+1),level+1,bounds_wkt(area_grid_def)) 

        # If subcell area is bigger than specified ratio of total extents then 
        # replace entire grid

        base_area=bounds_area(def1)
        write_log("Total area of subcells={0} - area of base cell={1}".format(
            total_area,base_area),level)

        if total_area > base_area*max_subcell_ratio:
            write_log("More than {0} of area needs splitting - splitting entire area".format(max_subcell_ratio),level)
            # try:
            #     os.remove(grid1.source)
            # except:
            #     write_log("Failed to remove {0}".format(grid1.source))
            grid1=None
            create_grids( gridlist, modeldef, patchpath, name, level+1, def2, subcellsize, grid2, parent=parent, extentfile=extentfile, subtract_grids=subtract_grids )
            return

    grid1file = grid1.source
    patchdef=PatchGridDef(level,grid1file,parent,extentfile)
    gridlist.append(patchdef)
    write_log("Created {0}".format(grid1file),level)
    write_grid_wkt("{0} grid {1}".format(name,grid1file),level,def1)

    # Now generate the grids for each subcell

    for i,grid_def in enumerate(subcell_grid_defs):
        subcellname=name
        if len(subcell_grid_defs) > 1:
            subcellname="{0}_P{1}".format(name,i+1)
        # Need to keep a record of the extents over which the subcell 
        # exceeds the bounds, as outside this the subcell is merged into the 
        # parent 
        extentfileroot=os.path.join(patchpath,"{0}_L{1}".format(subcellname,level+1))
        subset_extentfile = extentfileroot+".extent.wkt"
        subset_bufferfile = extentfileroot+".buffer.wkt"
        extentpoly=bounds_poly(grid_def)
        areas=[]
        with open(subset_extentfile,"w") as f:
            for area in subcell_areas:
                if extentpoly.contains(area) or extentpoly.intersects(area):
                    areas.append(area)
                    f.write(area.to_wkt())
                    f.write("\n")
        areas = buffered_polygon(MultiPolygon(areas),areabuffersize)
        write_wkt("{0} buffered subcell area {1}".format(name,i+1),level+1,areas.to_wkt()) 
        with open(subset_bufferfile,"w") as f:
            f.write(areas.to_wkt())
            f.write("\n")

        # Now process the subcell
        create_grids( gridlist, modeldef, patchpath, subcellname, level+1, grid_def, subcellsize, parent=patchdef, extentfile=extentfileroot, subtract_grids=subtract_grids )

def grid_file_version( gridfile, version ):
    path,ext=os.path.splitext(gridfile)
    return path+version+ext

def split_forward_reverse_patches( patchpath, trialgrid, hybrid_tol, gridlist ):
    write_log("Splitting patch grids into forward/reverse patches")
    write_log("Forward patch tolerated distortion (ppm): {0}".format(hybrid_tol))

    reverse_extents = regionsExceedingLevel(trialgrid, forward_patch_test_column,hybrid_tol)
    reverse_extents=MultiPolygon(reverse_extents)
    write_log("{0} regions found".format(len(reverse_extents)))
    write_wkt("Extents requiring reverse patch (base tolerance)",0,reverse_extents.to_wkt())
    reversewkt=os.path.join(patchpath,'reverse_patch_extents.wkt')
    with open(reversewkt,'w') as f:
        f.write(reverse_extents.to_wkt())

    #write_log("Grid list:")
    #for g in gridlist:
    #    write_log("{0}".format(g),prefix="    ")

    # First create forward patches. Current implementation will only use top
    # level grid for forward patch.

    splitgrids=[]
    
    for g in gridlist:
        if not g.isTopLevelGrid():
            continue
        family=[gc for gc in gridlist if gc.topLevelGrid() == g]
        family.sort(key=lambda x: (x.level, x.file))

        # If the grid does not intersect the reverse extents then it and
        # all its children belong in the forward patch
        if not g.extents().intersects(reverse_extents):
            fgrids={}
            for gc in family:
                write_log('Copying {0} as forward patch grid only'.format(gc.name))
                gcf=gc.update(isforward=True,parent=fgrids.get(gc.parent,None))
                fgrids[gc]=gcf
                splitgrids.append(gcf)
            continue

        # Create smoothed grid for forward patch, smoothing over reverse patch
        # extents
        #
        # Have tried different smoothing options (including datumgrid), but
        # none seems to offer particular advantage.  Could be an area for more
        # research

        fgridfile=grid_file_version(g.file,'_F')
        write_log('Creating smoothed forward patch grid for {0}'.format(g.name))
        commands=[
            gridtool,
            'read','maxcols','3',g.file,
            'smooth','linear','inside',reversewkt,
            'write',fgridfile
            ]
        call(commands)
        if not os.path.exists(fgridfile):
            raise RuntimeError('Cannot build smoothed forward patch '+fgridfile)
        write_log('Created forward patch for {0}'.format(fgridfile))

        fgrid=g.update(file=fgridfile,isforward=True)
        splitgrids.append(fgrid)

        # Now create the reverse patch files by subtracting the forward patch
        # from them

        rgrids={}
        for gc in family:
            write_log('Creating reverse patch grid for {0}'.format(gc.name))
            rgridfile=grid_file_version(gc.file,'_R')
            commands=[
                gridtool,
                'read','maxcols','3',gc.file,
                'subtract',fgridfile,
                'write',rgridfile]
            call(commands)
            if not os.path.exists(rgridfile):
                raise RuntimeError('Cannot build reverse patch '+rgridfile)
            rgrid=gc.update(file=rgridfile,parent=rgrids.get(gc.parent,None),isforward=False)
            rgrids[gc]=rgrid
            splitgrids.append(rgrid)

    #write_log("\nSplit grids")
    #for g in splitgrids:
    #    write_log("-------------------\n{0}".format(g))
    
    return splitgrids
    #raise NotImplementedError('split_forward_reverse_patches not implemented yet')

def build_deformation_grids( patchpath, patchname, modeldef, splitbase=True, hybrid_tol=None ):
    '''
    Function builds all the grids that are used.  Basically sets up a trial grid to determine
    the extents on which patches are required, then for each patch required within the extents
    calls create_grids to build the actual patches.  create_grids is a recursive function that 
    will break the grids down further if the grid resolution compromises the accuracy requirements

    The returned list of grid definitions includes a field called base.  If splitbase is true
    then this will identify the separate patch areas, otherwise it will be a single patch area.
    '''

    # Calculate a trial grid to determine extents of patches

    trialgriddef=bounds_grid_def(base_extents,base_size)
    trialgridfile=patchfile+"_trial.grid".format(patchfile)

    write_log("Building trial grid using {0} to {1}".format(
        grid_def_spec(trialgriddef),trialgridfile))
    write_wkt("Trial extents",0,bounds_wkt(base_extents))

    trialgrid = calc_grid( modeldef, trialgriddef, trialgridfile )

    # Determine the extents on which the base tolerance is exceed - returns
    # a list of polygons

    real_extents = regionsExceedingLevel(trialgrid, base_limit_test_column,base_limit_tolerance)
    real_extents=MultiPolygon(real_extents)
    write_log("{0} patch extents found".format(len(real_extents)))
    write_wkt("Extents requiring patch (base tolerance)",0,real_extents.to_wkt())

    # Bounds beyond which patch not required because guaranteed by local accuracy bounds
    bounds_extents = regionsExceedingLevel(trialgrid, base_limit_bounds_column,base_limit_bounds_tolerance)
    bounds_extents=MultiPolygon(bounds_extents)

    # Now want to find maximum shift outside extents of test.. 
    # Prepare a multipolygon for testing containment.. add a buffer to
    # handle points on edge of region where patch runs against base polygon

    dsmax = maxShiftOutsideAreas( trialgrid, real_extents )
    buffersize = dsmax/base_ramp_tolerance
    write_log("Maximum shift outside patch extents {0}".format(dsmax))
    write_log("Buffering patches by {0}".format(buffersize))

    buffered_extent=buffered_polygon( MultiPolygon(real_extents), buffersize )
    write_wkt("Buffered extents",0,buffered_extent.to_wkt())
    write_wkt("Absolute bounds on patch",0,bounds_extents.to_wkt())

    real_extents=buffered_extent.intersection( MultiPolygon(bounds_extents) )
    if 'geoms' not in dir(real_extents):
        real_extents = MultiPolygon([real_extents])
    write_wkt("Bounds before potential land intersection",0,real_extents.to_wkt())


    # Test areas against land definitions, buffer, and merge overlapping grid definitions

    if land_areas:
        write_log("Intersecting areas with land extents".format(len(real_extents)))
        valid_areas=[]
        for sc in real_extents:
            p=sc.intersection(land_areas)
            if 'geoms' not in dir(p):
                p=[p]
            for g in p:
                if type(g) == Polygon:
                    valid_areas.append(g)
        real_extents=MultiPolygon(valid_areas)
        write_wkt("Bounds after land area intersection",0,real_extents.to_wkt())

        # If limiting to land areas then also need to calculate sea areas where
        # lower tolerance applies

        sea_extents = regionsExceedingLevel(trialgrid, base_limit_test_column,base_limit_sea_tolerance)
        if len(sea_extents) == 0:
            write_log("No potential sea extents found")
        else:
            sea_extents=MultiPolygon(sea_extents)
            write_log("{0} patch sea extents found".format(len(sea_extents)))
            write_wkt("Sea extents requiring patch (base tolerance)",0,sea_extents.to_wkt())

            # Sea bounds beyond which patch not required because guaranteed by local accuracy bounds
            sea_bounds_extents = regionsExceedingLevel(
                trialgrid, base_limit_bounds_column,base_limit_sea_bounds_tolerance)


            # Now want to find maximum shift outside extents of test.. 
            # Prepare a multipolygon for testing containment.. add a buffer to
            # handle points on edge of region where patch runs against base polygon

            sea_dsmax = maxShiftOutsideAreas( trialgrid, sea_extents )
            sea_buffersize = sea_dsmax/base_ramp_sea_tolerance
            write_log("Maximum shift outside sea patch extents {0}".format(sea_dsmax))
            write_log("Buffering sea patches by {0}".format(sea_buffersize))

            sea_buffered_extent=buffered_polygon( MultiPolygon(sea_extents), sea_buffersize )
            write_wkt("Sea buffered extents",0,buffered_extent.to_wkt())

            if len(sea_bounds_extents) > 0:
                sea_bounds_extents=MultiPolygon(sea_bounds_extents)
                write_wkt("Sea absolute bounds on patch",0,sea_bounds_extents.to_wkt())
                sea_extents=sea_buffered_extent.intersection( sea_bounds_extents) 
                if 'geoms' not in dir(sea_extents):
                    sea_extents = MultiPolygon([sea_extents])

            write_wkt("Sea bounds before union with land extents",0,sea_extents.to_wkt())

            real_extents=real_extents.union(sea_extents)
            if 'geoms' not in dir(real_extents):
                real_extents=MultiPolygon([real_extents])
            write_wkt("Bounds after union with sea extents",0,real_extents.to_wkt())

    # Form merged buffered areas
    extent_defs=[]
    extents=[]
    for i, extent in enumerate(real_extents):
        write_wkt("Base required area {0}".format(i+1),0,extent.to_wkt()) 
        bounds = expanded_bounds( extent, cellbuffer=base_size )
        write_wkt("Expanded extents {0}".format(i+1),0,bounds_wkt(bounds))
        griddef = bounds_grid_def( bounds, base_size )
        grid_def_list_merge( extent_defs, griddef )
        extents.append(extent)
    extent_defs.sort(key=lambda x:x[0][1])

    write_log("{0} patch extents after buffering".format(len(extent_defs)))

    # Now process the extents...
    name = patchname
    gridlist=[]
    for i,griddef in enumerate(extent_defs):
        if len(extent_defs) > 1:
            name='{0}_P{1}'.format(patchname,i)
        write_wkt("Grid extents {0}".format(name),0,bounds_wkt(griddef))
        # Write the wkt definitions of the extents and buffered extents defining
        # the areas.  Note need to recalculate buffer in to ensure resulting buffers
        # are merged properly if they overlap
        extfile=os.path.join(patchpath,name)
        extentpoly=bounds_poly(griddef)
        with open(extfile+".extent.wkt","w") as f: 
            areas=[]
            for area in extents:
                if extentpoly.contains(area) or extentpoly.intersects(area):
                    f.write(area.to_wkt())
                    f.write("\n")
                    areas.append(area)

        areas = buffered_polygon(MultiPolygon(areas),buffersize)
        with open(extfile+".buffer.wkt","w") as f:
            f.write(areas.to_wkt())
            f.write("\n")

        patchgrids=[]
        create_grids(patchgrids, modeldef,patchpath,name,1,griddef,base_size,extentfile=extfile)
        if hybrid_tol:
            patchgrids=split_forward_reverse_patches( patchpath, trialgrid, hybrid_tol, patchgrids )
        for p in patchgrids:
            p.base=name if splitbase else patchname
            p.basename=patchname
            gridlist.append(p)
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
    write_log("Creating patch CSV files")
    for patch in patchlist:
        # Get the files used for this component of the patch

        gridfile = patch.file
        component = patch.name
        csvfile = re.sub(r'(\.grid)?$','.csv',gridfile,re.I)
        gdffile = re.sub(r'(\.grid)?$','.gdf',gridfile,re.I)
        write_log("Creating patch file file {0}".format(csvfile))
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
    
        write_log("Running "+
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


def write_log( message, level=0, prefix="" ):
    prefix='  '*level+prefix
    message="\n".join([prefix+m for m in message.split("\n")])
    if logfile:
        logfile.write(message)
        logfile.write("\n")
        logfile.flush()
    if verbose or level < 0:
        print message
        sys.stdout.flush()

def write_grid_wkt( name, level, griddef ):
    ncol=griddef[2][0]
    nrow=griddef[2][1]
    if wktgridlines:
        lines = grid_def_lines( griddef )
        for l in lines:
            wktgridfile.write("{0}|{1}|{2}|{3}|LINESTRING({4} {5},{6} {7})\n".format(
                name,level,ncol,nrow,l[0][0],l[0][1],l[1][0],l[1][1]))
    else:
        points=[
            [griddef[0][0],griddef[0][1]],
            [griddef[0][0],griddef[1][1]],
            [griddef[1][0],griddef[1][1]],
            [griddef[1][0],griddef[0][1]],
            [griddef[0][0],griddef[0][1]],
            ]
        crds=",".join(("{0} {1}".format(*p) for p in points))
        wktgridfile.write("{0}|{1}|{2}|{3}|POLYGON(({4}))\n".format(name,level,ncol,nrow,crds))


def write_wkt( item, level, wkt ):
    wktfile.write("{0}|{1}|{2}\n".format(item,level,wkt))

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
        write_log("Creating shift definition file {0}".format(shiftfile))
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

def _parseDate( datestr ):
    datestr=datestr.strip()
    for format in ('%d %b %Y','%d %B %Y','%Y-%m-%d'):
        try:
            return datetime.datetime.strptime(datestr,format)
        except:
            pass
    raise ValueError("Cannot parse date {0}".format(datestr))

def _buildTimeModel(time_model,modeldate):
    if time_model is None:
        return [PatchVersionDef.TimePoint(_parseDate(modeldate),1.0)]
    model=[]
    ntime=len(time_model)
    if ntime < 2:
        raise RuntimeError("Invalid or empty TimeModel")
    if ntime % 2 != 0:
        raise RuntimeError("TimeModel must consist of paired date and value")
    for i in range(0,ntime,2):
        fdate=_parseDate(time_model[i])
        try:
            factor=float(time_model[i+1])
        except:
            raise RuntimeError("Invalid scale factor {0} in TimeModel"
                               .format(time_model[i+1]))
        model.append(PatchVersionDef.TimePoint(fdate,factor))

    for i in range(len(model)-1):
        if model[i].date > model[i+1].date:
            raise RuntimeError("Dates out of order in TimeModel")
    return model


def get_patch_spec( modelfile ):
    modelname=os.path.basename(modelfile)
    modelname=re.sub(r'\..*','',modelname)
    with open(modelfile) as f:
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
        modeldate=modeldate.strftime('%Y-%m-%d')
    
    patchfile=os.path.splitext(modelfile)[0]+'.patch'
    patchversions=[]
    if os.path.exists(patchfile):
        version=None
        reverse=False
        hybrid=False
        hybrid_tol=default_forward_patch_max_distortion
        nested=False
        time_model=None
        in_model=False

        with open(patchfile) as f:
            for l in f:
                # Ignore comments and blank lines
                if re.match(r'^\s*(\#|$)',l):
                    continue
                cmdmatch=re.match(r'^\s*(\w+)\s*\:\s*(.*?)\s*$',l)
                if cmdmatch:
                    in_model=False
                    command=cmdmatch.group(1).lower()
                    value=cmdmatch.group(2)
                    if command == 'version':
                        if version is not None:
                            try:
                                patchversions.append(PatchVersionDef(
                                    version,
                                    event,
                                    _parseDate(modeldate),
                                    reverse=reverse,
                                    nested=nested,
                                    hybrid=hybrid,
                                    hybrid_tol=hybrid_tol,
                                    time_model=_buildTimeModel(time_model,modeldate)))
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
                        hybrid_tol=float(value)
                    elif command == 'subgridmethod':
                        if value.lower() == 'nested':
                            nested=True
                        elif value.lower() == 'indpendent':
                            nested=False
                        else:
                            raise RuntimeError("Invalid SubgridMethod - must be nested or independent in {0}"
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
            patchversions.append(PatchVersionDef(
                version,
                event,
                _parseDate(modeldate),
                reverse=reverse,
                nested=nested,
                hybrid=hybrid,
                hybrid_tol=hybrid_tol,
                time_model=_buildTimeModel(time_model,modeldate)))
        except Exception as ex:
            raise RuntimeError("Error in {0}: {1}"
                               .format(patchfile,ex.message))

        for i in range(len(patchversions)-1):
            if patchversions[i].version >= patchversions[i+1].version:
                raise RuntimeError("Deformation model versions out of order in {0}".format(patchfile))
            if patchversions[i].nested != patchversions[i+1].nested:
                raise RuntimeError("Inconsistent SubgridMethod option in {0}".format(patchfile))

    return patchversions

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

    write_log("Writing published model submodel {0}".format(comppath))

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
                    compdata=dict(
                        min_lon=round(np.min(gd.column('lon')),10),
                        max_lon=round(np.max(gd.column('lon')),10),
                        min_lat=round(np.min(gd.column('lat')),10),
                        max_lat=round(np.max(gd.column('lat')),10),
                        npoints1=gd.array.shape[1],
                        npoints2=gd.array.shape[0],
                        max_displacement=round(math.sqrt(np.max(
                            gd.column('de')*gd.column('de') +
                            gd.column('dn')*gd.column('dn') +
                            gd.column('du')*gd.column('du'))),5),
                        file1=csvname,
                        )
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

def load_land_areas( polygon_file ):
    global land_areas
    from shapely.wkt import loads
    try:
        with open(polygon_file) as laf:
            wkt=laf.read()
            land_areas=loads(wkt)
        write_log( "Using land area definition from "+polygon_file)
        return
    except:
        raise RuntimeError("Cannot load land area definition from "+polygon_file)

if __name__ == "__main__":

    # Process arguments

    parser=argparse.ArgumentParser(description='Build set of grids for deformation patch')
    parser.add_argument('patch_file',help='Base name used for output files')
    parser.add_argument('model_file',help='Model file(s) used to calculate deformation, passed to calc_okada',nargs='+')
    parser.add_argument('--shift-model-path',help="Create a linzshiftmodel in the specified directory")
    parser.add_argument('--submodel-path',help="Create publishable component in the specified directory")
    parser.add_argument('--clean-dir',action='store_true',help="Clean publishable component subdirectory")
    parser.add_argument('--subgrids-nest',action='store_false',help="Grid CSV files calculated to replace each other rather than total to deformation")
    parser.add_argument('--apply-time-model-scale',action='store_true',help="Scale by the time model final value")
    parser.add_argument('--max-level',type=int,help="Maximum number of split levels to generate (each level increases resolution by 4)")
    parser.add_argument('--base-tolerance',type=float,help="Base level tolerance - depends on base column")
    parser.add_argument('--split-base',action='store_true',help="Base deformation will be split into separate patches if possible")
    parser.add_argument('--no-trim-subgrids',action='store_false',help="Subgrids will not have redundant rows/columns trimmed")
    parser.add_argument('--precision',type=int,default=5,help="Precision (ndp) of output grid displacements")
    parser.add_argument('--land-area',help="WKT file containing area land area over which model must be defined")
    parser.add_argument('--parcel-shift',action='store_true',help="Configure for calculating parcel_shift rather than rigorous deformation patch")
    parser.add_argument('--test-settings',action='store_true',help="Configure for testing - generate lower accuracy grids")
    parser.add_argument('--write-grid-lines',action='store_true',help="xx.grid.wkt contains grid cell lines rather than polygon")

    args=parser.parse_args()
        
    patchfile = args.patch_file
    split_base=args.split_base
    shift_path=args.shift_model_path
    comp_path=args.submodel_path
    if comp_path and shift_path:
        raise RuntimeError("Cannot create published patch and shift model")
    additive=args.subgrids_nest
    trimgrid=args.no_trim_subgrids
    if args.parcel_shift: 
        write_log("Configuring model to use parcel shift parameters")
        configure_for_parcel_shift()
    if args.test_settings: 
        write_log("Configuring model with low accuracy testing settings")
        configure_for_testing()
    if comp_path: 
        args.apply_time_model_scale=False

    max_split_level=args.max_level if args.max_level else max_split_level
    base_limit_tolerance=args.base_tolerance if args.base_tolerance else base_limit_tolerance
    wktgridlines=args.write_grid_lines

    patchpath,patchname=os.path.split(args.patch_file)
    if not os.path.isdir(patchpath):
        os.makedirs(patchpath)

    logfile=open(patchfile+".build.log","w")

    log_configuration()

    try:
        if args.land_area:
            load_land_areas(args.land_area)
                    
        models=args.model_file
        moddefs=[]
        patchspec={}
        if not models:
            models.append(patchname)

        # Check model files are valid

        for i,f in enumerate(models):
            found = False
            for template in ('{0}','{0}.model','fault_models/{0}','fault_models/{0}.model'):
                mf = template.format(f) 
                if os.path.exists(mf): 
                    found = True
                    models[i]=mf
                    if args.apply_time_model_scale:
                        scale=1.0
                        patchversions=get_patch_spec( mf )
                        if patchversions:
                            scale=patchversions[-1].time_model[-1].factor
                        if scale != 1.0:
                            mf=str(scale)+'*'+mf
                    moddefs.append(mf)
                    break
                    
            if not found:
                print "Model file {0} doesn't exist".format(f)
                sys.exit()

        # Model definition for calc okada
        modeldef='+'.join(moddefs)
        modelname=' '.join([os.path.splitext(os.path.basename(x))[0] for x in models])

        hybrid_tol=None
        hybrid_ver=None
        if comp_path: 
            if len(moddefs) > 1:
                raise RuntimeError("Cannot create published patch with more than one fault model")
            patchversions=get_patch_spec( mf )
            isnested=patchversions[0].nested
            additive=not patchversions[0].nested
            for pv in patchversions:
                if pv.nested != isnested:
                    raise RuntimeError(
                        "Inconsistent submodel method in patch versions\n{0}\n-------\n{1}"
                        .format(patchversions[0],pv))
                if pv.hybrid:
                    if hybrid_tol is None:
                        hybrid_tol=pv.hybrid_tol
                        hybrid_ver=pv
                    elif pv.hybrid_tol != hybrid_tol:
                        raise RuntimeError(
                            "Inconsistent hybrid grid tolerance in patch definitions\n{0}\n------\n{1}"
                            .format(hybrid_ver,pv))


        # Initiallize 
        wktfile=open(patchfile+".build.wkt","w")
        wktgridfile=open(patchfile+".grids.wkt","w")
        logfile.write("Building deformation grids for model: {0}\n".format(modeldef))
        wktfile.write("description|level|wkt\n")
        wktgridfile.write("name|level|ncol|nrow|wkt\n")
        
        gridlist = build_deformation_grids( patchpath, patchname, modeldef, splitbase=split_base, hybrid_tol=hybrid_tol )

        build_linzgrid=bool(shift_path)

        create_patch_csv( gridlist, modelname, linzgrid=build_linzgrid, additive=additive, trimgrid=trimgrid, precision=args.precision )

        if shift_path:
            build_linzshift_model( gridlist, modeldef, additive, shiftpath=shift_path, splitbase=split_base )

        if comp_path:
            build_published_component( gridlist, modeldef, modelname, additive, comp_path, args.clean_dir )
    except Exception as ex:
        write_log('\n\nFailed with error: {0}'.format(ex.message),level=-1)
        raise

