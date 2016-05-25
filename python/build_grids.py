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
from subprocess import call
import accuracy_standards
import affinity
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

class PatchGridDef:

    def __init__( self, level, file, parent=None, extentfile=None ):
        self.level=level
        self.file=file
        self.name=re.sub(r'\.grid$','',os.path.basename(file),re.I)
        self.parent=parent
        self.extentfile=extentfile
        self.csv=None
        self.gdf=None

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
subcell_ramp_tolerance=accuracy_standards.NRF.lap*tolerance_factor

# Replace entire grid with subcells if more than this percentage 
# of area requires splitting

max_subcell_ratio=0.5

# Maximum split level

max_split_level=5

# Alternative values for creating a Landonline parcel shifting patch
def configure_for_parcel_shift():
    global base_limit_test_column, base_limit_tolerance, base_ramp_tolerance
    base_limit_test_column='ds'
    base_limit_tolerance= 0.05   
    base_ramp_tolerance=accuracy_standards.BGN.lap*tolerance_factor

# Reduce the accuracies to produce smaller data sets for testing...
def configure_for_testing():
    global base_limit_test_column, base_limit_tolerance, base_ramp_tolerance
    global subcell_resolution_tolerance, subcell_ramp_tolerance
    global max_split_level
    subcell_resolution_tolerance *= 100
    subcell_ramp_tolerance *= 100
    max_split_level -= 2

# Apply scale factor and convergence

apply_sf_conv=True

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

published_version=1
publish_reverse_patch=True
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
    file2
    description
    '''.split()

# Log file ...

runtime=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
runtime_version=datetime.datetime.now().strftime("%Y%m%d")
logfile=None
wktfile=None
wktgridfile=None

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

def calc_grid( modeldef, griddef, gridfile ):
    global calc_okada
    gridspec=grid_def_spec( griddef )
    # if verbose:
    #     print "Calculating deformation on {0}".format(gridfile)
    #     print "Model {0}".format(modeldef)
    #     print "Grid spec {0}".format(gridspec)

    params=[calc_okada,'-f','-x','-l','-s',modeldef,gridspec,gridfile]
    if apply_sf_conv:
        params.insert(1,'-f')
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
    return defgrid.defgrid(gridfile)

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

def create_grids( gridlist, modeldef, patchpath, name, level, grid_def, cellsize, grid1=None, parent=None, extentfile=None ):
    write_log("Building level {0} patch {1}".format(level,name),level)
    def1 = grid_def
    if not grid1:
        gridfile = os.path.join(patchpath,"{0}_L{1}.grid".format(name,level))
        grid1 = calc_grid(modeldef,def1,gridfile)

    subcell_areas=[]
    if level < max_split_level:
        # Try calculating subcells if necessary

        subcellsize = list(x/cell_split_factor for x in cellsize)
        def2 = bounds_grid_def( def1, subcellsize, cell_split_factor )
        gridfile2 = os.path.join(patchpath,"{0}_SL{1}.grid".format(name,level+1))
        grid2 = calc_grid(modeldef,def2,gridfile2)

        # Determine required resolution of subcells, and find the extents for which 
        # the subcell resolution is required (ie the required resolution is less than
        # the parent cell size)

        grid2r = grid2.calcResolution(subcell_resolution_tolerance,precision=0.0000001)
        parentsize=cellsize[1]/scale_factor[1]
        subcell_areas=grid2r.regionsExceedingLevel('reqsize',-parentsize,multiple=-1)

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
            subcell_areas.extend(grid2r.regionsExceedingLevel('reqsize',-parentsize,multiple=-1))

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
            create_grids( gridlist, modeldef, patchpath, name, level+1, def2, subcellsize, grid2, parent=parent, extentfile=extentfile )
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
        # supercedes the bounds, as outside this the subcell is merged into the 
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
        create_grids( gridlist, modeldef, patchpath, subcellname, level+1, grid_def, subcellsize, parent=patchdef, extentfile=extentfileroot )


def build_deformation_grids( patchpath, patchname, modeldef, splitbase=True ):
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

    real_extents = trialgrid.regionsExceedingLevel(base_limit_test_column,base_limit_tolerance)
    real_extents=MultiPolygon(real_extents)
    write_log("{0} patch extents found".format(len(real_extents)))
    write_wkt("Extents requiring patch (base tolerance)",0,real_extents.to_wkt())

    # Bounds beyond which patch not required because guaranteed by local accuracy bounds
    bounds_extents = trialgrid.regionsExceedingLevel(base_limit_bounds_column,base_limit_bounds_tolerance)
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

        sea_extents = trialgrid.regionsExceedingLevel(base_limit_test_column,base_limit_sea_tolerance)
        if len(sea_extents) == 0:
            write_log("No potential sea extents found")
        else:
            sea_extents=MultiPolygon(sea_extents)
            write_log("{0} patch sea extents found".format(len(sea_extents)))
            write_wkt("Sea extents requiring patch (base tolerance)",0,sea_extents.to_wkt())

            # Sea bounds beyond which patch not required because guaranteed by local accuracy bounds
            sea_bounds_extents = trialgrid.regionsExceedingLevel(
                base_limit_bounds_column,base_limit_sea_bounds_tolerance)


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
        for p in patchgrids:
            p.base=name if splitbase else patchname
            p.basename=patchname
            gridlist.append(p)
    return gridlist

def create_patch_csv( patchlist, modelname, additive=True, trimgrid=True, linzgrid=False ):
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
    can be turned of to keep grid boundaries aligned with parent grid (makes no
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


def write_log( message, level=0 ):
    prefix='  '*level
    if logfile:
        logfile.write(prefix)
        logfile.write(message)
        logfile.write("\n")
        logfile.flush()
    if verbose:
        print prefix+message
        sys.stdout.flush()

def write_grid_wkt( name, level, griddef ):
    lines = grid_def_lines( griddef )
    for l in lines:
        wktgridfile.write("{0}|{1}|LINESTRING({2} {3},{4} {5})\n".format(
            name,level,l[0][0],l[0][1],l[1][0],l[1][1]))

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

def get_model_spec( modelfile ):
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
    ramp=[]
    rampfile=os.path.splitext(modelfile)[0]+'.ramp'
    if os.path.exists(rampfile):
        ramp=[]
        with open(rampfile) as f:
            for l in f:
                if re.match(r'^\s*\#\s*$',l):
                    continue
                m = re.match(r'^\s*(\d\d\d\d\-\d\d\-\d\d)\s+(\d+(?:\.\d*))\s*$',l)
                if m:
                    ramp.append((m.group(1),float(m.group(2))))
                else:
                    raise ValueError("Invalid record "+l+" in ramp file "+rampfile)
    if not ramp:
        ramp=[(modeldate,1.0)]
    return modelname, description, modeldate, ramp

def build_published_component( gridlist, modeldef, additive, comppath, cleandir=False ):
    '''
    Creates the component.csv file and grid components used as a published
    component of a LINZ published deformation model
    '''

    modelname,modeldesc,modeldate,ramps = get_model_spec( modeldef )
    patchdir='patch_'+modelname
    # If the model date doesn't contain a date, then append it
    if not re.search(r'[12]\d\d\d[01]\d[0123]\d',patchdir):
        patchdir=patchdir+'_'+modeldate.replace('-','')

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

    finalramp=ramps[-1][1]
    compcsv=os.path.join(comppath,'component.csv')
    with open(compcsv,"w") as ccsvf:
        ccsv=csv.writer(ccsvf)
        ccsv.writerow(published_component_columns)
        csvdata=dict(
            version_added=published_version,
            version_revoked=0,
            reverse_patch='Y' if publish_reverse_patch else 'N',
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
            factor0=-1 if publish_reverse_patch else 0,
            time1=modeldate,
            factor1=0 if publish_reverse_patch else 1,
            decay=0,
            file1='',
            file2='',
            description=modeldesc
        )
        for priority, grid in enumerate(gridlist):
            gd=defgrid.defgrid(grid.csv)
            (gridpath,csvname)=os.path.split(grid.csv)
            csvname = re.sub(r'^(grid_)?','grid_',csvname)
            compcsv=os.path.join(comppath,csvname)
            if compcsv != grid.csv:
                shutil.copyfile(grid.csv,compcsv)
            compdata=dict(
                min_lon=np.min(gd.column('lon')),
                max_lon=np.max(gd.column('lon')),
                min_lat=np.min(gd.column('lat')),
                max_lat=np.max(gd.column('lat')),
                npoints1=gd.array.shape[1],
                npoints2=gd.array.shape[0],
                max_displacement=math.sqrt(np.max(
                    gd.column('de')*gd.column('de') +
                    gd.column('dn')*gd.column('dn') +
                    gd.column('du')*gd.column('du'))),
                file1=csvname,
                )
            for ir,r in enumerate(ramps):
                rdate=r[0]
                rvalue=r[1]
                if ir == 0:
                    if rvalue==0:
                        continue
                    compdata['time_function']='step'
                    compdata['time0']=rdate
                else:
                    compdata['time_function']='ramp'
                    compdata['time0']=ramps[ir-1][0]
                    rvalue -= ramps[ir-1][1]
                compdata['time1']=rdate
                compdata['factor0']=-rvalue if publish_reverse_patch else 0
                compdata['factor1']=0 if publish_reverse_patch else rvalue
                if not additive:
                    compdata['component']=ir+1
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
    parser.add_argument('--submodel-version',default=runtime_version,help="Deformation model version for which submodel first applies")
    parser.add_argument('--subgrids-nest',action='store_false',help="Grid CSV files calculated to replace each other rather than total to deformation")
    parser.add_argument('--parcel-shift',action='store_true',help="Configure for calculating parcel_shift rather than rigorous deformation patch")
    parser.add_argument('--apply-ramp-scale',action='store_true',help="Scale the grid by the ramp final value")
    parser.add_argument('--apply-sf-conv',action='store_true',help="Apply the projection scale factor and convergence to the model calculated displacements")
    parser.add_argument('--test-settings',action='store_true',help="Configure for testing - generate lower accuracy grids")
    parser.add_argument('--max-level',type=int,help="Maximum number of split levels to generate (each level increases resolution by 4)")
    parser.add_argument('--reverse',action='store_true',help="Published model will be a reverse patch")
    parser.add_argument('--base-tolerance',type=float,help="Base level tolerance - depends on base column")
    parser.add_argument('--split-base',action='store_true',help="Base deformation will be split into separate patches if possible")
    parser.add_argument('--no-trim-subgrids',action='store_false',help="Subgrids will not have redundant rows/columns trimmed")
    parser.add_argument('--clean-dir',action='store_true',help="Clean publishable component subdirectory")
    parser.add_argument('--land-area',help="WKT file containing area land area over which model must be defined")

    args=parser.parse_args()
        
    patchfile = args.patch_file
    split_base=args.split_base
    shift_path=args.shift_model_path
    comp_path=args.submodel_path
    published_version=args.submodel_version
    publish_reverse_patch=args.reverse
    additive=args.subgrids_nest
    trimgrid=args.no_trim_subgrids
    apply_sf_conv=args.apply_sf_conv
    if args.parcel_shift: 
        write_log("Configuring model to use parcel shift parameters")
        configure_for_parcel_shift()
    if args.test_settings: 
        write_log("Configuring model with low accuracy testing settings")
        configure_for_testing()
    max_split_level=args.max_level if args.max_level else max_split_level
    base_limit_tolerance=args.base_tolerance if args.base_tolerance else base_limit_tolerance

    patchpath,patchname=os.path.split(args.patch_file)
    if not os.path.isdir(patchpath):
        os.makedirs(patchpath)

    if args.land_area:
        load_land_areas(args.land_area)
                
    models=args.model_file
    moddefs=[]
    if not models:
        models.append(patchname)

    # Check model files are valid

    for i,f in enumerate(models):
        found = False
        for template in ('{0}','{0}.model','model/{0}','model/{0}.model'):
            mf = template.format(f) 
            if os.path.exists(mf): 
                found = True
                models[i]=mf
                if args.apply_ramp_scale:
                    modelname,modeldesc,modeldate,ramps = get_model_spec( mf )
                    rampscale=ramps[-1][1]
                    if rampscale != 1.0:
                        mf=str(rampscale)+'*'+mf
                moddefs.append(mf)
                break
                
        if not found:
            print "Model file {0} doesn't exist".format(f)
            sys.exit()

    # Model definition for calc okada
    modeldef='+'.join(moddefs)
    modelname=' '.join([os.path.basename(x) for x in models])

    logfile=open(patchfile+".build.log","w")
    wktfile=open(patchfile+".build.wkt","w")
    wktgridfile=open(patchfile+".grids.wkt","w")
    logfile.write("Building deformation grids for model: {0}\n".format(modeldef))
    wktfile.write("description|level|wkt\n")
    wktgridfile.write("name|level|wkt\n")
    
    gridlist = build_deformation_grids( patchpath, patchname, modeldef, splitbase=split_base )

    build_linzgrid=bool(shift_path)

    create_patch_csv( gridlist, modelname, linzgrid=build_linzgrid, additive=additive, trimgrid=trimgrid )

    if shift_path:
        build_linzshift_model( gridlist, modeldef, additive, shiftpath=shift_path, splitbase=split_base )

    if comp_path:
        build_published_component( gridlist, modeldef, additive, comp_path, args.clean_dir )
