#!/usr/bin/python
# Script to build a set of nested grids to meet tolerances etc.
#
# Not a pretty script, if this was being used a lot it could be refactored to use a lot more
# object type stuff (ie not a good example of python code!).  And broken up
# into smaller bits!
#

from collections import namedtuple
from shapely.geometry import MultiPolygon,Polygon,LineString,Point
from shapely.prepared import prep
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
# and to which grids are trimmed)
base_extents=((166,-47.5),(176,-39.5))

# Conversion degrees to metres
scale_factor=(1.0e-5/math.cos(math.radians(40)),1.0e-5)

# Factor by which accuracy specs are multiplied to provide grid tolerances
tolerance_factor=0.4

# Values for deformation patch generation
# Test column for determining extents of patch
base_limit_test_column='err'

# Minium value for test column within patch (note - column is in ppm)
base_limit_tolerance=accuracy_standards.NRF.lap*tolerance_factor*1.0e6

# Precision for scaling down from edge of patch
base_ramp_tolerance=accuracy_standards.NRF.lap*tolerance_factor

def configure_for_parcel_shift():
    # Alternative values for creating a Landonline parcel shifting patch
    global base_limit_test_column, base_limit_tolerance, base_ramp_tolerance
    base_limit_test_column='ds'
    base_limit_tolerance= 0.05   
    base_ramp_tolerance=accuracy_standards.BGN.lap*tolerance_factor

# Values controlling splitting cells into subgrids

cell_split_factor=4
subcell_resolution_tolerance=accuracy_standards.NRF.lac*tolerance_factor
subcell_ramp_tolerance=accuracy_standards.NRF.lap*tolerance_factor

# Replace entire grid with subcells if more than this percentage 
# of area requires splitting

max_subcell_ratio=0.5

# Maximum split level

max_split_level=5

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
NEGATIVE yes
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
    subcomponent
    priority
    min_lon
    max_lon
    min_lat
    max_lat
    spatial_complete
    min_date
    max_date
    temporal_complete
    npoints1
    npoints2
    displacement_type
    error_type
    max_displacement
    spatial_model
    temporal_model
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
    gridspec=grid_def_spec( griddef )
    # if verbose:
    #     print "Calculating deformation on {0}".format(gridfile)
    #     print "Model {0}".format(modeldef)
    #     print "Grid spec {0}".format(gridspec)

    params=['calc_okada','-x','-l','-s',modeldef,gridspec,gridfile]
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
    and scales back.
    '''
    if buffer > 0:
        polygon=affinity.scale(polygon,xfact=1.0/scale_factor[0],yfact=1.0/scale_factor[1],origin=origin)
        polygon=polygon.buffer(buffer)
        polygon=affinity.scale(polygon,xfact=scale_factor[0],yfact=scale_factor[1],origin=origin)
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
        grid2r=None


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
                if extentpoly.contains(area):
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
    write_log("{0} patch extents found".format(len(real_extents)))

    # Now want to find maximum shift outside extents of test.. 
    # Prepare a multipolygon for testing containment.. add a small buffer to
    # handle points on edge of region where patch runs against base polygon

    dsmax = maxShiftOutsideAreas( trialgrid, real_extents )
    buffersize = dsmax/base_ramp_tolerance

    write_log("Maximum shift outside patch extents {0}".format(dsmax))
    write_log("Buffering patches by {0}".format(buffersize))

    # Properly should apply same merging of grid definitions here as in create_grid (or
    # refactor so that they all come together...)

    name = patchname
    gridlist=[]
    for i,extent in enumerate(real_extents):
        if len(real_extents) > 1:
            name='{0}_P{1}'.format(patchname,i)
        write_log("Building patch {0}".format(name))
        write_wkt("Patch base extents {0}".format(name),0,extent.to_wkt())
        buffered_extent=buffered_polygon( extent, buffersize )
        write_wkt("Buffered base extents {0}".format(name),0,buffered_extent.to_wkt())
        bounds = expanded_bounds( buffered_extent, cellbuffer=base_size )
        write_wkt("Expanded extents {0}".format(name),0,bounds_wkt(bounds))
        griddef = bounds_grid_def( bounds, base_size )
        write_wkt("Grid extents {0}".format(name),0,bounds_wkt(griddef))
        extfile=os.path.join(patchpath,name)
        with open(extfile+".extent.wkt","w") as f:
            f.write(extent.to_wkt())
            f.write("\n")
        with open(extfile+".buffer.wkt","w") as f:
            f.write(buffered_extent.to_wkt())
            f.write("\n")
        patchgrids=[]
        create_grids(patchgrids, modeldef,patchpath,name,1,griddef,base_size,extentfile=extfile)
        for p in patchgrids:
            p.base=name if splitbase else patchname
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
        headers = dict(
            header1="Patch component "+component,
            header2="Model source "+modelname,
            header3=runtime
        )

        # Build a grid tool command to process the file
        
        # Read in the existing file (only need de, dn, du components)
        merge_commands='gridtool read maxcols 3 gridfile '

        # If have a parent, then subtract it as want to merge into it
        if parent: merge_commands = merge_commands + ' subtract csv parentcsv'

        # Set to zero outside area extents for which patch subcell resolution is required
        merge_commands=merge_commands + ' zero outside extents'

        # Ensure zero at edge
        merge_commands=merge_commands + ' zero edge 1'

        # Smooth out to buffer to avoid sharp transition
        merge_commands=merge_commands + ' smooth linear outside extents not outside buffer not edge 1'

        # Trim redundant zeros around the edge of the patch
        merge_commands=merge_commands + ' trim 1'

        # If LINZ grid required then write it out.  LINZ grid is always additive
        # so this comes before adding the parent command back on
        if linzgrid:
            merge_commands=merge_commands + ' write_linzgrid NZGD2000 header1 header2 header3 gdffile'

        # If we have a parent and are not making additive patches then add the parent back
        # on.
        if parent and not additive: merge_commands = merge_commands + ' add csv parentcsv'

        # Finally write out as a CSV file
        merge_commands=merge_commands + ' write csv csvfile'

        # Run the commands
        gridtoolcmd=[ 
            gridfile if x == 'gridfile' else
            csvfile if x == 'csvfile' else
            gdffile if x == 'gdffile' else
            parent.csv if x == 'parentcsv' else
            extfile if x == 'extents' else
            buffile if x == 'buffer' else
            headers[x] if x in headers else
            x
            for x in merge_commands.split()
            ]

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
    if verbose:
        print prefix+message

def write_grid_wkt( name, level, griddef ):
    lines = grid_def_lines( griddef )
    for l in lines:
        wktgridfile.write("{0}|{1}|LINESTRING({2} {3},{4} {5})\n".format(
            name,level,l[0][0],l[0][1],l[1][0],l[1][1]))

def write_wkt( item, level, wkt ):
    wktfile.write("{0}|{1}|{2}\n".format(item,level,wkt))

def build_linzshift_model( gridlist, shiftpath=None ):
    '''
    Creates shift definition files using list of grids.  Assumes that grids have
    been generated using linzgrid=True
    '''

    if not gridlist or not gridlist[0].gdf:
        raise RuntimeError("Cannot build shift model - haven't built LINZ grid format files")
    if not shiftpath:
        shiftpath = os.path.dirname(gridlist[0].gdf)
    if not os.path.isdir(shiftpath):
        os.makedirs(shiftpath)
    baselist=set([patch.base for patch in gridlist])
    for base in baselist:
        shiftdef=base+'_shift.def'
        shiftbin=base+'_shift.bin'
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

def get_model_spec( modeldef ):
    name=''
    date=None
    descriptions=''
    for modelfile in modeldef.split('+'):
        modelname=os.path.basename(modelfile)
        modelname=re.sub(r'\..*','',modelname)
        name = name+modelname
        with open(modelfile) as f:
            event = f.readline()
            model = f.readline()
            version = f.readline()
            descriptions = descriptions+event+model+version
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
            if not date:
                date = modeldate
            if modeldate != date:
                raise RuntimeError("Events forming model don't have a common date")
    return name, descriptions.strip(), date

def build_published_component( gridlist, modeldef, additive, comppath=None ):
    '''
    Creates the component.csv file and grid components used as a published
    component of a LINZ published deformation model
    '''

    modelname,modeldesc,modeldate = get_model_spec( modeldef )
    patchdir='patch_'+modelname+'_'+modeldate.replace('-','')

    if not gridlist or not gridlist[0].csv:
        raise RuntimeError("Cannot build shift model - haven't built grid CSV files")
    if not comppath:
        comppath = os.path.dirname(gridlist[0].gdf)
    comppath = os.path.join( comppath, 'model', patchdir )
    if not os.path.isdir(comppath):
        os.makedirs(comppath)

    write_log("Writing published model component {0}".format(comppath))

    compcsv=os.path.join(comppath,'component.csv')
    with open(compcsv,"w") as ccsvf:
        ccsv=csv.writer(ccsvf)
        ccsv.writerow(published_component_columns)
        csvdata=dict(
            version_added=published_version,
            version_revoked=0,
            reverse_patch='Y' if publish_reverse_patch else 'N',
            subcomponent=1,
            priority=1,
            min_lon=0,
            max_lon=0,
            min_lat=0,
            max_lat=0,
            spatial_complete='Y',
            min_date=0,
            max_date=0,
            temporal_complete='Y',
            npoints1=0,
            npoints2=0,
            displacement_type='3d',
            error_type='none',
            max_displacement=0,
            spatial_model='llgrid',
            temporal_model='step',
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
                npoints1=gd.array.shape[0],
                npoints2=gd.array.shape[1],
                max_displacement=np.max(
                    gd.column('de')*gd.column('de') +
                    gd.column('dn')*gd.column('dn') +
                    gd.column('du')*gd.column('du')),
                file1=csvname,
                )
            if not additive:
                compdata['priority']=priority
            csvdata.update(compdata)
            ccsv.writerow([csvdata[c] for c in published_component_columns])

if __name__ == "__main__":

    # Process arguments

    parser=argparse.ArgumentParser(description='Build set of grids for deformation patch')
    parser.add_argument('patch_file',help='Base name used for output files')
    parser.add_argument('model_file',help='Model file(s) used to calculate deformation, passed to calc_okada',nargs='+')
    parser.add_argument('--shift-model-path',help="Create a linzshiftmodel in the specified directory")
    parser.add_argument('--component-path',help="Create publishable component in the specified directory")
    parser.add_argument('--component-model-version',default=runtime_version,help="Deformation model version to include in published component")
    parser.add_argument('--subgrids-override',action='store_false',help="Grid CSV files calculated to replace each other rather than total to deformation")
    parser.add_argument('--parcel-shift',action='store_true',help="Configure for calculating parcel_shift rather than rigorous deformation patch")
    parser.add_argument('--max-level',type=int,default=max_split_level,help="Maximum number of split levels to generate (each level increases resolution by 4)")
    parser.add_argument('--split-base',action='store_true',help="Base deformation will be split into separate patches if possible")
    parser.add_argument('--no-trim-subgrids',action='store_false',help="Subgrids will not have redundant rows/columns trimmed")

    args=parser.parse_args()
        
    patchfile = args.patch_file
    max_split_level=args.max_level
    split_base=args.split_base
    shift_path=args.shift_model_path
    comp_path=args.component_path
    published_version=args.component_model_version
    additive=args.subgrids_override
    trimgrid=args.no_trim_subgrids
    if args.parcel_shift: 
        write_log("Configuring model to use parcel shift parameters")
        configure_for_parcel_shift()

    patchpath,patchname=os.path.split(args.patch_file)
    if not os.path.isdir(patchpath):
        os.makedirs(patchpath)
                
    models=args.model_file
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
                break
        if not found:
            print "Model file {0} doesn't exist".format(f)
            sys.exit()

    # Model definition for calc okada
    modeldef='+'.join(models)
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
        build_linzshift_model( gridlist, shiftpath=shift_path )

    if comp_path:
        build_published_component( gridlist, modeldef, additive, comp_path )
