#!/usr/bin/python
# Script to build a set of nested grids to meet tolerances etc.

import argparse
import os.path
import math
import sys
import defgrid
from subprocess import call
from shapely.geometry import MultiPolygon,MultiLineString,Polygon,Point
from shapely.prepared import prep
import affinity
import accuracy_standards
from collections import namedtuple
import re

class PatchGridDef:

    def __init__( self, level, file, parent=None, extentfile=None ):
        self.level=level
        self.file=file
        self.parent=parent
        self.extentfile=extentfile
        self.csv=None
    
# Grids will be calculated using integer representation based origin.
# Subgrids will be obtained by dividing basesize by 4

verbose=True

# Origin for calculating grids (integer offsets from this
origin=(166.0,-48)

# Base size for grids, grids are built to this size and smaller
base_size=(0.15,0.125)

# Extents over which test grid is built (ie initial search for affected area,
# and to which grids are trimmed)
base_extents=((166,-47.5),(176,-39.5))

# Conversion degrees to metres
scale_factor=(1.0e-5*math.cos(math.radians(40)),1.0e-5)

# Factor by which accuracy specs are multiplied to provide grid tolerances
tolerance_factor=0.4

# Values for deformation patch generation
# Test column for determining extents of patch
base_limit_test_column='err'
# Minium value for test column within patch
base_limit_tol=accuracy_standards.NRF.lap*tolerance_factor
# Precision for scaling down from edge of patch
base_ramp_precision=accuracy_standards.NRF.lap*tolerance_factor

# Values for parcel block splitting
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

max_split_level=4

# Log file ...

logfile=None
wktfile=None

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
    meta='|'.join(params)
    metafile=gridfile+'.metadata'
    built = False
    # Check if grid file is already built - mainly convenience for
    # script development
    if os.path.exists(gridfile) and os.path.exists(metafile):
        with open(metafile) as f:
            oldmeta=f.read()
            if oldmeta == meta:
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

def grid_def_wkt( griddef ):
    lonv = list(griddef[0][0]+(griddef[1][0]-griddef[0][0])*float(x)/griddef[2][0]
            for x in range(griddef[2][0]+1))
    latv = list(griddef[0][1]+(griddef[1][1]-griddef[0][1])*float(x)/griddef[2][1]
            for x in range(griddef[2][1]+1))
    lines=[]
    for lat in latv:
        lines.append(list((x,lat) for x in lonv))
    for lon in lonv:
        lines.append(list((lon,y) for y in latv))
    return MultiLineString(lines).to_wkt()

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
    total_area = prep(MultiPolygon(real_extents).buffer(0.0001))
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
        # the parent cell size

        grid2r = grid2.calcResolution(subcell_resolution_tolerance,precision=0.0000001)
        parentsize=cellsize[1]/scale_factor[1]
        subcell_areas=grid2r.regionsExceedingLevel('reqsize',-parentsize,multiple=-1)
        grid2r=None


    subcell_grid_defs=[]
    total_area = 0.0
    areabuffersize=max(x/scale_factor[i]*2 for i,x in enumerate(subcellsize))

    if subcell_areas:
        write_log("{0} areas identified requiring subcells".format(len(subcell_areas)),level)

        # Expand the areas up to grid sizes that will contain them

        write_log("Buffering areas by {0}".format(areabuffersize),level)
        
        for i, area in enumerate(subcell_areas):
            write_wkt("{0} subcell area {1}".format(name,i+1),level+1,area.to_wkt()) 
            area_bounds = expanded_bounds( area, buffer=areabuffersize, trim=def1 )
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
    write_wkt("{0} grid {1}".format(name,grid1file),level,grid_def_wkt(def1))

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
        with open(subset_bufferfile,"w") as f:
            f.write(areas.to_wkt())
            f.write("\n")

        # Now process the subcell
        create_grids( gridlist, modeldef, patchpath, subcellname, level+1, grid_def, subcellsize, parent=patchdef, extentfile=extentfileroot )

def create_patch_csv( patchlist, additive=True ):
    # Note: Assumes patchlist is sorted such that parents are before children
    write_log("Creating patch CSV files")
    for patch in patchlist:
        # Get the files used for this component of the patch

        gridfile = patch.file
        csvfile = os.path.basename(gridfile)
        csvfile = re.sub(r'(\.grid)?$','.csv',gridfile,re.I)
        write_log("Creating patch file file {0}".format(csvfile))
        parent = patch.parent
        extfile = patch.extentfile+'.extent.wkt'
        buffile = patch.extentfile+'.buffer.wkt'

        # Build a grid tool command to process the file
        
        # Read in the existing file (only need de, dn, du components)
        merge_commands='gridtool read maxcols 3 gridfile '

        # If have a parent, then subtract it as want to merge into it
        if parent: merge_commands = merge_commands + ' subtract csv parentcsv'

        # Set to zero outside area extents for which patch subcell resolution is required
        merge_commands=merge_commands + ' zero outside extents'

        # Smooth out to buffer to avoid sharp transition
        merge_commands=merge_commands + ' smooth linear outside extents not outside buffer'

        # If we have a parent and are not making additive patches then add the parent back
        # on.
        if parent and not additive: merge_commands = merge_commands + ' add csv parentcsv'

        # Trim redundant zeros around the edge of the patch
        merge_commands=merge_commands + ' trim 1'

        # Finally write out as a CSV file
        merge_commands=merge_commands + ' write csv csvfile'

        # Run the commands
        gridtoolcmd=[ 
            gridfile if x == 'gridfile' else
            csvfile if x == 'csvfile' else
            parent.csv if x == 'parentcsv' else
            extfile if x == 'extents' else
            buffile if x == 'buffer' else
            x
            for x in merge_commands.split()
            ]

        call(gridtoolcmd)

        if not os.path.exists( csvfile ):
            raise RuntimeError("Could not create CSV file "+csvfile)
        patch.csv=csvfile


def write_log( message, level=0 ):
    prefix='  '*level
    if logfile:
        logfile.write(prefix)
        logfile.write(message)
        logfile.write("\n")
    if verbose:
        print prefix+message

def write_wkt( item, level, wkt ):
    wktfile.write("{0}|{1}|{2}\n".format(item,level,wkt))

if __name__ == "__main__":

    # Process arguments

    parser=argparse.ArgumentParser(description='Build set of grids for deformation patch')
    parser.add_argument('patch_file',help='Base name used for output files')
    parser.add_argument('model_file',help='Model file(s) used to calculate deformation, passed to calc_okada',nargs='*')

    args=parser.parse_args()
    patchfile = args.patch_file
    patchpath,patchname=os.path.split(args.patch_file)
                
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

    logfile=open(patchfile+".build.log","w")
    wktfile=open(patchfile+".build.wkt","w")
    logfile.write("Building deformation grids for model: {0}\n".format(modeldef))
    wktfile.write("description|level|wkt\n")

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
        bounds = expanded_bounds( buffered_extent )
        write_wkt("Expanded extents {0}".format(name),0,bounds_wkt(bounds))
        griddef = bounds_grid_def( bounds, base_size )
        extfile=os.path.join(patchpath,name)
        with open(extfile+".extent.wkt","w") as f:
            f.write(extent.to_wkt())
            f.write("\n")
        with open(extfile+".buffer.wkt","w") as f:
            f.write(buffered_extent.to_wkt())
            f.write("\n")
        create_grids(gridlist, modeldef,patchpath,name,1,griddef,base_size,extentfile=extfile)

    create_patch_csv( gridlist )
    print gridlist




