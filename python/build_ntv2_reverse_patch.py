#!/usr/bin/python
# Script to build a patch definition file from an NZGD2000 deformation model.


import argparse
import math
import os
import os.path
import random
import re
import struct
import subprocess
import sys
from datetime import datetime,date
from shapely import wkt, affinity
from shapely.geometry import MultiPoint

from LINZ.DeformationModel.Model import Model
from LINZ.Geodetic.Ellipsoid import GRS80

gridtool=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'gridtool','gridtool')
if not os.path.exists(gridtool):
    raise RuntimeError('Cannot find gridtool program at '+gridtool)

class GridExtents( object ):

    defaultTolerance=0.0000001

    def __init__( self, compid, xmin, ymin, xmax, ymax, cellsize, gridfile ):
        self.compid=compid
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax
        self.cellsize=cellsize
        self.gridfiles=list(gridfile) if isinstance(gridfile,list) else [gridfile]
        self.id=None
        self.compiled_gridfile=None
        self.parent=None
        self.children=[]
        
    def overlaps( self, other, tolerance=defaultTolerance ):
        return not (self.xmin-tolerance > other.xmax or
                other.xmin-tolerance > self.xmax or
                self.ymin-tolerance > other.ymax or
                other.ymin-tolerance > self.ymax)

    def contains( self, other, tolerance=defaultTolerance ):
        return (self.xmin-tolerance <= other.xmin and
                self.xmax-tolerance >= other.xmax and
                self.ymin-tolerance <= other.ymin and
                other.ymax+tolerance >= other.ymax)
    
    def mergedWith( self, other ):
        if self.cellsize != other.cellsize:
            raise RuntimeError("Cannot merge grids with different cellsize")
        gridfiles=list(self.gridfiles)
        gridfiles.extend(other.gridfiles)
        return GridExtents( 
            self.compid,
            min(self.xmin,other.xmin),
            min(self.ymin,other.ymin),
            max(self.xmax,other.xmax),
            max(self.ymax,other.ymax),
            self.cellsize,
            gridfiles)

    def compiledGrid( self ):
        if self.compiled_gridfile is None:
            raise RuntimeError('Compiled grid file not defined..')
        return GridExtents(0,self.xmin,self.ymin,self.xmax,self.ymax,
                           self.cellsize, self.compiled_gridfile)

    def gridspec( self ):
        ngx=int(math.floor((self.xmax-self.xmin)/self.cellsize[0]+0.01))
        ngy=int(math.floor((self.ymax-self.ymin)/self.cellsize[1]+0.01))
        return ( 
            self.xmin, self.ymin,
            self.xmax, self.ymax,
            ngx, ngy
            )

    def __str__( self ):
        return ":".join((str(x) for x in self.gridspec()))

def compileExtents( levels ):
    # Levels is an array of arrays of GridExtent objects.  For each level
    # merge overlapping extents, then form parent child relationships between
    # merged grids at each level
    #
    # Combine adjoining grids or overlapping grid extents for each level
    for level,extents in enumerate(levels):
        updated=True
        while updated:
            newextents=[]
            updated=False
            for g0 in extents:
                merged=False
                for i,g1 in enumerate(newextents):
                    if g0.overlaps(g1):
                        newextents[i]=g0.mergedWith(g1)
                        updated=True
                        merged=True
                        break
                if not merged:
                    newextents.append(g0)
            extents=newextents
        levels[level]=newextents
    
    # Identify parent/child grids for each level and assign ids to the grids.
    nextid=0
    for level,extents in enumerate(levels):
        for e in extents:
            nextid += 1
            e.id=nextid
            for pextents in reversed(levels[:level]):
                for p in pextents:
                    if p.contains(e):
                        e.parent=p
                        p.children.append(e)
                        break
                if e.parent is not None:
                    break
    return levels

def printLevels( levels ):
    ''' Print list of grids for debugging purposes '''
    for extents in levels:
        print '\nCellsize: {0}'.format(extents[0].cellsize)
        for e  in extents:
            print "   {0} {1} {2}".format(e.id,e.compid,e.gridspec())
            if e.parent is not None:
                print "      Parent {0}".format(e.parent.id)
            for gf in e.gridfiles:
                fn=os.path.basename(gf)
                print "      {0}".format(fn)

def ntv2GridList( levels ):
    gridlist=[]
    def _addGrid( grid ):
        if grid in gridlist:
            return
        gridlist.append(grid)
        for c in grid.children:
            _addGrid(c)
    for level in levels:
        for grid in level:
            _addGrid(grid)
    return gridlist

calc_lat=None
calc_dedln=None
calc_dndlt=None

def calcDlonDlat( lon, lat, de, dn ):
    if lat != calc_lat:
        calc_dedln,calc_dndlt=GRS80.metres_per_degree(lat)
    return de/dedln, dn/dndlt

def runCommand( command ):
    subprocess.check_output(command)

def main():
    import argparse

    parser=argparse.ArgumentParser("Build NTv2 reverse patch definition file")
    parser.add_argument('model_dir',help="Model directory")
    parser.add_argument('build_dir',help="Patch build directory")
    parser.add_argument('patch_name',nargs="?",help="Patch file name")
    parser.add_argument('-p','--test-points',help="Test point file")
    parser.add_argument('-n','--n-test-points',type=int,default=10,help="Number of test points per grid")
    parser.add_argument('--version',help="Deformation model for which patch applies, default is current version")
    parser.add_argument('-c','--created',default=date.today().strftime("%Y%m%d"),help='Created date')
    parser.add_argument('-A','--australian',action='store_true',help='Build australian format binary')
    parser.add_argument('-B','--big-endian',action='store_true',help='File is big endian (default little endian)')
    parser.add_argument('--ogr2ogr-bug-workaround',action='store_true',help='Workaround for ogr2ogr nested grid bug (#177) - carefully ordered non-nested grids!')
    
    args=parser.parse_args()
    modeldir=args.model_dir
    builddir=args.build_dir
    test_point_file=args.test_points
    n_test_points=args.n_test_points

    if not os.path.isdir(builddir):
        raise RuntimeError('Build directory {0} does not exist or is not a directory'
                           .format(builddir))

    model=Model(modeldir)
    version=args.version
    if version is None:
        version=model.version()
    datumcode=model.metadata('datum_code')
    modelname=model.metadata('model_name')


    patchname=args.patch_name
    if patchname is None:
        patchname='{0}_patch_{1}'.format(datumcode,version)
    ntv2created=args.created

    ngrid=0
    ncomp=0
    levelgrids={}
    sourcegrids=[]

    # Find all grids contributing to the reverse patch
    revcomps=model.reversePatchComponents(version)
    nscaled=0

    for factor,c in revcomps:
        ncomp += 1
        spatial=c.spatialModel
        if spatial.displacement_type not in ('3d','horizontal'):
            continue
        compgrids={}
        for gm in spatial.models():
            grid=gm.model()
            if type(grid).__name__ != 'Grid':
                raise RuntimeError('Only grid models handled by build_patch.py')
            minlon,minlat,maxlon,maxlat,nlon,nlat=grid.gridSpec()
            cellsize=(round((maxlon-minlon)/(nlon-1),8),round((maxlat-minlat)/(nlat-1),8))
            if cellsize not in compgrids:
                compgrids[cellsize]=[]
            gridfile=model.getFileName(grid.gridFile())
            sourcegrids.append(gridfile)
            if factor != 1.0:
                nscaled+=1
                fgfile=os.path.join(builddir,'ntv2_tmp_scaled_{0}.csv'.format(nscaled))
                command=[gridtool,
                          'read','csv','maxcols','2',gridfile,
                          'multiply',str(factor),
                          'write','csv',fgfile]
                runCommand(command)
                if not os.path.exists(fgfile):
                    raise RuntimeError('Cannot create scaled component file {0}'.format(fgfile))
                gridfile=fgfile
            extents=GridExtents(ncomp,minlon,minlat,maxlon,maxlat,cellsize,gridfile)
            compgrids[cellsize].append(extents)

        levels=[compgrids[c] for c in reversed(sorted(compgrids))]
        compileExtents(levels)

        # Compile component gridfiles for each level, trimming zeros and aligning
        # to parent grid. Final grids have parent subtracted, so are zero 
        # outside extents of grid. This will allow them to be added together 
        # at the next step.

        for level,extents in enumerate(levels):
            for e in extents:
                cfile=os.path.join(builddir,'ntv2_tmp_{0}_L{1}_{2}.csv'.format(e.compid,level,e.id))
                pfile=e.parent.compiled_gridfile if e.parent is not None else None
                spec=[str(x) for x in e.gridspec()]
                command=[gridtool,'create','gridspec']
                command.extend(spec)
                command.extend(('columns','de+dn'))
                if pfile:
                    command.extend(('add','csv',pfile))
                for gf in e.gridfiles:
                    command.extend(('replace','csv',gf))
                if pfile:
                    command.extend((
                        'subtract','csv',pfile,
                        'trim','tolerance','0.000001','noexpand','1'))
                command.extend(('write','csv',cfile))
                runCommand(command)
                if not os.path.exists(cfile):
                    raise RuntimeError("Failed to create compiled grid {0}".format(cfile))
                e.compiled_gridfile=cfile
                cellsize=e.cellsize
                if cellsize not in levelgrids:
                    levelgrids[cellsize]=[]
                levelgrids[cellsize].append(e.compiledGrid())
                # Remove temporary files
                for gf in e.gridfiles:
                    if gf not in sourcegrids:
                        os.unlink(gf)
            
    # Now compile by adding component grids for each level

    levels=[levelgrids[c] for c in reversed(sorted(levelgrids))]
    compileExtents(levels)
    compgrids={}
    for level,extents in enumerate(levels):
        for e in extents:
            cfile=os.path.join(builddir,'ntv2_tmp_L{1}_{2}.csv'
                               .format(e.compid,level,e.id))
            pfile=e.parent.compiled_gridfile if e.parent is not None else None
            spec=[str(x) for x in e.gridspec()]
            command=[gridtool,'create','gridspec']
            command.extend(spec)
            command.extend(('columns','de+dn'))
            for gf in e.gridfiles:
                command.extend(('add','csv',gf))
            if pfile:
                command.extend((
                    'subtract','csv',pfile,
                    'trim','tolerance','0.000001','noexpand','1',
                    'add','csv',pfile))
            command.extend(('write','csv',cfile))
            runCommand(command)
            if not os.path.exists(cfile):
                raise RuntimeError("Failed to create compiled grid {0}".format(cfile))
            e.compiled_gridfile=cfile
            cellsize=e.cellsize
            if cellsize not in compgrids:
                compgrids[cellsize]=[]
            compgrids[cellsize].append(e.compiledGrid())

    levels=[compgrids[c] for c in reversed(sorted(compgrids))]

    # Create test point files for each grid

    if test_point_file:
        with open(test_point_file,'w') as tpf:
            tpf.write('lon,lat\n')
            for g in ntv2GridList(levels):
                for i in range(n_test_points):
                    lon=random.uniform(g.xmin,g.xmax)
                    lat=random.uniform(g.ymin,g.ymax)
                    tpf.write("{0:.5f},{1:.5f}\n".format(lon,lat))

    # Could build in option for splitting and trimming here - would reduce
    # grid size significantly for Kaikoura deformation with strong SW-NE trend


    # Now generate NTv2 grids
    gridlist=ntv2GridList(levels)
    for i,g in enumerate(gridlist):
        g.id=i+1

    if args.ogr2ogr_bug_workaround:
        print "Applying ogr2ogr bug (proj issue #177) bug fix"
        gridlist.reverse()

    endian = '>' if args.big_endian else '<'
    sformat=endian+'8s8s'
    iformat=endian+'8si' if args.australian else endian+'8si4x'
    dformat=endian+'8sd'
    gsformat=endian+'ffff'

    ntv2created=(ntv2created+'       ')[:8]
    ntv2version=(version+'       ')[:8]

    with open(patchname+'.asc','w') as nt,open(patchname+'.gsb','wb') as nb:
        nt.write("{0:8s}{1:3d}\n".format('NUM_OREC',11))
        nt.write("{0:8s}{1:3d}\n".format('NUM_SREC',11))
        nt.write("{0:8s}{1:3d}\n".format('NUM_FILE',len(gridlist)))
        nt.write("{0:8s}{1:8s}\n".format('GS_TYPE','SECONDS'))
        nt.write("{0:8s}{1:8s}\n".format('VERSION',ntv2version))
        nt.write("{0:8s}{1:8s}\n".format('SYSTEM_F','GRS80'))
        nt.write("{0:8s}{1:8s}\n".format('SYSTEM_T','GRS80'))
        nt.write("{0:8s}{1:12.3f}\n".format('MAJOR_F',GRS80.a))
        nt.write("{0:8s}{1:12.3f}\n".format('MINOR_F',GRS80.b))
        nt.write("{0:8s}{1:12.3f}\n".format('MAJOR_T',GRS80.a))
        nt.write("{0:8s}{1:12.3f}\n".format('MINOR_T',GRS80.b))

        nb.write(struct.pack(iformat,'NUM_OREC',11))
        nb.write(struct.pack(iformat,'NUM_SREC',11))
        nb.write(struct.pack(iformat,'NUM_FILE',len(gridlist)))
        nb.write(struct.pack(sformat,'GS_TYPE ','SECONDS '))
        nb.write(struct.pack(sformat,'VERSION ',ntv2version))
        nb.write(struct.pack(sformat,'SYSTEM_F','GRS80   '))
        nb.write(struct.pack(sformat,'SYSTEM_T','GRS80   '))
        nb.write(struct.pack(dformat,'MAJOR_F ',GRS80.a))
        nb.write(struct.pack(dformat,'MINOR_F ',GRS80.b))
        nb.write(struct.pack(dformat,'MAJOR_T ',GRS80.a))
        nb.write(struct.pack(dformat,'MINOR_T ',GRS80.b))

        ntgridname={}
        for grid in gridlist:
            name="GRID{0:02d}  ".format(grid.id)
            ntgridname[grid]=name
            pname=ntgridname.get(grid.parent,'NONE    ')
            lon0,lat0,lon1,lat1,nlon,nlat=grid.gridspec()
            lon0 *= 3600.0
            lat0 *= 3600.0
            lon1 *= 3600.0
            lat1 *= 3600.0
            loninc=(lon1-lon0)/nlon
            latinc=(lon1-lon0)/nlat
            npt=(nlon+1)*(nlat+1)

            gridfiles=grid.gridfiles
            if len(gridfiles) != 1:
                raise RuntimeError('Invalid grid extents for NTv2 file - grids not compiled')
            gridfile=gridfiles[0]

            nt.write("{0:8s}{1:8s}\n".format('SUB_NAME',name))
            nt.write("{0:8s}{1:8s}\n".format('PARENT',pname))
            nt.write("{0:8s}{1:8s}\n".format('CREATED',ntv2created))
            nt.write("{0:8s}{1:8s}\n".format('UPDATED',ntv2created))
            nt.write("{0:8s}{1:15.6f}\n".format('S_LAT',lat0))
            nt.write("{0:8s}{1:15.6f}\n".format('N_LAT',lat1))
            nt.write("{0:8s}{1:15.6f}\n".format('E_LONG',-lon1))
            nt.write("{0:8s}{1:15.6f}\n".format('W_LONG',-lon0))
            nt.write("{0:8s}{1:15.6f}\n".format('LAT_INC',latinc))
            nt.write("{0:8s}{1:15.6f}\n".format('LONG_INC',loninc))
            nt.write("{0:8s}{1:6d}\n".format('GS_COUNT',npt))

            nb.write(struct.pack(sformat,'SUB_NAME',name))
            nb.write(struct.pack(sformat,'PARENT  ',pname))
            nb.write(struct.pack(sformat,'CREATED ',ntv2created))
            nb.write(struct.pack(sformat,'UPDATED ',ntv2created))
            nb.write(struct.pack(dformat,'S_LAT   ',lat0))
            nb.write(struct.pack(dformat,'N_LAT   ',lat1))
            nb.write(struct.pack(dformat,'E_LONG  ',-lon1))
            nb.write(struct.pack(dformat,'W_LONG  ',-lon0))
            nb.write(struct.pack(dformat,'LAT_INC ',latinc))
            nb.write(struct.pack(dformat,'LONG_INC',loninc))
            nb.write(struct.pack(iformat,'GS_COUNT',npt))
        
            points=[]
            with open(gridfile) as lgf:
                lgf.readline()
                for line in lgf:
                    if ',' not in line:
                        break
                    dln,dlt=[float(f)*3600 for f in line.split(',')][2:4]
                    points.append([dln,dlt])

            if len(points) != npt:
                raise RuntimeError('Number of points {3} in grid {0} does not match grid spec ({1},{2})'
                        .format(gridfile,nlon+1,nlat+1,len(points)))


            for ilat in range(nlat+1):
                ilatn=ilat*nlon
                for ilon in reversed(range(nlon+1)):
                    dln,dlt = points[ilatn+ilon]
                    nt.write("{0:10.6f}{1:10.6f}{2:10.6f}{3:10.6f}\n".format(dlt,-dln,-1.0,-1.0))
                    nb.write(struct.pack(gsformat,dlt,-dln,-1.0,-1.0))
                        
        nt.write("{0:8s}\n".format('END'))
        nb.write(struct.pack(sformat,'END     ','\0\0\0\0\0\0\0\0'))
                
if __name__=='__main__':
    main()
