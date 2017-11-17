#!/usr/bin/python
# Script to build a patch definition file from an NZGD2000 deformation model.
# 
# This creates representation of the vertical deformation only as a set of bin (geoid format) files
# The files are ordered from the finest to the coarsest grid so that the finest grids are used in 
# preference to coarsergrids.
#
# Based on FME multiple grid capability using GDC file described in
#  https://knowledge.safe.com/articles/29325/creating-vertical-adjustment-grid-files-for-use-wi.html
# Binary grid file defined in 
#  https://www.nrcan.gc.ca/sites/www.nrcan.gc.ca/files/earthsciences/pdf/gpshgrid_e.pdf

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
from ellipsoid import grs80

from LINZ.DeformationModel.Model import Model

gridtool=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'gridtool','gridtool')
if not os.path.exists(gridtool):
    raise RuntimeError('Cannot find gridtool program at '+gridtool)

class GridExtents( object ):

    defaultTolerance=0.0000001

    @staticmethod
    def fromFile( gridfile ):
        with open(gridfile) as gf:
            header=gf.readline()
            xmin=None
            ymin=None
            xmax=None
            ymax=None
            xinc=None
            yinc=None
            for l in gf:
                l=l.strip()
                if l:
                    x,y=(float(x) for x in l.split(',')[:2])
                    if xmin is None:
                        xmin=x
                        xmax=x
                        ymin=y
                        ymax=y
                        x0=x
                        y0=y
                    else:
                        if x < xmin: xmin=x
                        elif x > xmax: xmax=x
                        if y < ymin: ymin=y
                        elif y > ymax: ymax=y
                        if xinc is None and x != x0: xinc=x-x0
                        if yinc is None and y != y0: yinc=y-y0
        if xinc is None or yinc is None:
            return None
        cellsize=(round(xinc,8),round(yinc,8))
        return GridExtents(0,xmin,ymin,xmax,ymax,cellsize,gridfile)

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

def linearGridList( levels ):
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

def runCommand( command ):
    subprocess.check_output(command)

def main():
    import argparse

    parser=argparse.ArgumentParser("Build gdc (NGA format geoid file list) vertical reverse patch definition file")
    parser.add_argument('model_dir',help="Model directory")
    parser.add_argument('build_dir',help="Patch build directory")
    parser.add_argument('patch_name',nargs="?",help="Patch file name")
    parser.add_argument('--version',help="Deformation model for which patch applies, default is current version")
    parser.add_argument('-B','--big-endian',action='store_true',help='File is big endian (default little endian)')
    
    args=parser.parse_args()
    modeldir=args.model_dir
    builddir=args.build_dir

    if not os.path.isdir(builddir):
        raise RuntimeError('Build directory {0} does not exist or is not a directory'
                           .format(builddir))

    model=Model(modeldir)
    version=args.version
    if version is None:
        version=model.version()
    datumcode=model.metadata('datum_code')

    patchname=args.patch_name
    patchdir=''
    if patchname is not None:
        if os.path.isdir(patchname):
            patchdir=patchname
            patchname=None
        else:
            patchdir=os.path.dirname(patchname)
            patchname=os.path.basename(patchname)
    if not patchname:
        patchname="{0}_{1}_V".format(datumcode,version)
    if not patchdir:
        patchdir='.'

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
        if spatial.displacement_type not in ('3d','vertical'):
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
                fgfile=os.path.join(builddir,'gcd_tmp_scaled_{0}.csv'.format(nscaled))
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
        # outside extents of grid and can be added to build the total deformation.
        # This allows separate levels to be easily combined, eg co + post seismic.

        pfiles={}
        for level,extents in enumerate(levels):
            for e in extents:
                cfile=os.path.join(builddir,'gdc_tmp_{0}_L{1}_{2}.csv'.format(e.compid,level,e.id))
                rfile=cfile
                pfile=pfiles.get(e.parent)
                spec=[str(x) for x in e.gridspec()]
                command=[gridtool,'create','gridspec']
                command.extend(spec)
                command.extend(('columns','du'))
                if pfile:
                    command.extend(('add','csv',pfile))
                for gf in e.gridfiles:
                    command.extend(('replace','csv',gf))
                command.extend(('write','csv',cfile))
                if pfile:
                    rfile=os.path.join(builddir,'gdc_tmp_{0}_L{1}_{2}_r.csv'.format(e.compid,level,e.id))
                    command.extend((
                        'subtract','csv',pfile,
                        'trim','tolerance','0.000001','noexpand','1',
                        'write','csv',rfile))
                runCommand(command)
                if not os.path.exists(rfile):
                    raise RuntimeError("Failed to create compiled grid {0}".format(cfile))
                pfiles[e]=cfile
                e.compiled_gridfile=cfile
                cellsize=e.cellsize
                if cellsize not in levelgrids:
                    levelgrids[cellsize]=[]
                compgrid=GridExtents.fromFile(rfile)
                levelgrids[cellsize].append(compgrid)
                # Remove temporary files
                for gf in e.gridfiles:
                    if gf not in sourcegrids:
                        os.unlink(gf)
            
    # Now compile by adding component grids for each level


    levels=[levelgrids[c] for c in reversed(sorted(levelgrids))]
    compileExtents(levels)
    compgrids={}
    gridmap={}
    for level,extents in enumerate(levels):
        for e in extents:
            cfile=os.path.join(builddir,'gdc_tmp_L{1}_{2}.csv'
                               .format(e.compid,level,e.id))
            pfile=None
            pgrid=gridmap.get(e.parent)
            if pgrid is not None:
                pfile=pgrid.gridfiles[0]
            spec=[str(x) for x in e.gridspec()]
            command=[gridtool,'create','gridspec']
            command.extend(spec)
            command.extend(('columns','du'))
            for gf in e.gridfiles:
                command.extend(('add','csv',gf))
            if pfile:
                command.extend((
                    'trim','tolerance','0.000001','noexpand','1',
                    'alignto','csv',pfile,
                    'add','csv',pfile))
            command.extend(('write','csv',cfile))
            runCommand(command)
            if not os.path.exists(cfile):
                raise RuntimeError("Failed to create compiled grid {0}".format(cfile))
            cellsize=e.cellsize
            if cellsize not in compgrids:
                compgrids[cellsize]=[]
            compgrid=GridExtents.fromFile(cfile)
            compgrid.parent=pgrid
            gridmap[e]=compgrid
            compgrids[cellsize].append(compgrid)

    levels=[compgrids[c] for c in reversed(sorted(compgrids))]

    # Could build in option for splitting and trimming here - would reduce
    # grid size significantly for Kaikoura deformation with strong SW-NE trend

    # Now generate geoid grids
    gridlist=linearGridList(levels)
    gridlist.reverse()

    endian = '>' if args.big_endian else '<'
    iformat=endian+'i'
    dformat=endian+'d'
    hformat=endian+'f'

    
    with open(os.path.join(patchdir,patchname+'.gdc'),'w') as gdch:
        for id,grid in enumerate(gridlist):
            gridname="g{0}z{1:02d}.bin".format(version[:4],id+1)
            gdch.write("{0}\n".format(gridname))
            with open(os.path.join(patchdir,gridname),'wb') as nb:
                lon0,lat0,lon1,lat1,nlon,nlat=grid.gridspec()
                loninc=(lon1-lon0)/nlon
                latinc=(lat1-lat0)/nlat
                nb.write(struct.pack(dformat,lat0))
                nb.write(struct.pack(dformat,lon0))
                nb.write(struct.pack(dformat,latinc))
                nb.write(struct.pack(dformat,loninc))
                nb.write(struct.pack(iformat,nlat+1))
                nb.write(struct.pack(iformat,nlon+1))
                nb.write(struct.pack(iformat,1))

                gridfiles=grid.gridfiles
                if len(gridfiles) != 1:
                    raise RuntimeError('Invalid grid extents for grid file - grids not compiled')
                gridfile=gridfiles[0]
            
                npt=(nlon+1)*(nlat+1)
                nread=0
                with open(gridfile) as lgf:
                    lgf.readline()
                    for line in lgf:
                        if ',' not in line:
                            break
                        lon,lat,du=[float(f) for f in line.split(',')][0:4]
                        nb.write(struct.pack(hformat,du))
                        nread += 1

                if nread != npt:
                    raise RuntimeError('Number of points {3} in grid {0} does not match grid spec ({1},{2})'
                            .format(gridfile,nlon+1,nlat+1,len(points)))
                            
if __name__=='__main__':
    main()
