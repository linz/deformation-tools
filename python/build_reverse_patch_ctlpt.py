#!/usr/bin/python
# Script to build a set of control points defining a reverse patch.  The control points are defined 
# by the grid points at each level.  Grid points are discarded if all the grid points in a parent
# cell have values insignificantly different from their parent.

import argparse
import math
import os
import os.path
import random
import re
import struct
import subprocess
import sys
import numpy as np
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
        print('\nCellsize: {0}'.format(extents[0].cellsize))
        for e  in extents:
            print("   {0} {1} {2}".format(e.id,e.compid,e.gridspec()))
            if e.parent is not None:
                print("      Parent {0}".format(e.parent.id))
            for gf in e.gridfiles:
                fn=os.path.basename(gf)
                print("      {0}".format(fn))

def ctlptGridList( levels ):
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
        calc_dedln,calc_dndlt=grs80.metres_per_degree(lon,lat)
    return de/calc_dedln, dn/calc_dndlt

def runCommand( command, verbose=False ):
    try:
        result=subprocess.check_output(command)
        if verbose:
            print(result)
    except Exception as ex:
        print("Command parameters:")
        for c in command:
            print("  {0}".format(c))
        raise ex

used_ordinates={
    ('3d','3d'): ('de','dn','du'),
    ('3d','horizontal'): ('de','dn'),
    ('3d','vertical'): ('du',),
    ('horizontal','3d'): ('de','dn'),
    ('horizontal','horizontal'): ('de','dn'),
    ('vertical','3d'): ('du',),
    ('vertical','vertical'): ('du',),
    }

output_ordinates={
    '3d': ('de','dn','du'),
    'horizontal': ('de','dn'),
    'vertical': ('du',),
    }



def main():
    import argparse

    parser=argparse.ArgumentParser("Build NTv2 reverse patch definition file")
    parser.add_argument('model_dir',help="Model directory")
    parser.add_argument('build_dir',help="Patch build directory")
    parser.add_argument('control_point_name',nargs="?",help="Control point file name")
    parser.add_argument('--ordinates',choices=('3d','horizontal','vertical'),default='3d',help="Ordinates required")
    parser.add_argument('--version',help="Deformation model for which patch applies, default is current version")
    parser.add_argument('-v','--verbose',action='store_true',help="More output")
    
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
    modelname=model.metadata('model_name')


    ordval='' if args.ordinates != '3d' else '_'+args.ordinates[0].upper()
    ctlpt_file=args.control_point_name
    if ctlpt_file is None:
        ctlpt_file='{0}_{1}{2}.csv'.format(datumcode,version,ordval)

    if os.path.exists(ctlpt_file):
        os.remove(ctlpt_file)

    ctlpt_ordinates=output_ordinates[args.ordinates]

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
        ordinates=used_ordinates.get((args.ordinates,spatial.displacement_type))
        if ordinates is None:
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
                fgfile=os.path.join(builddir,'ctlpt_tmp_scaled_{0}{1}.csv'.format(nscaled,ordval))
                command=[gridtool,
                          'read','csv',gridfile,
                          'multiply',str(factor),
                          'write','csv',fgfile]
                runCommand(command,args.verbose)
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
                cfile=os.path.join(builddir,'ctlpt_tmp_{0}_L{1}{2}_{3}.csv'
                                   .format(e.compid,level,ordval,e.id))
                rfile=cfile
                pfile=pfiles.get(e.parent)
                spec=[str(x) for x in e.gridspec()]
                command=[gridtool,'create','gridspec']
                command.extend(spec)
                command.extend(('columns','+'.join(ctlpt_ordinates)))
                if pfile:
                    command.extend(('add','csv',pfile))
                for gf in e.gridfiles:
                    command.extend(('replace','matching_from','csv',gf))
                command.extend(('write','csv',cfile))
                if pfile:
                    rfile=os.path.join(builddir,'ctlpt_tmp_{0}_L{1}{2}_{3}_r.csv'
                                       .format(e.compid,level,ordval,e.id))
                    command.extend((
                        'subtract','csv',pfile,
                        'trim','tolerance','0.000001','noexpand','1',
                        'write','csv',rfile))
                runCommand(command,args.verbose)
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
            cfile=os.path.join(builddir,'ntv2_tmp_L{1}{2}_{3}.csv'
                               .format(e.compid,level,ordval,e.id))
            pfile=None
            pgrid=gridmap.get(e.parent)
            if pgrid is not None:
                pfile=pgrid.gridfiles[0]
            spec=[str(x) for x in e.gridspec()]
            command=[gridtool,'create','gridspec']
            command.extend(spec)
            command.extend(('columns','+'.join(ctlpt_ordinates)))
            for gf in e.gridfiles:
                command.extend(('add','csv',gf))
            if pfile:
                command.extend((
                    'trim','tolerance','0.000001','noexpand','1',
                    'alignto','csv',pfile,
                    'add','csv',pfile))
            command.extend(('write','csv',cfile))
            runCommand(command,args.verbose)
            if not os.path.exists(cfile):
                raise RuntimeError("Failed to create compiled grid {0}".format(cfile))
            cellsize=e.cellsize
            if cellsize not in compgrids:
                compgrids[cellsize]=[]
            compgrid=GridExtents.fromFile(cfile)
            compgrid.parent=pgrid
            gridmap[e]=compgrid
            compgrids[cellsize].append(compgrid)

    # Build up final control point file, starting with the finest grid

    testcols='|'+'|'.join(ctlpt_ordinates)+'|'
    append='write'
    for cellsize in sorted(compgrids):
        for e in compgrids[cellsize]:
            efile=e.gridfiles[0]
            pfile=e.parent.gridfiles[0] if e.parent else None
            if True:
                efname=os.path.basename(efile)
                pfname=os.path.basename(pfile) if pfile else ''
                print("Adding {0}: {1}".format(efname,pfname))
            command=[gridtool,'read','csv',efile]
            if pfile:
                command.extend(('subtract','csv',pfile))
            command.extend(('mark','where',testcols,'>','0.0001','plus','expand','1.5'))
            if append == 'append':
                command.extend(('minus','nearest_to',ctlpt_file,'within','0.00000001'))
            if pfile:
                command.extend(('add','csv',pfile))
            command.extend((append,'dos','csv',ctlpt_file,'where','marked'))
            runCommand(command,args.verbose)
            append='append'
                
if __name__=='__main__':
    main()
