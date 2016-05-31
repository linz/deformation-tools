#!/usr/bin/python
#
# Script to compile zip files of the reverse patch deformation grids applied to Landonline
# in a format suitable for distributing to 3rd parties.
#
# Also compiles a list of test points for testing the grid
#

import sys
import os
import os.path
import re
import numpy as np
import ellipsoid
from defgrid import DeformationGrid
from random import uniform, seed
from subprocess import call

ntestpergrid=20

land_area_tolerance=0.01
land_areas=None
seed('Consistent random seed')

ntv2py=os.path.join(os.path.dirname(os.path.realpath(sys.argv[0])),'make_ntv2.py')

def load_land_areas( polygon_file ):
    global land_areas
    if polygon_file:
        from shapely import wkt
        from shapely.geometry import shape, MultiPolygon, Polygon

        print "Loading land area definition from",polygon_file
        with open(polygon_file) as f:
            pgnwkt=f.read()
            land_areas=wkt.loads(pgnwkt)

def build_dump( def_file ):
    components=[]
    component=None
    basedir=os.path.dirname(def_file)
    testpts=[]
    with open(def_file) as f:
        for l in f:
            m=re.match(r'\s*SHIFT_MODEL_COMPONENT\s+(\w+\.\w+)\s*$',l)
            if m:
                gfile=m.group(1)
                if basedir:
                    gfile=os.path.join(basedir,gfile)
                    if not os.path.exists(gfile):
                        raise RuntimeError('Cannot find component '+gfile)
                component={'file':gfile,'negative':False}
                components.append(component)
                continue
            if not component:
                continue
            m=re.match(r'\s*MODEL_TYPE\s+(\w+)\s*$',l.upper())
            if m:
                if m.group(1) != 'GRID':
                    raise RuntimeError('Can only handle GRID type components')
                component['type']=m.group(1)
                continue
            m=re.match(r'\s*NEGATIVE\s+(\w+)\s*$',l.upper())
            if m:
                component['negative']=m.group(1)=='YES'
                continue

    csvref=re.sub(r'(\.def)?$','.csv',def_file,re.I)
    filelist=[csvref]
    gridfiles=[]
    with open(csvref,'w') as sf:
        sf.write('gridfile,minlon,maxlon,nlon,minlat,maxlat,nlat,wkt\r\n')
        for c in components:
            gfile=c['file']
            cfile=re.sub(r'(\.gdf)?$','.csv',gfile,re.I)
            cname=os.path.basename(cfile)
            c['csvfile']=cname
            factor = -1.0 if c['negative'] else 1.0
            xmin=0.0
            xmax=0.0
            ymin=0.0
            ymax=0.0
            ngrdx=0
            ngrdy=0
            dx=0.0
            dy=0.0
            headers=False
            started=False
            npt=0
            with open(gfile) as gf:
                with open(cfile,'w') as cf:
                    for l in gf:
                        if not started:
                            m=re.match(r'^\s*(\w+)\:\s+(.*?)\s*$',l)
                            if m:
                                record=m.group(1)
                                value=m.group(2)
                                if record == 'NGRDX':
                                    ngrdx=int(value)
                                elif record == 'NGRDY':
                                    ngrdy=int(value)
                                elif record == 'XMIN':
                                    xmin=float(value)
                                elif record == 'XMAX':
                                    xmax=float(value)
                                elif record == 'YMIN':
                                    ymin=float(value)
                                elif record == 'YMAX':
                                    ymax=float(value)
                                elif record == 'NDIM':
                                    if value=='2':
                                        cf.write('lon,lat,de,dn\r\n')
                                    elif value=='3':
                                        cf.write('lon,lat,de,dn,du\r\n')
                                    else:
                                        raise RuntimeError('Cannot handle dimensions other than 2 and 3')
                                    headers=True
                                elif record == 'VALUES':
                                    if value != 'REAL':
                                        raise RuntimeError('Can only use "REAL" data values')
                                    continue
                        m=re.match(r'^\s*V(\d+)\,(\d+)\:\s*(.*?)\s*$',l)
                        if m:
                            if not started:
                                dx=(xmax-xmin)/(ngrdx-1)
                                dy=(ymax-ymin)/(ngrdy-1)
                                if not headers:
                                    raise RuntimeError('Oops, didn\'t find NDIM record')
                                started=True
                            nx=int(m.group(1))-1
                            ny=int(m.group(2))-1
                            values=[xmin+nx*dx,ymin+ny*dy]
                            values.extend((float(p)*factor for p in m.group(3).split()))
                            cf.write(','.join(str(v) for v in values))
                            cf.write('\r\n')
                            npt+=1
            if not started:
                raise RuntimeError('No values read from gdf file')
            if npt != ngrdx*ngrdy:
                raise RuntimeError('Number of points {0} does\'t match expected {1}'.format(
                    npt,ngrdx*ngrdy))
            wkt='"POLYGON(({0} {2},{1} {2},{1} {3},{0} {3},{0} {2}))"'.format(xmin,xmax,ymin,ymax)
            sf.write(','.join(str(v) for v in (cname,xmin,xmax,ngrdx,ymin,ymax,ngrdy,wkt)))
            sf.write('\r\n')
            filelist.append(cfile)
            gridfiles.append(cfile)

            for i in range(ntestpergrid):
                tries=20
                while True:
                    lon= round(uniform(xmin,xmax),5)
                    lat= round(uniform(ymin,ymax),5)
                    if land_areas is not None:
                        from shapely.geometry import Point
                        if not land_areas.contains(Point(lon,lat)):
                            tries -= 1
                            if tries > 0:
                                continue
                            break
                    testpts.append(( lon, lat, round(uniform(0.0,1000.0),3)))
                    break

    testpts=np.array(testpts)
    results=np.zeros((testpts.shape[0],9))
    results[:,0:3]=testpts

    for gfile in gridfiles:
        g=DeformationGrid(gfile)
        gcomp=g.bilinear(testpts[:,0],testpts[:,1])
        ncomp=gcomp.shape[1]
        results[:,3:ncomp+1] += gcomp[:,2:ncomp]

    elp=ellipsoid.grs80
    for p in results:
        lon,lat,hgt=p[0:3]
        de,dn,du=p[3:6]
        xyz=elp.xyz(lon,lat)
        dxyz=elp.enu_axes(lon,lat).transpose().dot([de,dn,du])
        xyz += dxyz
        lon,lat,h0=elp.geodetic(xyz)
        p[6:9]=[lon,lat,hgt+du]
        
    testfile=re.sub(r'(\.def)?$','_testpts.csv',def_file,re.I)
    with open(testfile,"wb") as tf:
        tf.write("lon,lat,hgt,de,dn,du,lon1,lat1,hgt1\r\n")
        np.savetxt(tf,results,
                   fmt='%.5f %.5f %.3f %.4f %.4f %.4f %.8f %.8f %.4f'.split(),
                   delimiter=',',newline='\r\n')
    filelist.append(testfile)

    # Build NTv2 grid files (.asc, .gsb)

    ntv2base=re.sub(r'(\.def)?$','',def_file,re.I)
    call([
        'python',
        ntv2py,
        '-b',basedir,
        '-v','20130801',
        ntv2base
        ])
    filelist.append(ntv2base+'.asc')
    filelist.append(ntv2base+'.gsb')
    
    zipfile=re.sub(r'(\.def)?$','.zip',def_file,re.I)
    zipfile=os.path.basename(zipfile)
    if os.path.exists(zipfile):
        os.remove(zipfile)
    command=['zip','-m','-j',zipfile]
    command.extend(filelist)
    call(command)

if __name__=='__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser('Create publishable shift zip files from reverse_patch files')
    parser.add_argument('def_files',nargs='+',help='List of shift .def files from which to construct zip files')
    parser.add_argument('--patch-dir','-p',default='reverse_patch',help='Directory to scan for .def files')
    parser.add_argument('--land-areas','-l',help='File containing a WKT mutipolygon definition of the land areas - used to restrict test points')
    parser.add_argument('--n-test-points','-n',type=int,default=ntestpergrid,help='Number of test points per grid area')

    args=parser.parse_args()

    if args.land_areas:
        load_land_areas( args.land_areas )
    ntestpergrid=args.n_test_points

    deffiles=[]
    if len(args.def_files) > 0:
        deffiles.extend(args.def_files)
    else:
        for dirpath,dirnames,filenames in os.walk(args.patch_dir):
            for f in filenames:
                if f.endswith('.def') and not f.startswith('.'):
                    deffiles.append(os.path.join(dirpath,f))
    for f in deffiles:
        print "Building shift grids for "+f
        build_dump(f)
