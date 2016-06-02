import numpy as np
import matplotlib.mlab as ml
import matplotlib.pyplot as plt
import math
from collections import namedtuple

'''
Functions for processing deformation grid to extract useful information!
'''

radius_of_earth=6400000.0
deg_to_metres=radius_of_earth*math.pi/180.0

class ContourLine(namedtuple('ContourLine','level line')):

    def wkt(self):
        from shapely.geometry import LineString
        return LineString(self.line).to_wkt()

class ContourPoly(namedtuple('ContourPoly','level poly')):

    def wkt(self):
        from shapely.geometry import Polygon
        return Polygon(self.poly[0],self.poly[1]).to_wkt()


def stats( arr, noprint=False, percentiles=None):
    '''
    Print summary statistics for supplied values
    '''
    percentiles = percentiles or (1,2.5,5.10,25,50,75,90,95,97.5,99 )
    results = [
        ('Mean',np.mean(arr)),
        ('Std.dev',np.std(arr)),
        ('Minimum',np.min(arr)),
        ('Maximum',np.max(arr)),
    ]
    pv = np.percentile(arr,percentiles)
    for p,v in zip(percentiles,pv):
        results.append(('Percentile {0:.1f}%'.format(p),v))
    if not noprint:
        for p, v in results:
            print '{0:20s} {1:12f}'.format(p+':',v)
    return results

class Grid( object ):
    '''
    Base Grid class, provides some plotting and analysis functions 
    for a grid
    '''

    def __init__( self, griddata, columns=None, source='' ):
        if len(griddata.shape) != 3:
            raise ValueError('Invalid shape of array for grid data (must have 3 dimensions)')

        self.array = griddata
        self.columns = columns
        self.source = source
    
    def getIndex( self, col ):
        if type(col) == int:
            if col < 0:
                col = self.array.shape[-1]+col
            return col
        if isinstance(col,basestring) and self.columns:
            return self.columns.index(col)

    def _calcLevels( self, data, levels, percentiles ):
        if percentiles:
            if not isinstance(percentiles,(list,tuple)):
                percentiles = list(np.arange(0.0,100.0,percentiles))
                if percentiles[-1] < 100.0:
                    percentiles.append(100.0)
            levels=np.percentile(data,percentiles)
        if not levels:
            levels = (0.003,0.01,0.03,0.1,0.3,1.0,3.0)
        return levels

    def contour( self, col=-1, colours=None, levels=None, percentiles=None, returnContours=False ):
        index = self.getIndex(col)
        label=self.columns[index] if self.columns else 'Col '+str(index)
        data = self.array[:,:,index]
        levels = self._calcLevels( data, levels, percentiles )
        if not colours:
            colours = ('cyan','blue')
        lon=self.array[:,:,0]
        lat=self.array[:,:,1]
        cs=plt.contour(lon,lat,data,levels,colors=colours,label=label)
        result = levels
        if returnContours:
            id = 0
            result = []
            for i, line in enumerate(cs.collections):
                level = levels[i]
                for path in line.get_paths():
                    result.append(ContourLine(level,path.vertices))
        return result

    def contourf( self, col=-1, colours=None, levels=None, percentiles=None, multiple=None, returnContours=False ):
        index = self.getIndex(col)
        label=self.columns[index] if self.columns else 'Col '+str(index)
        data = self.array[:,:,index]
        if multiple:
            data = data * multiple
        levels = self._calcLevels( data, levels, percentiles )
        if not colours:
            colours = ('cyan','blue')
        lon=self.array[:,:,0]
        lat=self.array[:,:,1]
        cs=plt.contourf(lon,lat,data,levels,colors=colours,label=label)
        result = levels
        if returnContours:
            id = 0
            result = []
            for i, c in enumerate(cs.collections):
                level = levels[i]
                for path in c.get_paths():
                    path.should_simplify=False
                    poly = path.to_polygons()
                    if len(poly[0]) < 3: 
                        continue
                    result.append(ContourPoly(level,[poly[0],poly[1:]]))
        return result

    def dumpwkt( self, filename, polygon=True ):
        '''
        Write the grid file to a WKT file as a set of rectangular polygons or line strings
        '''
        if polygon:
            start='POLYGON(('
            end='))\n'
        else:
            start='LINESTRING('
            end=')\n'

        with open(filename,'w') as f:
            f.write('wkt\n')
            for i in range(self.array.shape[0]-1):
                for j in range(self.array.shape[1]-1):
                    y0,y1 = self.array[i:i+2,j,1]
                    x0,x1 = self.array[i,j:j+2,0]
                    f.write(start)
                    f.write('{0} {1}'.format(x0,y0))
                    f.write(',{0} {1}'.format(x0,y1))
                    f.write(',{0} {1}'.format(x1,y1))
                    f.write(',{0} {1}'.format(x1,y0))
                    f.write(',{0} {1}'.format(x0,y0))
                    f.write(end)

    def writeWkt( self, contours, wktfile ):
        '''
        Write contours generated by contour or contourf to a wkt file
        '''
        id = 0
        with open(wktfile,'w') as f:
            f.write("id|level|wkt\n")
            for contour in contours:
                id += 1;
                f.write("{0:d}|{1:f}|{2:s}\n".format(id,contour.level,contour.wkt()))

    def regionsExceedingLevel( self, col, limit, multiple=None ):
        from shapely.geometry import Polygon
        data = self.column(col)
        if multiple:
            cmax=np.max(data*multiple)
        else:
            cmax=np.max(data)
        cmax = cmax*2
        contours = self.contourf(col,levels=[limit,cmax],multiple=multiple, returnContours=True)
        result = []
        for c in contours:
            if c.level == limit:
                result.append(Polygon(c.poly[0],[]))
        return result

    def column( self,col ):
        '''
        Returns the data for a column (specified either by column name 
        or number)
        '''
        return self.array[:,:,self.getIndex(col)]

    def colstats( self, *cols, **params ):
        '''
        Print summary stats for one or more columns. Params are passed
        to stats function (most useful is percentiles).
        '''
        if self.source:
            print self.source
        if not cols:
            cols = range(2,self.array.shape[2])
        for col in cols:
            index = self.getIndex(col)
            label=self.columns[index] if self.columns else 'Col '+str(index)
            print "Statistics for "+label 
            stats(self.array[:,:,index],**params)

    def writecsv( self, filename, delim=',' ):
        '''
        Dump the grid as a CSV file
        '''
        with open(filename,'w') as f:
            f.write(delim.join(self.columns))
            f.write('\n');
            for r in self.array:
                for c in r:
                    f.write(delim.join((str(v) for v in c)))
                    f.write('\n')

    def nodes( self ):
        shape=self.array.shape
        return self.array.reshape(shape[0]*shape[1],shape[2])


class DeformationGrid( Grid ):
    '''
    DeformationGrid class is a grid loaded from a grid file generated by calc_okada
    with some additional functions defining for analysing the 
    modelled deformation.

    Assumes first columns are lon,lat,de,dn,du
    '''

    def __init__( self, filename ):
        '''
        Load the grid from a specified file
        '''
        delimiter=None
        with open(filename,'r') as f:
            line=f.readline()
            if "," in line:
                delimiter=","
            elif "\t" in line:
                delimiter="\t"
            columns=line.lower().strip().split(delimiter)
        if columns[:2] != ['lon','lat']:
            raise RuntimeError('DeformationGrid: input file '+filename+' must have first two columns lon, lat')
        g = np.loadtxt(filename,skiprows=1,dtype=float,delimiter=delimiter)
        rows = ml.find(g[:,0]==g[0,0])
        cols = g.shape[1]
        g.shape=(rows.shape[0],rows[1],cols)
        Grid.__init__(self,g,columns=columns,source=filename)
        self.dln=(g[0,-1,0]-g[0,0,0])/(g.shape[1]-1)
        self.dlt=(g[-1,0,1]-g[0,0,1])/(g.shape[0]-1)
        self.extents=np.array([g[0,0,0:2],g[-1,-1,0:2]])



    def bilinear( self, x, y ):
        '''
        Evaluate the grid at a set of x y values where x,y can each be
        lists of longitudes, latitudes
        '''
        grdx=(np.array(x)-self.extents[0,0])/self.dln
        grdy=(np.array(y)-self.extents[0,1])/self.dlt
        g=self.array
        valid=(grdx >= 0.0) & (grdx < (g.shape[1]-1)) & (grdy >= 0) & (grdy < (g.shape[0]-1))
        nx=np.where(valid,grdx.astype(int),0)
        ny=np.where(valid,grdy.astype(int),0)
        # print "extents",self.extents
        # print "shape",g.shape
        # print "nx",nx
        # print "ny",ny
        # print "valid",valid
        # print "gx",grdx
        # print "gy",grdy
        # print "nodes 1",g[nx,ny]
        # print "nodes 2",g[nx+1,ny]
        # print "nodes 3",g[nx,ny+1]
        # print "nodes 4",g[nx+1,ny+1]
        grdx -= nx
        grdy -= ny
        grdx=grdx.reshape(grdx.shape[0],1)
        grdy=grdy.reshape(grdy.shape[0],1)
        valid=valid.reshape(valid.shape[0],1)
        interp=((g[ny,nx]*(1-grdx)+g[ny,nx+1]*grdx)*(1-grdy)+
               (g[ny+1,nx]*(1-grdx)+g[ny+1,nx+1]*grdx)*grdy)
        empty=np.zeros((grdx.shape[0],g.shape[2]))
        # print "interp",interp
        # print "empty",empty
        return np.where(valid,
                        (g[ny,nx]*(1-grdx)+g[ny,nx+1]*grdx)*(1-grdy)+
                        (g[ny+1,nx]*(1-grdx)+g[ny+1,nx+1]*grdx)*grdy,
                        empty)

    def strainComponents( self ):
        '''
        Calculate horizontal strain components: returns a with columns
           dilatation (linear)
           rotation 
           shear
           distortion

        Max and min scale change are dilatation +/- shear
        Max and min bearing change are rotation +/- shear
        Max angle change = 2*shear

        "Distortion" is the maximum length of vector change to a unit vector.

        All are in units of ppm

        Calculation is done at the centre of each grid cell based on the movement calculated
        at the midpoints of each side.  
        '''
        from ellipsoid import grs80
        g=self.array
        lats=g[:,0,1]
        midlats=(lats[:-1]+lats[1:])/2.0
        de,dn=grs80.metres_per_degree(0,midlats)
        de *= self.dln
        dn *= self.dlt
        de=de.reshape((de.size,1,1))
        dn=dn.reshape((dn.size,1,1))
        den=g[:,:,2:4]
        # Values at middle of E/W sides of each grid cell, then take difference
        # and divide by width of cell
        ddx=(den[:-1,:,:]+den[1:,:,:])/2
        ddx=ddx[:,1:,:]-ddx[:,:-1,:]
        ddx /= de
        dxx=ddx[:,:,0]
        dyx=ddx[:,:,1]
        # Same for latitude (N/S sides of each grid cell)
        ddy=(den[:,:-1,:]+den[:,1:,:])/2
        ddy=ddy[1:,:,:]-ddy[:-1,:,:]
        ddy /= dn
        dxy=ddy[:,:,0]
        dyy=ddy[:,:,1]
        dil=(dxx+dyy)/2.0
        rot=(dxy-dyx)/2.0
        shear=np.sqrt(((dxx-dyy)/2.0)**2 + ((dxy+dyx)/2)**2)
        A=dxx*dxx+dyx*dyx
        B=dxy*dxy+dyy*dyy
        C=dxx*dxy+dyx*dyy
        distortion=np.sqrt((A+B)/2+np.sqrt(((A-B)/2)**2+C**2))
        newshape=(dil.shape[0],dil.shape[1],1)
        dil=dil.reshape(newshape)
        rot=rot.reshape(newshape)
        shear=shear.reshape(newshape)
        distortion=distortion.reshape(newshape)
        dil *= 1000000.0
        rot *= 1000000.0
        shear *= 1000000.0
        distortion *= 1000000.0

        midltln=(g[1:,1:,:2]+g[1:,:-1,:2]+g[:-1,1:,:2]+g[:-1,:-1,:2])/4
        result=np.concatenate((midltln,dil,rot,shear,distortion),axis=2)
        return result

    def calcResolution( self, tolerance, maxsize=100000.0,precision=0.0001,margin=1 ):
        '''
        For all the internal nodes in the grid calculate the misfit 
        of the node from a value calculated from the four adjacent 
        nodes, and the size of grid cells required to keep this less
        than tolerance (up to maxsize).

        Precision is the numeric precision of the grid data. This is used to manage
        the impact of rounding errors on calculating the grid spacing.
        '''

        coslat=math.cos(math.radians(((self.extents[0,0]+self.extents[0,1])/2)))
        gridsize = np.hypot(self.dln*coslat,self.dlt)*deg_to_metres
        minerror = tolerance/(maxsize/gridsize)**2

        # Calculate maximum size that we can reliably determine given a precision
        # of 0.0001 in the deformation (so want difference of 0.0005 to get 
        # reasonable estimate of grid size).

        ms = np.sqrt(tolerance/(5*precision))*gridsize;
        offsets=[margin]

        maxmargin = min(self.array.shape)/4-2
        while margin < maxsize/ms and margin < maxmargin: 
            margin *= 2
            offsets.append(margin)
        # print "Offset required to calculate resolution",offsets[-1]," (gridsize ",gridsize*offsets[-1],")"

        gsgrid=np.empty((self.array.shape[0]-2*margin,self.array.shape[1]-2*margin,4))
        gsgrid[:,:,0:2]=self.array[margin:-margin,margin:-margin,0:2]
    
        en = self.array[:,:,2:4]
        rmax = en.shape[0];
        cmax = en.shape[1];
        gse = gsgrid[:,:,2]
        gsr = gsgrid[:,:,3]
        gsif = -1.0
        
        endiffa = np.empty(gsgrid.shape)
        endiff = endiffa[:,:,0:2]
        error = endiffa[:,:,2]
        res = endiffa[:,:,3]

        for ofs in offsets:
            endiff[:,:,:] = (en[margin:-margin,margin:-margin]-(
                en[margin-ofs:-(margin+ofs),margin:-margin]+
                en[(margin+ofs):(rmax+ofs-margin),margin:-margin]+
                en[margin:-margin,margin-ofs:-(margin+ofs)]+
                en[margin:-margin,(margin+ofs):(cmax+ofs-margin)]
                )/4)
            error[:]=np.hypot(endiff[:,:,0],endiff[:,:,1])/(ofs*ofs)
            res[:]=gridsize*np.sqrt(tolerance/np.maximum(error,minerror))
            if ofs == offsets[0]:
                gse[:] = error
                gsr[:] = res
            else:
                gse[:] = np.where(error < gsif,error,gse)
                gsr[:] = np.where(error < gsif,res,gsr)
            gsif=precision*5/(ofs*ofs)
        return Grid(gsgrid,columns=('lon','lat','error','reqsize'),
            source='Grid size to achieve tolerance '+str(tolerance))

    def grid_definition( self ):
        nlon,nlat=self.array.shape[0:2]
        return ("grid:"+":".join([str(x) for x in self.extents.reshape((4,))])+':'+
                str(self.array.shape[1]-1) + ":" + str(self.array.shape[0]-1))

        

if __name__=="__main__":
    import sys
    for file in sys.argv[1:]:
        g=DeformationGrid(file)
        print file,g.grid_definition()

    
