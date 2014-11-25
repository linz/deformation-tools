from __future__ import with_statement

from CsvFile import CsvFile
import numpy as np
import math
import os.path

from DeformationList import DeformationList
from Error import ModelDefinitionError, UndefinedValueError, OutOfRangeError

class TIN( object ):
    '''
    Calculates a triangulated irregular network (TIN) deformation model.
    '''
    pts_fields= ['id int','lon float','lat float']
    trg_fields= ['ids[]=id1 int','ids[]=id2 int','ids[]=id3 int']
    
    def __init__( self, model, trigptsfile,trigtrgfile, minlon, maxlon, minlat, maxlat, npt, ntrg, columns, name=None ):
        '''
        Load the definition of a TIN model.  The model is not loaded until either it is required to
        calculate a value or it is explicitly loaded using the load function.
        '''
        if not os.path.exists(model.getFileName(trigptsfile)):
            raise ModelDefinitionError("Invalid trig point filename "+str(trigptsfile))
        if not os.path.exists(model.getFileName(trigtrgfile)):
            raise ModelDefinitionError("Invalid trig triangulation filename "+str(trigtrgfile))
        self._model = model
        self._trigptsfile = trigptsfile
        self._trigtrgfile = trigtrgfile
        name = name if name else trigptsfile
        self._name = name

        self._columns = columns
        self._dimension = len(columns)
        self._fields=list(TIN.pts_fields)
        self._fields.extend('data[]='+x+' float' for x in columns)

        self._npt = int(npt)
        self._ntrg = int(ntrg)
        if self._npt < 2 or self._ntrg < 1:
            raise ModelDefinitionError("Invalid number of triangulation points or triangles in deformation model definition for "+name)

        self._minlon = float(minlon)
        self._maxlon = float(maxlon)
        if self._minlon >= self._maxlon:
            raise ModelDefinitionError("Invalid longitude range "+str(minlon)+" - "+str(maxlon)+" in deformation model definition for "+name)

        self._minlat = float(minlat)
        self._maxlat = float(maxlat)
        if self._minlat >= self._maxlat:
            raise ModelDefinitionError("Invalid latitude range "+str(minlat)+" - "+str(maxlat)+" in deformation model definition for "+name)

        self._loaded = False
        self._valid = False
        self._data = None
        self._points = None
        self._trg = None
        self._centroids = None
        self._adjacent = None
        self._edgevec = None


    def load( self ):
        '''
        Load the TIN deformation model
        '''
        if self._loaded:
            return
        self._points = np.empty([self._npt+1,2], float)
        self._trg = np.empty([self._ntrg,3],int)
        self._data = DeformationList( self._columns, self._npt+1 )
        self._data.addPoint([0]*self._dimension)
        npt = 0
        ntrg = 0
        name = self._name
        model=self._model
        try:
            trigptsmetadata = [npt]
            data = self._model.cacheData( self._trigptsfile, trigptsmetadata )
            if data:
                self._data.setData(data)
            else:
                with CsvFile('trig_pts',model.getFileName(self._trigptsfile),self._fields) as f:
                    npt = 0
                    for tpt in f:
                        npt += 1
                        if npt > self._npt:
                            raise ModelDefinitionError("Too many points in triangulation model"+','.join(line)+' for '+name)
                        if tpt.id != npt:
                            raise ModelDefinitionError("TIN point id out of sequence: "+','.join(line)+' for '+name)
                        if (tpt.lon < self._minlon or tpt.lon > self._maxlon or
                            tpt.lat < self._minlat or tpt.lat > self._maxlat):
                            raise ModelDefinitionError("TIN latitude/longitude out of range: "+','.join(line)+' for '+name)
                        self._data.addPoint(tpt.data)
                        self._points[npt:] = [tpt.lon,tpt.lat]
                self._data.checkValid()
                self._model.setCacheData(self._data.data(), self._trigptsfile, trigptsmetadata )

            trigtrgmetadata = [npt]
            data = self._model.cacheData( self._trigtrgfile, trigtrgmetadata )
            if data:
                if data.shape != (self._ntrg,3):
                    raise ModelDefinitionError("Cached triangulation has wrong shape")
                self._trg = data
            else:
                with CsvFile('trig_trg',model.getFileName(self._trigtrgfile),TIN.trg_fields) as f:
                    ntrg = 0
                    for trg in f:
                        if ntrg >= self._ntrg:
                            raise ModelDefinitionError("Too many triangles in triangulation model: "+','.join(line)+' for '+name)
                        for id in trg.ids:
                            if id < 1 or id > self._npt:
                                raise ModelDefinitionError("Invalid triangle point id "+str(id)+': '+','.join(line))
                        self._trg[ntrg:]=trg.ids
                        ntrg += 1
                if ntrg != self._ntrg:
                    raise ModelDefinitionError("Not enough triangle definitions in trig file - expected "+str(self._ntrg)+' found '+ str(ntrg) + ' for '+name)
                self._model.setCacheData( self._trg, self._trigtrgfile, trigtrgmetadata )

            self._setupTriangulation()
            self._valid = True
        finally:
            self._loaded = True

    def calcDeformation( self, x, y ):
        '''
        Calculate the deformation at a specific point in the TIN
        '''
        if not self._loaded:
            self.load()
        if not self._valid:
            raise ModelDefinitionError("Cannot use invalid TIN component - see previous errors")

        x0 = x
        while x < self._minlon:
            x += 360
        if (x < self._minlon or x > self._maxlon or 
            y < self._minlat or y > self._maxlat ):
            raise OutOfRangeError(str(x0)+','+str(y)+' is out of range of TIN for '+self._name)
        ntrg, rows, factors = self.findTriangle(x,y)
        return self._data.calcDeformation(rows,factors)

    def _setupTriangulation( self ):
        # Have we cached this?

        metadata = [len(self._points),len(self._trg)]
        files = [self._trigptsfile, self._trigtrgfile]
        self._centroids = self._model.cacheData(self._trigptsfile+'.centroids',metadata,files=files)
        self._edgevec = self._model.cacheData(self._trigtrgfile+'.edgevec',metadata,files=files)
        self._adjacent = self._model.cacheData(self._trigtrgfile+'.adjacent',metadata,files=files)
        if self._centroids and self._edgevec and self._adjacent:
            return

        # First check that triangles are all anticlockwise

        pts = self._points
        trg = self._trg
        name = self._name
        areas = np.cross(pts[trg[:,1]]-pts[trg[:,0]], pts[trg[:,2]]-pts[trg[:,0]])

        badtrg = np.flatnonzero(areas < 0)
        if len(badtrg):
            raise ModelDefinitionError(str(len(badtrg))+' of '+str(len(areas))+
                  ' triangles are clockwise eg '+str(trg[badtrg[0],:])+' with area '+str(areas[badtrg[0]])+
                  ' in '+name)


        # Now check edges
        # Also construct adjacent triangle list - identifying adjacent triangle
        # across each edge of triangle opposite the corresponding node
        # or -1 if there is none.

        self._adjacent = np.zeros_like(trg)-1

        edgelist = {}
        for nt in range(len(trg)):
            t = trg[nt,:]
            for i in range(0,3):
                v=str(t[(i+1)%3])+' '+str(t[(i+2)%3])
                if v in edgelist:
                    raise ModelDefinitionError("Edge "+v+" repeated in triangulation definition")
                edgelist[v]=[nt,i if i >= 0 else 2 ]
        boundary = {}
        nedge = 0
        start = 0
        for e in edgelist:
            vr = ' '.join(reversed(e.split()))
            if vr in edgelist:
                nt,ne = edgelist[e]
                self._adjacent[nt,ne] = edgelist[vr][0]
            else:
                nedge += 1
                f,t = map(int,e.split())
                start = f
                boundary[f] = t
        del edgelist
        nloop = 0
        start0 = start
        p0 = start
        while True:
            p1 = boundary[p0]
            if p1 not in boundary:
                raise ModelDefinitionError("Triangle boundary error at node "+str(p1))
            p2 = boundary[p1]
            area = np.cross(pts[p1]-pts[p2],pts[p1]-pts[p0])
            if area < 0:
                raise ModelDefinitionError("Triangulation boundary concave at node "+str(p1))
            p0 = p1
            nloop += 1
            if p0 == start:
                break
            if nloop >= nedge:
                raise ModelDefinitionError("Invalid triangulation boundary")
        if nloop < nedge:
            raise ModelDefinitionError("Triangulation is not a single convex polygon")
        del boundary
        
        # Now create an array of triangle centroids and scaled vectors normal to the edges

        self._centroids = (pts[trg[:,0]]+pts[trg[:,1]]+pts[trg[:,2]])/3
        self._edgevec = (pts[trg[:,[2,0,1]]] - pts[trg[:,[1,2,0]]])/np.reshape(areas,(-1,1,1))
        self._model.setCacheData(self._centroids, self._trigptsfile+'.centroids',metadata=metadata,files=files)
        self._model.setCacheData(self._edgevec, self._trigtrgfile+'.edgevec',metadata=metadata,files=files)
        self._model.setCacheData(self._adjacent, self._trigtrgfile+'.adjacent',metadata=metadata,files=files)

    def findTriangle( self, x, y ):
        self.load()
        checked = []
        pt = np.array([x,y])
        # Must be a better way!
        start = np.argmin(np.hypot(self._centroids[:,0]-x,self._centroids[:,1]-y))
        
        while True:
            checked.append(start)
            weights = np.cross(self._edgevec[start],pt-self._centroids[start])+1.0/3.0
            next = -1
            for i in range(0,3):
                if weights[i] < 0:
                    next = self._adjacent[start,i]
                    if next < 0:
                        raise OutOfRangeError(str(x)+','+str(y)+' is out of range of triangulation')
                    break
            if next < 0 or next in checked:
                break
            start = next

        rows = self._trg[start]
        return start,rows, weights

    def triangles( self ):
        '''
        Iterate over triangles in the network.  Returns a dictionary of id (triangle id),
        pt_ids (list of point ids), points (list of point coordinates)
        '''
        self.load()
        for id in range(len(self._trg)):
            yield dict( id=id, pt_ids=self._trg[id], points=[self._points[pt] for pt in self._trg[id]])

    def containsPoint( self, x, y ):
        try:
            self.findTriangle(x,y)
            return True
        except:
            return False
            

