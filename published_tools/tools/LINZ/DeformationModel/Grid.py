from __future__ import with_statement

from CsvFile import CsvFile
import numpy as np
import math
import os.path

from DeformationList import DeformationList
from Error import ModelDefinitionError, UndefinedValueError, OutOfRangeError

class Grid( object ):
    '''
    Module to define and use a grid deformation model.
    '''
    def __init__( self, model, gridfile, minlon, maxlon, minlat, maxlat, nlon, nlat, columns, name=None ):
        if not os.path.exists(model.getFileName(gridfile)):
            raise ModelDefinitionError("Invalid grid filename "+str(gridfile))
        self._model = model
        self._gridfile = gridfile
        name = name if name else gridfile
        self._name = name
        self._columns = list(columns)
        self._dimension=len(columns)
        self._fields=['lon float','lat float']
        self._fields.extend('data[]='+c+' ?float' for c in columns)

        self._nlon = int(nlon)
        self._nlat = int(nlat)
        if self._nlon < 2 or self._nlat < 2:
            raise ModelDefinitionError("Invalid number of grid rows or columns in deformation model definition for "+name)

        self._minlon = float(minlon)
        self._maxlon = float(maxlon)
        self._dlon = (self._maxlon-self._minlon)/(self._nlon-1)
        if self._dlon < 0:
            raise ModelDefinitionError("Invalid longitude range "+str(minlon)+" - "+str(maxlon)+" in deformation model definition for "+name)

        self._minlat = float(minlat)
        self._maxlat = float(maxlat)
        self._dlat = (self._maxlat-self._minlat)/(self._nlat-1)
        if self._dlat < 0:
            raise ModelDefinitionError("Invalid latitude range "+str(minlat)+" - "+str(maxlat)+" in deformation model definition for "+name)

        self._npt = self._nlon * self._nlat
        self._loaded = False
        self._valid = False
        self._data = DeformationList( columns, self._npt )

    def gridSpec( self ):
        '''
        Returns min lon, min lat, max lon, max lat, nlon, nlat
        '''
        return self._minlon, self._minlat, self._maxlon, self._maxlat, self._nlon, self._nlat

    def resolution( self ):
        '''
        Returns the longitude and latitude increments
        '''
        return self._dlon, self._dlat

    def gridFile( self ):
        return self._gridfile

    def data( self ):
        return self._data

    def load( self ):
        if self._loaded:
            return
        try:
            gridmetadata = [self._nlon,self._nlat]
            gridmetadata.extend(self._columns)
            data = self._model.cacheData( self._gridfile, gridmetadata )
            if data:
                self._data.setData(data)
                self._valid = True
            else:
                with CsvFile('grid',self._model.getFileName(self._gridfile),self._fields) as f:
                    nc = -1
                    nr = 0
                    lontol = self._dlon/10000.0
                    lattol = self._dlat/10000.0
                    xc = self._minlon-self._dlon
                    yc = self._minlat
                    for gpt in f:
                        nc += 1
                        xc += self._dlon
                        if nc >= self._nlon:
                            nc = 0
                            xc = self._minlon
                            nr += 1
                            yc += self._dlat
                            if nr > self._nlat:
                                raise ModelDefinitionError("Too many grid points in "+self._name)
                        if (math.fabs(gpt.lon-xc) > lontol or math.fabs(gpt.lat-yc) > lattol):
                            raise ModelDefinitionError("Grid latitude/longitude out of sequence: ("+str(gpt.lon)+','+str(gpt.lat)+') should be ('+str(xc)+','+str(yc)+') in '+self._name)

                        self._data.addPoint(gpt.data)
                    self._data.checkValid()
                    self._valid = True
                    self._model.setCacheData( self._data.data(), self._gridfile, gridmetadata )
        finally:
            self._loaded = True

    def calcDeformation( self, x, y ):
        '''
        Calculate the deformation at a specific cpoint
        '''
        if not self._loaded:
            self.load()
        if not self._valid:
            raise ModelDefinitionError("Cannot use invalid grid component - see previous errors")

        x0 = x
        while x < self._minlon:
            x += 360
        if (x > self._maxlon or y < self._minlat or y > self._maxlat ):
            raise OutOfRangeError(str(x0)+','+str(y)+' is out of range of grid in '+self._name)
        wx = (x-self._minlon)/self._dlon
        wy = (y-self._minlat)/self._dlat
        nx = int(wx)
        ny = int(wy)
        if nx >= self._nlon-1: nx = self._nlon-2
        if ny >= self._nlat-1: ny = self._nlat-2
        wx -= nx
        wy -= ny
        ny *= self._nlon
        rows=(nx+ny,nx+ny+1,nx+ny+self._nlon,nx+ny+self._nlon+1)
        factors=((1-wx)*(1-wy),wx*(1-wy),(1-wx)*wy,wx*wy)
        return self._data.calcDeformation(rows,factors)
