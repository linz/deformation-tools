from __future__ import with_statement

import numpy as np
import csv
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
        self._columns[0:0]=['lon','lat']

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


    def load( self ):
        if self._loaded:
            return
        try:
            gridmetadata = [self._nlon,self._nlat].extend(self._columns)
            data = self._model.cacheData( self._gridfile, gridmetadata )
            if data:
                self._data.setData(data)
                self._valid = True
            else:
                with open(self._model.getFileName(self._gridfile),"rb") as f:
                    c = csv.reader(f)
                    header = c.next()
                    if header != self._columns:
                        raise ModelDefinitionError("Invalid grid model header "+','.join(header)+' in '+self._name + ' (expected '+','.join(self._columns)+')')
                    nc = -1
                    nr = 0
                    lontol = self._dlon/10000.0
                    lattol = self._dlat/10000.0
                    xc = self._minlon-self._dlon
                    yc = self._minlat
                    for line in c:
                        nc += 1
                        xc += self._dlon
                        if nc >= self._nlon:
                            nc = 0
                            xc = self._minlon
                            nr += 1
                            yc += self._dlat
                            if nr > self._nlat:
                                raise ModelDefinitionError("Too many grid points: "+','.join(line)+' in '+self._name)
                        parts = [float(x) if x != '' else None for x in line]
                        x, y = parts[0:2]
                        if (math.fabs(x-xc) > lontol or math.fabs(y-yc) > lattol):
                            raise ModelDefinitionError("Grid latitude/longitude out of sequence: "+','.join(line)+' should be ('+str(xc)+','+str(yc)+') in '+self._name)

                        dvec = parts[2:]
                        if len(dvec) != self._dimension:
                            raise ModelDefinitionError("Missing grid displacement components in record: "+','.join(line)+' in '+self._name)
                        self._data.addPoint(dvec)
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

        if (x < self._minlon or x > self._maxlon or 
            y < self._minlat or y > self._maxlat ):
            raise OutOfRangeError(str(x)+','+str(y)+' is out of range of grid in '+self._name)
        wx = (x-self._minlon)/self._dlon
        wy = (y-self._minlat)/self._dlat
        nx = int(wx)
        ny = int(wy)
        if nx >= self._nlon: nx = self._nlon-1
        if ny >= self._nlat: ny = self._nlat-1
        wx -= nx
        wy -= ny
        ny *= self._nlon
        rows=(nx+ny,nx+ny+1,nx+ny+self._nlon,nx+ny+self._nlon+1)
        factors=((1-wx)*(1-wy),wx*(1-wy),(1-wx)*wy,wx*wy)
        return self._data.calcDeformation(rows,factors)
