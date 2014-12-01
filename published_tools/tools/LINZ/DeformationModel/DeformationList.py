import numpy as np

from Error import ModelDefinitionError, UndefinedValueError

class DeformationList( object ):
    '''
    A list of deformation values from which deformations are calculated by
    taking a weighted sum of values. 
    
    The model may have 1 to 5 components.  This may include missing 
    (ie undefined) values.
    '''

    DeformationColumns = ['de','dn','du','eh','ev']

    def __init__( self, columns, size, data=None ):
        for c in columns:
            if c not in DeformationList.DeformationColumns:
                raise ModelDefinitionError('Invalid column '+c+' of deformation data')
        if not c:
            raise ModelDefinitionError('No deformation columns defined')
        if size < 1:
            raise ModelDefinitionError('Invalid size of deformation data list')

        self._nread = 0
        self._columns = list(columns)
        self._columnmapping=[columns.index(c) if c in columns else -1 for c in DeformationList.DeformationColumns]
        self._errorcolumns = []
        if 'eh' in columns: self._errorcolumns.append(columns.index('eh'))
        if 'ev' in columns: self._errorcolumns.append(columns.index('ev'))

        self._size = size
        self._dimension = len(columns)
        self._data = None
        if data:
            self.setData(data)

    def addPoint( self, values ):
        if self._data is None:
            self._data = np.empty( [self._size,self._dimension], float )
        if self._nread >= self._size:
            raise ModelDefinitionEror('Too many data points supplied')
        if len(values) != self._dimension:
            raise ModelDefinitionError('Incorrect number of components at data point')
        self._data[self._nread,:] = values
        self._nread += 1
        # Convert errors to variances as these will be used in calculations
        if self._nread == self._size:
            for ic in self._errorcolumns:
                self._data[:,ic] *= self._data[:,ic]

    def setData( self, data ):
       if data.shape != (self._size,self._dimension):
            raise ModelDefinitionError("Deformation list: supplied data the wrong shape")
       self._data = data
       self._nread=self._size

    def data( self ):
        self.checkValid()
        return self._data;

    def checkValid( self):
        if self._nread != self._size:
            raise ModelDefinitionError('Too few data points supplied ('+str(self._nread)+' instead of '+str(self._size)+')')

    def calcDeformation( self, rows, factors ):
        '''
        Calculate the deformation as a weighted sum of rows of the model.
        rows    - a list of row numbers to use
        factors - a list of weights to apply to the values
        '''
        assert self._nread == self._size
        value = self._data[rows[0],:] * factors[0];
        for i in range(1,len(rows)):
            value += self._data[rows[i],:] * factors[i];
        for v in value:
            if np.isnan(value).any():
                raise UndefinedValueError('The deformation is undefined at this location')
        return [value[c] if c >= 0 else 0.0 for c in self._columnmapping]

