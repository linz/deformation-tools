

class Cache( object ):
    '''
    A binary cache for array data to improve load time for deformation models.

    Uses HDF5 storage implemented by PyTables.  If this is not availabe, then it will
    just create a null cache and all data will be reloaded from the ASCII files.

    Each cached array is stored using its filename (relative to the model base
    directory) as an identifier.  It also stores a metadata string which is used to test
    if the cached data is still valid.  Sensibly this will be a combination of file
    attributes (size, last accessed date) and significant metadata (number of rows,
    column names, etc) to ensure that it is still valid.
    '''

    def __init__( self, filename ):
        self._locked = True
        self._h5file = None
        # Try opening with write access ...
        try:
            import tables
            import warnings
            warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)
            self._h5file=tables.openFile(filename,mode="a",title="Deformation model cache")
            self._locked = False
        except ImportError:
            return 
        except IOError:
            pass
        # Failing that, try read only
        if not self._h5file:
            try:
                self._h5file=tables.openFile(filename,mode="r",title="Deformation model cache")
            except IOError:
                pass

    def get( self, filename, metadata ):
        if not self._h5file:
            return None
        grid=None
        try:
            grid=self._h5file.getNode('/'+filename)
            if grid.getAttr('cache_metadata') != metadata:
                grid=None
        except:
            grid=None
        return grid

    def set( self, filename, metadata, value ):
        if self._locked:
            return None
        hf=self._h5file
        if '/'+filename in hf:
            hf.removeNode('/',filename)
        parts=filename.split('/')
        name = parts[-1]
        path = '/'+'/'.join(parts[:-1])
        grid=hf.createArray(path,name,value,createparents=True)
        grid.setAttr('cache_metadata',metadata)
        hf.flush()

    def close( self ):
        if self._h5file:
            self._h5file.close()
            self._h5file = None




    


