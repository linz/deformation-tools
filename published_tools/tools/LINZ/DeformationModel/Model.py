
import os.path
import sys
import logging
import datetime
import time
import math

from Error import ModelDefinitionError, OutOfRangeError, UndefinedValueError
from CsvFile import CsvFile
from Time import Time
from TimeModel import TimeModel
from Grid import Grid
from TIN import TIN
from Cache import Cache


def _buildHash(x,attrs):
    return ':'.join((str(getattr(x,a)) for a in attrs))

class TimeComponent( object ):
    '''
    The time component of a deformation model.  A wrapper around a time model managing also the 
    range of valid values for the model and for caching the calculated value
    '''

    hashattr=['temporal_model','factor0','time0','factor1','time1','decay']

    @staticmethod
    def hashKey( compdef ):
        return _buildHash( compdef, TimeComponent.hashattr )

    @staticmethod
    def get( models, compdef ):
        hash=TimeComponent.hashKey(compdef)
        if hash not in models:
            models[hash] = TimeComponent(compdef)
        return models[hash]

    def __init__( self, compdef ):
        self.hash = self.hashKey(compdef)
        self.min_date = compdef.min_date
        self.max_date = compdef.max_date
        self.temporal_model = compdef.temporal_model
        self.temporal_complete = compdef.temporal_complete
        self.description = compdef.description
        self.calc_date = None
        self.base_calc_date = None
        self.calc_value = None
        self.calc_error_value = None
        self.calc_valid = True
        self.time0=compdef.time0
        self.time1=compdef.time1
        self.factor0=compdef.factor0
        self.factor1=compdef.factor1
        self.decay=compdef.decay
        self._model = TimeModel(compdef.temporal_model,compdef.factor0,compdef.time0,compdef.factor1,compdef.time1,compdef.decay)

    def calcFactor( self, date, baseDate=None ):
        if date != self.calc_date or baseDate != self.calc_base_date:
            self.calc_date=date
            self.calc_base_date = None
            self.calc_value = 0.0
            self.calc_error_value = 0.0
            ok = True

            for d in (baseDate, date):
                if d == None: continue
                factor = 0.0
                if (self.min_date and d < self.min_date) or (self.max_date and d > self.max_date):
                    if not self.temporal_complete:
                        ok = False
                        break
                else:
                    factor = self._model.calcFactor(d)
                    logging.info("Time factor %s calculated at %s for %s",factor,d,self._model)
                self.calc_value = factor - self.calc_value

            if not ok:
                self.calc_value = None
                self.calc_error_value = None
            else:
                self.calc_error_value = abs(self.calc_value)
                if self._model.squareVarianceFactor: 
                    self.calc_error_value *= self.calc_error_value

        if self.calc_value == None:
            raise OutOfRangeError('Date outside valid range')
        return (self.calc_value, self.calc_error_value)

    def model(self):
        return self._model

    def __str__( self ):
        return str(self._model)


class SpatialComponent( object ):
    '''
    Manages the spatial component, TIN or grid.
    
    A spatial component may be invoked by multiple components in the 
    model definition component.csv file.  The hashKey function is used 
    to determine if two components use the same spatial model.  If so 
    then the model just adds time components.  This means the potentially
    expensive spatial component calculation needs only be performed once
    if necessary.

    Spatial and time components are added using a compdef - basically 
    a row from the component.csv file.

    Spatial and time calculations are cached to make calculating a set
    of deformations at the same time or a time series at a fixed location
    efficient.
    '''

    hashattr = ['spatial_model','file1','file2']
    checkattr = ['min_lon','min_lat','max_lon','max_lat','spatial_complete','npoints1','npoints2',
                 'displacement_type']

    @staticmethod
    def hashKey( component, compdef ):
        return component+':'+_buildHash(compdef,SpatialComponent.hashattr)

    @staticmethod
    def compatibleDefinition( compdef1, compdef2 ):
        for attr in SpatialComponent.checkattr:
            if getattr(compdef1,attr) != getattr(compdef2,attr):
                return False
        return True

    @staticmethod
    def get( model, models, component, compdef, load=False ):
        hash = SpatialComponent.hashKey(component, compdef)
        if hash not in models:
            models[hash] = SpatialComponent( model, component, compdef, load )
        else:
            if not SpatialComponent.compatibleDefinition(models[hash],compdef):
                raise ModelDefinitionError('Inconsistent usage of grid/TIN file '+compdef.file1+' in '+component+' component.csv')
        return models[hash]

    def __init__( self, model, component, compdef, load=False ):
        self.hash = SpatialComponent.hashKey(component, compdef)
        self.min_lon=compdef.min_lon
        self.min_lat=compdef.min_lat
        self.max_lon=compdef.max_lon
        self.max_lat=compdef.max_lat
        self.spatial_complete=compdef.spatial_complete
        self.npoints1=compdef.npoints1
        self.npoints2=compdef.npoints2
        self.displacement_type=compdef.displacement_type
        self.error_type=compdef.error_type
        self.columns=[]
        if self.displacement_type in ['horizontal','3d']: self.columns.extend(['de','dn'])
        if self.displacement_type in ['vertical','3d']: self.columns.append('du')
        if self.error_type in ['horizontal','3d']: self.columns.append('eh')
        if self.error_type in ['vertical','3d']: self.columns.append('ev')
        self.spatial_model=compdef.spatial_model
        self.file1=os.path.join(component, compdef.file1)
        self.file2=os.path.join(component ,compdef.file2) if compdef.file2 else None
        self._model = None
        name = os.path.join(component,compdef.file1)
        self._name = name
        if compdef.spatial_model == 'llgrid':
            self._model = Grid(model,self.file1,self.min_lon,self.max_lon,self.min_lat,self.max_lat,self.npoints1,self.npoints2,self.columns,name=name)
            self._description = "Grid model using "+self._name
        elif compdef.spatial_model == 'lltin':
            if not self.file2:
                raise ModelDefinitionError('file2 is not defined for TIN model')
            self.file2 = os.path.join(component,compdef.file2)
            self._model = TIN(model,self.file1,self.file2,self.min_lon,self.max_lon,self.min_lat,self.max_lat,self.npoints1,self.npoints2,self.columns,name=name)
            self._description = "TIN model using "+self._name
        else:
            raise ModelDefinitionError('Invalid spatial model type '+compdef.spatial_model)

        # Cached calculations
        self._xy=(None,None)
        self._xydisp = [0.0]*5
        self._xyInRange = True
        self._xyRangeError=None
        self._defUndefinedError=None

        self._modelError = None

        if load:
            self.load()

    def __str__(self):
        return self._description

    # So that this can be used as a spatial component set
    def models( self ):
        yield self

    def load( self ):
        '''
        Spatial models are loaded on demand by default.  
        
        The load method forces loading of the spatial component 
        for validation.
        '''
        try:
            self._model.load()
        except ModelDefinitionError:
            self._modelError = str(sys.exc_info()[1])
            raise

    def name( self ):
        return self._name

    def model( self ):
        return self._model

    def containsPoint( self, x, y):
        if x < self.min_lon or x > self.max_lon: return False
        if y < self.min_lat or y > self.max_lat: return False
        if self.spatial_model=='lltin': return self._model.containsPoint(x,y)
        return True

    def calcDeformation( self, x, y):
        logging.info("Calculating spatial component %s at (%f,%f)",self._name,x,y)
        if self._modelError:
            raise ModelDefinitionError( self._modelError )
        try:
            if (x,y) != self._xy:
                try:
                    self._xy=(x,y)
                    self._xyRangeError = None
                    self._defUndefinedError = None
                    self._xydisp = self._model.calcDeformation(x,y)
                    self._xyInRange = True
                except OutOfRangeError:
                    if self.spatial_complete:
                        self._xydisp = [0.0]*5
                        self._xyInRange = False
                    else:
                        self._xyRangeError = sys.exc_info()[1]
                        raise
                except UndefinedValueError:
                    self._defUndefinedError = sys.exc_info()[1]
                    raise
                logging.info("Spatial component calculated as %s",self._xydisp)
            elif self._xyRangeError:
                raise OutOfRangeError( self._xyRangeError )
            elif self._defUndefinedError:
                raise UndefinedValueError( self._defUndefinedError )
            else:
                logging.info("Using cached spatial component %s",self._xydisp)
    
            return self._xydisp, self._xyInRange

        except ModelDefinitionError:
            self._modelError = str(sys.exc_info()[1])
            raise

class SpatialComponentSet( object ):
    '''
    Manages a prioritised set of spatial models, such as a nested grid.  

    This provides the same interface as SpatialComponent.
    '''

    checkattr = ['version_added','version_revoked','displacement_type','error_type']

    @staticmethod
    def compatibilityHash( compdef ):
        return (_buildHash(compdef,SpatialComponentSet.checkattr) + ':'
                + _buildHash(compdef,TimeComponent.hashattr))

    def __init__( self, component, model, compdef ):
        self._component = component
        self._subcomponent = compdef.subcomponent
        self._models = [model]
        self._priorities = [compdef.priority]
        self._sortedModels = [model]
        self._baseModel = model
        self._checkhash = SpatialComponentSet.compatibilityHash(compdef)
        self.min_lon=model.min_lon
        self.min_lat=model.min_lat
        self.max_lon=model.max_lon
        self.max_lat=model.max_lat

        # Cached calculations
        self._xy=(None,None)
        self._xydisp = [0.0]*5
        self._xyInRange = True
        self._xyRangeError=None
        self._defUndefinedError=None
        self._modelError = None

    def addComponent( self, model, compdef ):
        if SpatialComponentSet.compatibilityHash(compdef) != self._checkhash:
            raise ModelDefinitionError('Subcomponent '+str(self.subcomponent)+' of '+component+' uses inconsistent versions, time models or displacement/error components')

        self.min_lon = min(self.min_lon,model.min_lon)
        self.min_lat = min(self.min_lat,model.min_lat)
        self.max_lon = max(self.max_lon,model.max_lon)
        self.max_lat = max(self.max_lat,model.max_lat)
        self._models.append(model)
        self._priorities.append(compdef.priority)
        self._sortedModels = [
            self._models[i] for i in 
            sorted(range(len(self._models)),key=lambda j:self._priorities[j],reverse=True)
            ]
        self._baseModel = self._sortedModels[-1]

    def __str__(self):
        if len(self._models) == 1:
            return self._baseModel._description
        description="Nested models:"
        for m in reversed(self._sortedModels):
            description=description+"\n"+str(m)
        return description

    def models(self):
        for m in reversed(self._sortedModels):
            yield m

    def load( self ):
        '''
        Spatial models are loaded on demand by default.  
        
        The load method forces loading of the spatial component 
        for validation.
        '''
        try:
            for model in self._sortedModels:
                model.load()
        except ModelDefinitionError:
            self._modelError = str(sys.exc_info()[1])
            raise

    def name( self ):
        return self._baseModel._name + ' and subcomponents'

    def containsPoint( self, x, y):
        if x < self.min_lon or x > self.max_lon: return False
        if y < self.min_lat or y > self.max_lat: return False
        for model in reversed(self._sortedModels):
            if not model.containsPoint(x,y):
                return False
        return True

    def calcDeformation( self, x, y):
        logging.info("Calculating spatial component %s at (%f,%f)",self.name(),x,y)
        if self._modelError:
            raise ModelDefinitionError( self._modelError )
        try:
            if (x,y) != self._xy:
                try:
                    self._xy=(x,y)
                    self._xyRangeError = None
                    self._defUndefinedError = None
                    for m in self._sortedModels:
                        self._xydisp, self._xyInRange = m.calcDeformation(x,y)
                        if self._xyInRange:
                            break
                except OutOfRangeError:
                    self._xyRangeError = sys.exc_info()[1]
                    raise
                except UndefinedValueError:
                    self._defUndefinedError = sys.exc_info()[1]
                    raise
                logging.info("Spatial component %s calculated as %s",self.name(),self._xydisp)
            elif self._xyRangeError:
                raise OutOfRangeError( self._xyRangeError )
            elif self._defUndefinedError:
                raise UndefinedValueError( self._defUndefinedError )
            else:
                logging.info("Using cached %s spatial component %s",self.name(),self._xydisp)
    
            return self._xydisp, self._xyInRange

        except ModelDefinitionError:
            self._modelError = str(sys.exc_info()[1])
            raise

class ModelComponent( object ):
    '''
    A model component combines a spatial and time component with a range of valid versions.
    The list of model components is used to compile the deformation model for any required version.
    '''

    def __init__( self, component, compdesc, compdef, spatialComp, timeComp ):
        self.component = component
        self.compdesc = compdesc
        self.description = compdef.description
        self.versionAdded = compdef.version_added
        self.versionRevoked = compdef.version_revoked
        self.subcomponent = compdef.subcomponent
        self.priority = compdef.priority
        self.spatialComponent = spatialComp
        self.name = component+'/'+spatialComp.name()
        self.timeComponent = timeComp
        self.factor = 1.0
        self.timeFactor = 0.0
        self.timeErrorFactor = 0.0

    def appliesForVersion( self, version ):
        '''
        Test if the component applies to a specific version.
        '''
        return self.versionAdded <= version and (self.versionRevoked == '0' or self.versionRevoked > version)

    def setFactor( self, factor ):
        self.factor = factor

    def setDate( self, date, baseDate=None ):
        '''
        Calculate the deformation as a triple [de,dn,du] at a specific time and 
        location.  
        
        The time calculation is cached as the most common usage will be for
        many calculations at the same date.
        '''
        logging.info("Setting component %s date %s (base date %s)",self.name,date,baseDate)
        
        self.timeFactor, self.timeErrorFactor = self.timeComponent.calcFactor( date, baseDate )
        self.timeFactor *= self.factor
        self.timeErrorFactor *= self.factor
        logging.info("Time factor calculated as %s",self.timeFactor)

    def calcDeformation( self, x, y ):
        '''
        Calculate the deformation [de,dn,du,vh,vv] at a specified location.  
        Note that vh and vv are the variances horizonal and vertical if defined
        Assumes that the time component has already been calculated by a call to setDate
        '''
        logging.info("Calculating component %s for location (%s,%s)",self.name,x,y)

        # If the time factor is 0 then don't need to do any more
        t0 = self.timeFactor
        if t0 == 0.0:
            logging.info("Time factor = 0.0 - spatial not calculated")
            return [0.0,0.0,0.0,0.0,0.0]

        t1 = self.timeErrorFactor
        value = self.spatialComponent.calcDeformation(x,y)[0]
        return [value[0]*t0,value[1]*t0,value[2]*t0,value[3]*t1,value[4]*t1]

class Model( object ):
    '''
    Defines a deformation model which may have multiple versions and multiple components.  

    The model is loaded by specifying a base directory containing the files defining the model.  
    It can be used to calulate the deformation from any version of the model and at
    any time.  Also it can calculate the difference between two versions.

    The model is defined by a set of CSV (comma separated value) files.
    '''

    versionspec = CsvFile.FieldSpec('version',[
        'version \\d{8}',
        'release_date datetime',
        'reverse_patch boolean',
        'reason unicode'
        ])

    modelspec = CsvFile.FieldSpec('model',[
        'component \\w+',
        'version_added \\d{8}',
        'version_revoked (\\d{8}|0)',
        'reverse_patch boolean',
        'description unicode'
        ])

    metadataspec = CsvFile.FieldSpec('metadata',[
        'item \\w+',
        'value unicode'
        ])

    componentspec = CsvFile.FieldSpec('component',[
        'version_added \\d{8}',
        'version_revoked (\\d{8}|0)',
        'reverse_patch boolean',
        'subcomponent int',
        'priority int',
        'min_lon float',
        'max_lon float',
        'min_lat float',
        'max_lat float',
        'spatial_complete boolean',
        'min_date datetime',
        'max_date datetime',
        'temporal_complete boolean',
        'npoints1 int',
        'npoints2 int',
        'displacement_type (horizontal|vertical|3d|none)',
        'error_type (horizontal|vertical|3d|none)',
        'max_displacement float',
        'spatial_model (llgrid|lltin)',
        'temporal_model (velocity|step|ramp|decay)',
        'time0 ?datetime',
        'factor0 ?float',
        'time1 ?datetime',
        'factor1 ?float',
        'decay ?float',
        'file1 \w+\.csv',
        'file2 ?\w+\.csv',
        'description unicode',
        ])

    metadataitems = '''
        model_name
        description
        version
        datum_code
        datum_name
        datum_epoch
        datum_epsg_srid
        ellipsoid_a
        ellipsoid_rf
        authority
        authority_website
        authority_address
        authority_email
        source_url
        '''.split()

    def __init__( self, basedir, version=None, baseVersion=None, loadComponent=None, load=False, useCache=True, clearCache=False ):
        '''
        Loads the deformation model located at the specified base directory (the
        directory holding the model.csv file).  If load=True then the spatial components
        of all models are preloaded (ie TINs and grids).  Otherwise they are loaded only
        when they are required for calculations.

        The version and baseVersion for deformation calculations can be specified, otherwise
        the latest version is used by default.

        By default all components are loaded, but just one individual componenent can be
        selected.
        '''
        logging.info("Loading deformation model from %s",basedir)
        self._basedir = basedir
        if not os.path.isdir(basedir):
            raise ModelDefinitionError("Invalid deformation model base directory "+basedir)
        modfile = os.path.join(basedir,'model.csv')
        verfile = os.path.join(basedir,'version.csv')
        mtdfile = os.path.join(basedir,'metadata.csv')
        for f in (modfile, verfile, mtdfile ):
            if not os.path.isfile(f):
                raise ModelDefinitionError("File " + modfile + " is missing from deformation model")

        versions = {}
        curversion = None
        for ver in CsvFile('version',verfile,self.versionspec):
            if ver.version in versions:
                raise ModelDefinitionError('Version '+ver.version+' repeated in '+verfile)
            versions[ver.version] = ver
            if curversion == None or ver.version > curversion:
                curversion = ver.version
        self._versions = versions
        self._curversion = curversion

        metadata = {}
        for mtd in CsvFile('metadata',mtdfile,self.metadataspec):
            metadata[mtd.item] = mtd.value
        self._metadata = metadata

        for item in self.metadataitems:
            if item not in metadata:
                raise ModelDefinitionError('Metadata item '+item+' missing in '+mtdfile)

        mtdversion = str(metadata['version'])
        if mtdversion not in versions:
            raise ModelDefinitionError('Version '+mtdversion+' from metadata is not defined in version.csv')
        elif mtdversion != curversion:
            raise ModelDefinitionError('Version '+mtdversion+' from metadata is not most recent version in version.csv')


        self._name = str(metadata['model_name'])
        self._datumcode = str(metadata['datum_code'])
        self._datumname = str(metadata['datum_name'])
        self._ellipsoid = None
        try:
            self._datumsrid = int(metadata['datum_epsg_srid'])
        except:
            raise ModelDefinitionError("Invalid datum EPSG srid - must be an integer")
        try:
            self._datumepoch=Time.Parse(metadata['datum_epoch'])
        except:
            message = str(sys.exc_info()[1])
            raise ModelDefinitionError("Invalid datum epoch in "+mtdfile+": "+message)

        # List of model components, and hash of spatial files used to identify which have
        # already been loaded

        self._components = []
        self._spatial_models={}
        self._temporal_models={}
        self._cache=None

        cacheFile = os.path.join(self._basedir,'cache.h5')
        if clearCache and os.path.exists(cacheFile):
            os.remove(cacheFile)
        if useCache:
            self._cache = Cache(cacheFile)

        # Components to use.  Default is all.  Specific components can be selected
        # as "compononent+...+component", or "-component+component+....+component"

        componentList=[]
        useList=True
        if loadComponent:
            if loadComponent.startswith('-'):
                useList=False
                loadComponent=loadComponent[1:]
            componentList=loadComponent.split('+')

        for mdl in CsvFile('model',modfile,self.modelspec):
            component = mdl.component
            if componentList:
                if component in componentList:
                    if not useList:
                        continue
                elif useList:
                    continue

            if mdl.version_added not in versions:
                raise ModelDefinitionError("Model component "+mdl.component+" version_added "+mdl.version_added+" is not in version.csv")
            if mdl.version_revoked != '0' and mdl.version_revoked not in versions:
                raise ModelDefinitionError("Model component "+mdl.component+" version_revoked "+mdl.version_revoked+" is not in version.csv")
            compbase = os.path.join(basedir,component)
            if not os.path.isdir(compbase):
                raise ModelDefinitionError("Model component "+mdl.component+" directory is missing")
            compfile = os.path.join(compbase,'component.csv')
            if not os.path.isfile(compfile):
                raise ModelDefinitionError("Model component "+mdl.component+" component.csv file is missing")
            compname = os.path.join(component,'component.csv')

            filehashcheck = {}
            subcomponents = {}
            for compdef in CsvFile('component',compfile,self.componentspec):
                if compdef.version_added not in versions:
                    raise ModelDefinitionError("Component version_added " + compdef.version_added + " in " + compname + "is not in version.csv")
                if compdef.version_revoked != '0' and compdef.version_revoked not in versions:
                    raise ModelDefinitionError("Component version_revoked "+compdef.version_revoked+" in "+compname+" is not in version.csv")
                if compdef.displacement_type == 'none' and compdef.error_type == 'none':
                    raise ModelDefinitionError("Component in "+compname+" has displacement_type and error_type as none")

                spatial = SpatialComponent.get( self, self._spatial_models, component, compdef, load )
                temporal = TimeComponent.get( self._temporal_models, compdef )

                subcomponentid = compdef.subcomponent

                if subcomponentid > 0:
                    if subcomponentid in subcomponents:
                        subcomponents[subcomponentid].addComponent(spatial,compdef)
                        continue
                    else:
                        subcomp = SpatialComponentSet(component,spatial,compdef)
                        subcomponents[subcomponentid] = subcomp
                        spatial = subcomp

                self._components.append(
                    ModelComponent(component,mdl.description,compdef,spatial,temporal ))

        self._stcomponents = []
        self._version = ''
        self._baseVersion = ''
        self._versionName = ''
        self.setVersion( version, baseVersion )

    def __enter__( self ):
        return self

    def __exit__( self, exc_type, exc_value, traceback ):
        self.close()
        return False

    def close( self ):
        '''
        Release resources to avoid circular links that will prevent garbage collection.
        '''
        self._spatial_models = None
        self._temporal_models = None
        self._components = None
        self._stcomponents = None
        if self._cache:
            self._cache.close()
            self._cache = None

    def setVersion( self, version=None, baseVersion=None ):
        '''
        Reset the version used for calculations.  If the version is None, then the current version 
        is used.
        
        If baseVersion is specified then the calculation will be the difference between the 
        version and the baseVersion. 
        '''
        self._stcomponents = []
        if version==None:
            version = self._curversion
        else:
            version = str(version) 
            if version not in self._versions:
                raise ValueError("Requested version "+version+" of deformation model is not defined")
            
        if baseVersion:
            baseVersion = str(baseVersion)
            if baseVersion not in self._versions:
                raise ValueError("Requested base version "+baseVersion+" of deformation model is not defined")
        for c in self._components:
            factor = 0
            if c.appliesForVersion(version):
                factor = 1
            if baseVersion and c.appliesForVersion(baseVersion):
                factor -= 1
            if factor != 0:
                c.setFactor(factor)
                self._stcomponents.append(c)
            
        vername = version
        if baseVersion:
            vername = version + '-' + baseVersion
        self._version = version
        self._baseVersion = baseVersion
        self._versionName = vername
        self._date = None
        self._baseDate = None
        self._timeRangeError=None

    def getFileName( self, *parts ):
        return os.path.join(self._basedir,*parts)

    def _cacheMetadata( self, metadata, files ):
        metadataparts=[]
        for f in files:
            mtime=os.path.getmtime(self.getFileName(f))
            mtimestr=time.strftime("%Y%m%d%H%M%S",time.localtime(mtime))
            metadataparts.append(f.replace('\\','/'))
            metadataparts.append(mtimestr)
        if metadata:
            metadataparts.extend([str(x) for x in metadata])
        return ':'.join(metadataparts)

    def cacheData( self, file, metadata=None, files=None ):
        if not self._cache:
            return None
        files = files or [file]
        metadata=self._cacheMetadata(metadata,files)
        file = file.replace('\\','/')
        return self._cache.get(file,metadata)

    def setCacheData( self, data, file, metadata=None, files=None ):
        if not self._cache:
            return
        files = files or [file]
        metadata=self._cacheMetadata(metadata,files)
        file = file.replace('\\','/')
        self._cache.set(file,metadata,data)

    def setDate( self, date=None, baseDate=None ):
        '''
        Set the date when the deformation will be calculated, and optionally 
        a baseDate to calculate the difference in deformation between the date
        and baseDate.
        '''
        if date == None:
            date = datetime.datetime.now()
        if date != self._date or baseDate != self._baseDate:
            self._date = date
            self._baseDate = baseDate
            self._timeRangeError=None
            try:
                for comp in self._stcomponents:
                    comp.setDate( date, baseDate )
            except OutOfRangeError:
                self._timeRangeError = sys.exc_info()[1]

    def calcDeformation( self, x, y, date=None, baseDate=None ):
        '''
        Calculate the deformation at a specified location.  The date and 
        baseDate can be set at the same time if required, otherwise the 
        values set with setDate will be used.
        '''
        if not self._date or date != None or baseDate != None:
            self.setDate( date, baseDate )
        if self._timeRangeError:
            raise OutOfRangeError( self._timeRangeError )

        result = [0.0,0.0,0.0,0.0,0.0]
        for comp in self._stcomponents:
            compvalue = comp.calcDeformation( x, y )
            for i in range(5):
                result[i] += compvalue[i]
        result[3]=math.sqrt(abs(result[3]))
        result[4]=math.sqrt(abs(result[4]))
        return result

    def ellipsoid( self ):
        if self._ellipsoid is None:
            from ..Geodetic import ellipsoid
            a=float(self._metadata['ellipsoid_a'])
            rf=float(self._metadata['ellipsoid_rf'])
            self._ellipsoid=ellipsoid.ellipsoid(a,rf)
        return self._ellipsoid

    def applyTo( self, lon, lat=None, hgt=None, date=None, baseDate=None, subtract=False ):
        '''
        Applies the deformation to longitude/latitude coordinates.
        
        Input can be one of 
           lon,lat
           lon,lat,hgt
           [lon,lat],
           [lon,lat,hgt],
           [[lon,lat],[lon,lat]...]
           [[lon,lat,hgt],[lon,lat,hgt]...]

        For the first four cases returns a single latitude/longitude/height
        For the other cases returns an array of [lon,lat,hgt]

        The deformation is added to the coordinates unless subtract is True,
        in which case it is removed.
        '''
        import numpy as np
        ell=self.ellipsoid()
        if lat is None:
            crds=lon
            if not isinstance(crds,np.ndarray):
                crds=np.array(crds)
            single=len(crds.shape)==1
            if single:
                crds=crds.reshape((1,crds.size))
        else:
            single=True
            crds=[[lon,lat,hgt or 0]]

        results=[]
        factor=-1 if subtract else 1
        for crd in crds:
            ln,lt=crd[:2]
            ht=crd[2] if len(crd) > 2 else 0
            deun=self.calcDeformation(ln,lt,date,baseDate)[:3]
            dedln,dndlt=ell.metres_per_degree(ln,lt)
            results.append([ln+factor*deun[0]/dedln,lt+factor*deun[1]/dndlt,ht+factor*deun[2]])
        return results[0] if single else np.array(results)

    def name( self ):
        '''
        Return the name of this model
        '''
        return self._name

    def versionName( self ):
        '''
        Returns the name of the version of the model currently set
        '''
        return self._versionName

    def version(self):
        '''
        Returns the current version number
        '''
        return self._version

    def baseVersion(self):
        '''
        Returns the current version base number if calculating a difference, or None
        if not
        '''
        return self._baseVersion

    def currentVersion( self ):
        '''
        Returns the current version of the model (ie the latest version)
        '''
        return self._curversion

    def versions(self):
        '''
        Returns a list of versions available
        '''
        return sorted(self._versions.keys())

    def versionInfo( self, version ):
        '''
        Returns information from the versions.csv file for a specific version of the model
        '''
        return self._versions[version]

    def datumName( self ):
        '''
        Returns the name of the datum defined in the model metadata
        '''
        return self._datumname
    def datumCode( self ):
        '''
        Returns the datum code defined in the model metadata
        '''
        return self._datumcode

    def datumEpoch( self ):
        '''
        Returns the datum epoch defined in the model metadata
        '''
        return self._datumepoch

    def datumEpsgSrid( self ):
        '''
        Returns the EPSG srid (spatial reference id) of the datum
        latitude/longitude coordinate system
        '''
        return self._datumsrid

    def components( self, allversions=False ):
        compkey=lambda c: (0 if c.component == 'ndm' else 1,c.versionAdded,c.component)
        for c in sorted(self._components,key=compkey):
            if allversions or c.appliesForVersion(self.version()):
                yield c

    def description( self, allversions=False, components=True ):
        '''
        Return a description of the model
        '''
        import StringIO
        outs = StringIO.StringIO()
        mtd = self._metadata
        outs.write("Deformation model: "+mtd['model_name']+"\n")
        outs.write("Datum: "+mtd['datum_name']+" (reference epoch "+mtd['datum_epoch']+")\n")
        outs.write("Version: "+self.version()+"\n")
        outs.write("\n")
        outs.write(mtd['description']+"\n")

        if allversions:
            outs.write("\nVersions available:\n")
            for version in self.versions():
                v = self.versionInfo(version)
                outs.write("    "+v.version+
                          " released "+v.release_date.strftime('%d-%b-%Y')+
                          ": "+v.reason+"\n")

        if components:
            compcount = {}
            for c in self.components(allversions): 
                compcount[c.component] = compcount[c.component]+1 if c.component in compcount else 1

            outs.write("\nComponents:\n")
            lastcomponent=None
            for c in self.components( allversions ):
                if c.component != lastcomponent:
                    description = c.compdesc.strip().replace("\n","\n        ")
                    outs.write("\n    Component: "+c.component+": " + description +"\n")
                    lastcomponent=c.component
                prefix="    "
                if compcount[c.component] > 1:
                    description = c.description.strip().replace("\n","\n            ")
                    outs.write("        Sub-component: "+description+"\n")
                    prefix="        "
                if allversions:
                    outs.write(prefix+"    Version added: "+c.versionAdded)
                    if c.versionRevoked != '0':
                        outs.write(" revoked: "+c.versionRevoked)
                    outs.write("\n")
                outs.write(prefix+"    Time model: "+str(c.timeComponent)+"\n")
                description = str(c.spatialComponent).replace("\n","\n        "+prefix)
                outs.write(prefix+"    Spatial model: "+description+"\n")

        description = outs.getvalue()
        outs.close()
        return description

    def __str__( self ):
        return self.description(allversions=True)
