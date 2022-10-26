
import os.path
import csv
import inspect
import re
from datetime import datetime

class CsvResults( object ):
    ''' 
    Store results in CSV file

    
    Results are stored indexed by on or more key fields.  
    '''

    def __init__( self, filename=None, key=None, fields=None ):
        '''
        Create the CsvResults object. 

        filename is the name of the file, the default is caller.csv
        key is a list of fields used to index the results.  Only one value
        can be stored for each key.  Supplied as a list or space separated
        string.
        fields is a list of fields in the file (any fields in a call to
        save will be added to this.  Also a default field _date will always
        be added.
        '''
        if not filename:
            filename=inspect.getframeinfo(inspect.stack()[1][0])[0]
            filename=re.sub(r'\.pyc?','',filename,re.I)
            filename=filename+'.csv'
        self._filename=filename
        if not key:
            key=[]
        self._key = key if type(key)==list else key.split()
        if not fields:
            fields=[]
        self._fields = fields if type(fields)==list else fields.split()
        for k in self._key:
            if k not in self._fields:
                self._fields.append(k)
        self._data={}
        self._load()

    def _datakey(self,data):
        return ':'.join([str(data.get(k)) for k in self._key]) if self._key else len(self._data)+1

    def _load(self):
        if os.path.exists( self._filename ):
            with open(self._filename) as f:
                c=csv.DictReader(f)
                for f in c.fieldnames: 
                    if f not in self._fields:
                        self._fields.append(f)
                data={}
                for r in c:
                    key=self._datakey(r)
                    data[key]=r
                self._data=data

    def _save(self):
        with open(self._filename,'w') as f:
            self._fields.remove('_date')
            self._fields.append('_date')
            c=csv.DictWriter(f,self._fields)
            c.writeheader()
            for k in sorted(self._data.keys()):
                c.writerow(self._data[k])

    def save( self, **data ):
        '''
        Save data - data is supplied as a list of key/value pairs
        eg

        CsvResults results
        results.save(id=42,name='Pooh',preference='Honey')
        '''
        if not isinstance(data,dict):
            raise ValueError('CsvResults can only save dict objects')
        data=dict(**data)
        data['_date'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        for k in list(data.keys()):
            if k not in self._fields:
                self._fields.append(k)
        self._data[self._datakey(data)]=data
        self._save()

    def recall( self, **key ):
        if not self._key:
            raise ValueError("Cannot use CsvResult.recall - no key defined for store")
        key = self._datakey(key)
        return self._data[key]









        
