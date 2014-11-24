#!/usr/bin/python
import sys
import datetime

__version__='1.1'

# Setup the ITRF transformation

from LINZ.DeformationModel.Time import Time
from LINZ.DeformationModel import Model
from LINZ.DeformationModel.Error import ModelDefinitionError, OutOfRangeError, UndefinedValueError

from LINZ.Geodetic import ellipsoid
from LINZ.Geodetic import ITRF_transformation

class Transformation( object ):
    '''
    Class for transforming between ITRF's and NZGD2000. The usage is

        tfm=ITRF_NZGD2000.Transformation('ITRF2008',toNZGD2000=True)
        nzlon,nzlat,nzhgt=tfm(itrflon,itrflat,itrfhgt,date)

        itrf=tfm.itrf
        model=tfm.model
        version=tfm.version

    '''

    itrf=None
    toNZGD2000=None
    model=None
    modeldir=None
    version=None
    transform=None

    def __init__( self, 
                 itrf='ITRF2008', 
                 toNZGD2000=True, 
                 modeldir=None, 
                 version=None,
                 usecache=True,
                 clearcache=False):

        '''
        Set up an ITRF to NZGD2000 transformation or reverse transformation

        Arguments:

            itrf         The itrf to transform from/to eg 'ITRF2008'
            toNZGD2000   If false then the transformation is NZGD2000->ITRF
            modeldir     The base directory of the deformation model
            version      The version of the deformation model to use
            usecache     If true then use the binary cached model
            clearcache   If true then delete the cached model

        '''

        if not modeldir:
            from os.path import dirname, abspath, join
            modeldir = join(dirname(dirname(abspath(__file__))),'model')

        model = Model.Model(modeldir,useCache=usecache,clearCache=clearcache )
        if version == None:
            version = model.currentVersion()
        model.setVersion( version )
            
        try:
            itrf=itrf.upper()
            itrf_src=ITRF_transformation.transformation(from_itrf='ITRF96',to_itrf=itrf)
            if toNZGD2000:
                itrf_src = itrf_src.reversed()
            itrf_tfm=itrf_src.transformLonLat
        except:
            raise RuntimeError( "Invalid ITRF "+itrf )

        if toNZGD2000:
            def transform( lon, lat, hgt, date ):
                if type(date) != Time:
                    date=Time(date)
                llh=itrf_tfm( lon, lat, hgt, date=date.asYear() )
                llh=model.applyTo( llh, date=date, subtract=True )
                return llh
        else:
            def transform( lon, lat, hgt, date ):
                if type(date) != Time:
                    date=Time(date)
                llh=model.applyTo( lon, lat, hgt, date=date.asYear() )
                llh=itrf_tfm( llh, date=date.asYear() )
                return llh

        self.itrf=itrf
        self.toNZGD2000=toNZGD2000
        self.model=model
        self.modeldir=modeldir
        self.version=version
        self.transform=transform

    def __call__( self, lon, lat, hgt, date ):
        return self.transform( lon, lat, hgt, date )

    def __del__( self ):
        self.close()

    def close(self):
        if self.model:
            self.model.close()

