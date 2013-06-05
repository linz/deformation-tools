import numpy as np
import math

class ellipsoid( object ):

    @staticmethod
    def _cossin( angle ):
        angle=np.radians(angle)
        return np.cos(angle),np.sin(angle)

    @staticmethod
    def enu_axes( lon, lat ):
        cln,sln = ellipsoid._cossin(lon)
        clt,slt = ellipsoid._cossin(lat)
        ve=np.array([-sln,cln,0])
        vn=np.array([-cln*slt,-sln*slt,clt])
        vu=np.array([clt*cln,clt*sln,slt])
        return np.vstack((ve,vn,vu))

    def __init__( self, a, rf ):
        '''
        Initiallize an ellipsoid based on semi major axis and inverse flattening
        '''
        self.a=float(a)
        self.rf=float(rf)
        self.b=a-a/rf if rf else a
        self.a2=a*a
        self.b2=self.b*self.b
        self.a2b2=self.a2-self.b2

    def xyz( self, lon, lat ):
        cln,sln = ellipsoid._cossin(lon)
        clt,slt = ellipsoid._cossin(lat)
        bsac=np.hypot(self.b*slt,self.a*clt)
        p = self.a2*clt/bsac
        xyz=[p*cln,p*sln,self.b2*slt/bsac]
        return xyz

grs80 = ellipsoid(6378160.0,298.25)
