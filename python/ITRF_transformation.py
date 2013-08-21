
import math
import numpy as np

secstorads = math.radians(1.0/3600.0)

class ITRF_transformation( object ):

    def __init__( self, rffrom, rfto, params, rates=None, refdate=None, source=None ):
        self.rffrom=rffrom
        self.rfto=rfto
        self.params=None if params is None else list(params)
        self.rates=None if rates is None else list(rates)
        self.refdate=refdate
        self.source=source
        self._tf=None

    def setRates( self, rates ):
        self.rates=rates

    def __str__(self):
        return ("Transformation from "+self.rffrom+" to "+self.rfto+"\n"+
                (
                " Reference date {0:.1f}\n".format(self.refdate)
                    if self.rates and self.refdate else '') +
                "   Translations {0:.2f}  {1:.2f}  {2:.2f} mm\n".format(*self.params[0:3])+
                (
                "          rates {0:.4f}  {1:.4f}  {2:.4f} mm/yr\n".format(*self.rates[0:3])
                    if self.rates and self.refdate else '') +
                "      Rotations {0:.2f}  {1:.2f}  {2:.2f} mas\n".format(*self.params[4:7])+
                (
                "          rates {0:.4f}  {1:.4f}  {2:.4f} mas/yr\n".format(*self.rates[4:7])
                    if self.rates and self.refdate else '') +
                "          Scale {0:.2f}  ppb\n".format(self.params[3])+
                (
                "          rates {0:.4f} ppb/yr\n".format(self.rates[3])
                    if self.rates and self.refdate else ''))

    def reversed( self ):
        return ITRF_transformation(
            self.rfto,
            self.rffrom,
            [-p for p in self.params],
            None if self.rates is None else [-r for r in self.rates],
            self.refdate,
            self.source
            )

    def atDate( self, date ):
        p=list(self.params)
        refdate=None
        if self.rates and self.refdate:
            diff = date - self.refdate
            for i,r in enumerate(self.rates):
                p[i]=self.params[i] + r*diff
            refdate = date
        return ITRF_transformation(
            self.rffrom,
            self.rfto,
            p,
            self.rates,
            refdate,
            self.source )

    def add( self, other ):
        if self.rffrom==other.rfto:
            rffrom=other.rffrom
            rfto=self.rfto
        elif self.rfto==other.rffrom:
            rffrom=self.rffrom
            rfto=other.rfto
        else:
            raise RuntimeError("Cannot join incompatible transformations (must have common start/end reference frame")

        refdate=self.refdate if self.refdate is not None else other.refdate
        if refdate and refdate != other.refdate:
            other=other.atDate(refdate)

        return ITRF_transformation(
            rffrom,
            rfto,
            [p1+p2 for p1,p2 in zip(self.params,other.params)],
            (self.rates if other.rates is None 
             else other.rates if self.rates is None 
             else [p1+p2 for p1,p2 in zip(self.rates,other.rates)]),
            refdate,
            source if self.source == other.source else None
            )

    def subtract( self, other ):
        return self.add( other.reversed())

    def transFunc( self, date ):
        params=self.params
        if self.rates and self.refdate and date:
            diff = date - self.refdate
            params=[p+r*diff for p,r in zip(params,self.rates)]
        txyz=np.array([[params[0],params[1],params[2]]])*0.001
        scale=params[3]*1.0e-9
        rotscale=secstorads*0.001
        rx=params[4]*rotscale
        ry=params[5]*rotscale
        rz=params[6]*rotscale
        rxyz=np.transpose(np.array([[scale,-rz,ry],[rz,scale,-rx],[-ry,rx,scale]]))
        def tf( coords ):
            if not isinstance(coords,np.ndarray):
                coords=np.array(coords)
            single=len(coords.shape)==1
            if single:
                coords=coords.reshape((1,coords.size))
            coords=coords+txyz+coords.dot(rxyz)
            if single:
                coords=coords.reshape((coords.size))
            return coords
        return tf

    def transform( self, xyz, date=None ):
        if self._tf is None or date != self._tfdate:
            self._tf=self.transFunc(date)
            self._tfdate=date
        return self._tf(xyz)


scalefactors = (0.001, 0.001, 0.001, 1.0e-9, secstorads*0.001, secstorads*0.001, secstorads*0.001 )

# Data from http://itrf.ensg.ign.fr/doc_ITRF/Transfo-ITRF2008_ITRFs.txt

# Note : These parameters are derived from those already published in the IERS
# Technical Notes and Annual Reports. The transformation parameters should be
# used with the standard model (1) given below and are valid at the indicated
# epoch.
# 
# 
# : XS :    : X :   : Tx :   :  D   -Rz   Ry : : X :
# :    :    :   :   :    :   :               : :   :
# : YS :  = : Y : + : Ty : + :  Rz   D   -Rx : : Y :                       (1)
# :    :    :   :   :    :   :               : :   :
# : ZS :    : Z :   : Tz :   : -Ry   Rx   D  : : Z :
# 
# 
# Where X,Y,Z are the coordinates in ITRF2008 and XS,YS,ZS are the coordinates in
# the other frames.

iers_base='ITRF2008'

iers_data='''
  ITRF2005       -2.0     -0.9     -4.7      0.94      0.00      0.00      0.00    2000.0
       rates      0.3      0.0      0.0      0.00      0.00      0.00      0.00
  ITRF2000       -1.9     -1.7    -10.5      1.34      0.00      0.00      0.00    2000.0
       rates      0.1      0.1     -1.8      0.08      0.00      0.00      0.00
  ITRF97          4.8      2.6    -33.2      2.92      0.00      0.00      0.06    2000.0
       rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
  ITRF96          4.8      2.6    -33.2      2.92      0.00      0.00      0.06    2000.0
       rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
  ITRF94          4.8      2.6    -33.2      2.92      0.00      0.00      0.06    2000.0
       rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
  ITRF93        -24.0      2.4    -38.6      3.41     -1.71     -1.48     -0.30    2000.0
       rates     -2.8     -0.1     -2.4      0.09     -0.11     -0.19      0.07
  ITRF92         12.8      4.6    -41.2      2.21      0.00      0.00      0.06    2000.0
       rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
  ITRF91         24.8     18.6    -47.2      3.61      0.00      0.00      0.06    2000.0
       rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
  ITRF90         22.8     14.6    -63.2      3.91      0.00      0.00      0.06    2000.0
       rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
  ITRF89         27.8     38.6   -101.2      7.31      0.00      0.00      0.06    2000.0
       rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
  ITRF88         22.8      2.6   -125.2     10.41      0.10      0.00      0.06    2000.0
       rates      0.1     -0.5     -3.2      0.09      0.00      0.00      0.02
       '''
# From 
# Transforming Positions and Velocities between the International Terrestrial Reference Frame of 2000 # and North American Datum of 1983
# Tomas Soler, M.ASCE, and Richard A. Snay
# JOURNAL OF SURVEYING ENGINEERING (c) ASCE / MAY 2004 / P49
#
# Referenced to 
# Springer, T. A., Kouba, J., and Mireault, Y. 2000. "1999 analysis coor-
# dinator report." 1999 Tech. Rep., International GPS Service for Geo-
# dynamics, Jet Propulsion Laboratory, Pasadena, Calif., 15-55.

igs_base='ITRF97'
igs_data='''
   ITRF96 -2.07 -0.21 9.95 -0.93496 +0.12467 -0.22355 -0.06065 1997.0
     rates 0.69 -0.10 1.86 -0.19201 +0.01347 -0.01514 +0.00027 
'''

def read_trans( data, base, source, transformations ):
    for l in data.split('\n'):
        parts=l.split()
        if len(parts) < 7:
            continue
        code=parts[0]
        params=[float(x) for x in parts[1:8]]
        if code.startswith('ITRF'):
            trans=ITRF_transformation(base,code,params,refdate=float(parts[8]),source=source)
            transkey=base+'-'+code
            transformations[transkey]=trans
        elif code == 'rates':
            trans.setRates(params)


iers_transformations={}
igs_transformations={}

read_trans(iers_data, iers_base, 'IERS', iers_transformations )
read_trans(igs_data, igs_base, 'IGS', igs_transformations )

itrf2008_nzgd2000=iers_transformations['ITRF2008-ITRF97'].add(igs_transformations['ITRF97-ITRF96']).atDate(2000)

if __name__=="__main__":

    print "IGS ITRF97 to ITRF96"
    print igs_transformations['ITRF97-ITRF96']

    print "IGS ITRF97 to ITRF96 at 2000.0"
    print igs_transformations['ITRF97-ITRF96'].atDate(2000.0)

    print "ITRF2008 to ITRF96"
    print iers_transformations['ITRF2008-ITRF96']

    print "ITRF2008 to NZGD2000"
    print itrf2008_nzgd2000

    crd=[-5115333.2235,   477886.9008, -3767147.4737]
    crd2=[-5153430.4382, 513057.5295, -3710656.1942]
    crds=[crd,crd2]
    result=itrf2008_nzgd2000.transform(crd)
    print result
    print itrf2008_nzgd2000.params
    print itrf2008_nzgd2000.transform(crds)
    print itrf2008_nzgd2000.transform(crds,date=2013.5)


