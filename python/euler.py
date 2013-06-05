
import numpy as np
import math
import ellipsoid

class euler_rotation( object ):

    @staticmethod
    def _cossin( angle ):
        angle=np.radians(angle)
        return np.cos(angle),np.sin(angle)

    def __init__( self, ellipsoid, lon, lat, rate ):
        '''
        Initiallize a euler rotation definition,
        Longitude and latitude in degrees
        Rotation in radians / 10**7 years
        '''
        self._ellipsoid = ellipsoid
        self._eaxes = euler_rotation.axes( lon, lat )
        self._rate = rate * 1E-7;

    @staticmethod
    def axes( lon, lat ):
        cln,sln = euler_rotation._cossin(lon)
        clt,slt = euler_rotation._cossin(lat)
        vz=np.array([clt*cln,clt*sln,slt])
        if abs(slt) > 0.5:
            rv=np.array([1,0,0])
        else:
            rv=np.array([0,0,1])
        vx = np.cross(vz,rv)
        vx /= np.linalg.norm(vx)
        vy = np.cross(vz,vx)
        return np.vstack((vx,vy,vz))

    def velocity( self, lon, lat ):
        xyz=self._ellipsoid.xyz(lon,lat)
        x=self._eaxes[0].dot(xyz) * self._rate
        y=self._eaxes[1].dot(xyz) * self._rate
        dxyz=self._eaxes[1]*x - self._eaxes[0]*y
        dxyz=self._ellipsoid.enu_axes(lon,lat).dot(dxyz)
        return (dxyz[0],dxyz[1])

if __name__ == '__main__':
    import argparse
    import pointsource

    parser = argparse.ArgumentParser(description='Calculate displacement from Euler rotation')
    parser.add_argument('euler',type=float,nargs=3,help='Euler rotation - lon lat rate.\nlon and lat in degrees, rate in radians/10**7 years')
    parser.add_argument('point_source',help='Input points, csv file with lon, lat columns, or grid:min_lon:min_lat:max_lon:max_lat:nlon:nalat')
    parser.add_argument('output_file',help='Output file')

    args = parser.parse_args()

    rotation = euler_rotation( ellipsoid.grs80, *args.euler )

    output=open(args.output_file,"w")
    output.write('lon,lat,de,dn\n')

    for crd in pointsource.points(args.point_source):
        v=rotation.velocity(*crd)
        output.write("{0:.6f},{1:.6f},{2:.6f},{3:.6f}\n".format(crd[0],crd[1],v[0],v[1]))

