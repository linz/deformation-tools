import numpy as np
import math

nuvel_1a=dict(
    PCFC=[-0.0907,0.2902,-0.5976],
    AFRC=[0.0532,-0.1856,0.2348],
    ANTA=[-0.0494,-0.1018,0.2218],
    ARAB=[0.4003,-0.0311,0.4049],
    AUST=[0.4695,0.3072,0.3762],
    CARB=[-0.0109,-0.2027,0.0945],
    COCO=[-0.6249,-1.2944,0.6544],
    EURA=[-0.0590,-0.1434,0.1887],
    NAZC=[-0.0921,-0.5138,0.5756],
    NOAM=[0.0152,-0.2155,-0.0094],
    SOAM=[-0.0624,-0.0906,-0.0523],
    JUFU=[0.2995,0.4805,-0.2936],
    INDI=[0.3995,0.0026,0.4066],
    PHIL=[0.5913,-0.4412,-0.5976],
    )

def plates():
    '''
    Return a list of plates
    '''
    return sorted(nuvel_1a.keys())

def euler_rate( plate ):
    '''
    Return the rotation rate as an euler pole (lon,lat in degrees)
    and rotation rate (radians/10**7 years)
    '''
    rvec = nuvel_1a[plate]
    # Rate expressed as 
    rate = math.radians(np.linalg.norm(rvec)*10)
    lon = math.degrees(math.atan2(rvec[1],rvec[0]))
    lat = math.degrees(math.atan2(rvec[2],math.hypot(rvec[0],rvec[1])))
    return lon,lat,rate

if __name__ == '__main__':
    import sys
    import ellipsoid
    import euler
    import argparse
    import pointsource

    parser = argparse.ArgumentParser(description='Calculate displacement on NNR-Nuvel-1A')
    parser.add_argument('plate',help='Tectonic plate')
    parser.add_argument('point_source',help='Input points, csv file with lon, lat columns, or grid:min_lon:min_lat:max_lon:max_lat:nlon:nalat')
    parser.add_argument('output_file',help='Output file')

    args = parser.parse_args()

    plate = args.plate.upper()
    if not plate in nuvel_1a:
        print 'Plate '+plate+' not valid - must be one of '+', '.join(sorted(nuvel_1a.keys()))
        sys.exit()

    rotation = euler.euler_rotation( ellipsoid.grs80, *euler_rate(plate) )

    if args.output_file == '-':
        output = sys.stdout
    else:
        output=open(args.output_file,"w")
    output.write('lon,lat,de,dn\n')

    for crd in pointsource.points(args.point_source):
        v=rotation.velocity(*crd)
        output.write("{0:.6f},{1:.6f},{2:.6f},{3:.6f}\n".format(crd[0],crd[1],v[0],v[1]))

