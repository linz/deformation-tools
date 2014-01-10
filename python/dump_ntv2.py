
import argparse
import os
import os.path
import re
from collections import namedtuple
import subprocess
from math import *
from datetime import date
import struct


parser=argparse.ArgumentParser('Dump NTv2 to ASCI - input ntv2 file, ascii file')
parser.add_argument('-A','--australian',action='store_true',help='Build australian format binary')
parser.add_argument('-B','--big-endian',action='store_true',help='File is big endian (default little endian)')
parser.add_argument('-d','--data',action='store_true',help='Include data values in output')
parser.add_argument('-c','--csv',help='CSV output file for data')
parser.add_argument('ntv2_filename',help='Name of the input binary NTv2 grid file')
parser.add_argument('ascii_filename',help='Name of the output ASCII grid file')

args=parser.parse_args()

australian=args.australian
bigendian=args.big_endian
ntfile=args.ntv2_filename
asciifile=args.ascii_filename
showdata=args.data

print "Dumping NTv2 file",ntfile

endian = '>' if bigendian else '<'
sformat=endian+'8s8s'
iformat=endian+'8si' if australian else endian+'8si4x'
dformat=endian+'8sd'
gsformat=endian+'ffff'
iflen=12 if australian else 16

if args.csv:
    csvfile=open(args.csv,'w')
    csvfile.write("grid,ptno,lon,lat,dlon,dlat,parent\n")

with open(ntfile,'rb') as nt, open(asciifile,'w') as outf:

    buffer=nt.read(iflen)
    item,head_len=struct.unpack(iformat,buffer)
    assert item=='NUM_OREC'

    buffer=nt.read(iflen)
    item,grid_head_len=struct.unpack(iformat,buffer)
    assert item=='NUM_SREC'

    assert head_len==11
    assert grid_head_len==11

    buffer=nt.read(iflen)
    item,nfile=struct.unpack(iformat,buffer)
    assert item=='NUM_FILE'

    buffer=nt.read(16)
    item,gs_type=struct.unpack(sformat,buffer)
    assert item=='GS_TYPE '

    buffer=nt.read(16)
    item,version=struct.unpack(sformat,buffer)
    assert item=='VERSION '

    buffer=nt.read(16)
    item,from_system=struct.unpack(sformat,buffer)
    assert item=='SYSTEM_F'

    buffer=nt.read(16)
    item,to_system=struct.unpack(sformat,buffer)
    assert item=='SYSTEM_T'

    buffer=nt.read(16)
    item,from_a=struct.unpack(dformat,buffer)
    assert item=='MAJOR_F '

    buffer=nt.read(16)
    item,from_b=struct.unpack(dformat,buffer)
    assert item=='MINOR_F '

    buffer=nt.read(16)
    item,to_a=struct.unpack(dformat,buffer)
    assert item=='MAJOR_T '

    buffer=nt.read(16)
    item,to_b=struct.unpack(dformat,buffer)
    assert item=='MINOR_T '

    outf.write("NUM_FILE: {0}\n".format(nfile))
    outf.write("GS_TYPE: {0}\n".format(gs_type))
    outf.write("VERSION: {0}\n".format(version))
    outf.write("SYSTEM_F: {0}\n".format(from_system))
    outf.write("SYSTEM_T: {0}\n".format(to_system))
    outf.write("MAJOR_F: {0}\n".format(from_a))
    outf.write("MINOR_F: {0}\n".format(from_b))
    outf.write("MAJOR_T: {0}\n".format(to_a))
    outf.write("MINOR_T: {0}\n".format(to_b))

    scale=1.0/3600.0 if gs_type.strip() == 'SECONDS' else 1.0;
               

    ntgridname={}
    for i in range(nfile):

            buffer=nt.read(16)
            item,name=struct.unpack(sformat,buffer)
            assert item=='SUB_NAME'

            buffer=nt.read(16)
            item,pname=struct.unpack(sformat,buffer)
            assert item=='PARENT  '

            buffer=nt.read(16)
            item,created=struct.unpack(sformat,buffer)
            assert item=='CREATED '

            buffer=nt.read(16)
            item,updated=struct.unpack(sformat,buffer)
            assert item=='UPDATED '

            buffer=nt.read(16)
            item,lat0=struct.unpack(dformat,buffer)
            assert item=='S_LAT   '

            buffer=nt.read(16)
            item,lat1=struct.unpack(dformat,buffer)
            assert item=='N_LAT   '

            buffer=nt.read(16)
            item,lon1=struct.unpack(dformat,buffer)
            assert item=='E_LONG  '

            buffer=nt.read(16)
            item,lon0=struct.unpack(dformat,buffer)
            assert item=='W_LONG  '

            buffer=nt.read(16)
            item,latinc=struct.unpack(dformat,buffer)
            assert item=='LAT_INC '

            buffer=nt.read(16)
            item,loninc=struct.unpack(dformat,buffer)
            assert item=='LONG_INC'

            buffer=nt.read(iflen)
            item,ngridpt=struct.unpack(iformat,buffer)
            assert item=='GS_COUNT'

            lon0=-lon0;
            lon1=-lon1;

            outf.write("\nGRID {0}:\n".format(i+1))
            outf.write("  SUB_NAME: {0}\n".format(name))
            outf.write("  PARENT: {0}\n".format(pname))
            outf.write("  CREATED: {0}\n".format(created))
            outf.write("  UPDATED: {0}\n".format(updated))
            outf.write("  S_LAT: {0}\n".format(lat0*scale))
            outf.write("  N_LAT: {0}\n".format(lat1*scale))
            outf.write("  E_LONG: {0}\n".format(lon0*scale))
            outf.write("  W_LONG: {0}\n".format(lon1*scale))
            outf.write("  LAT_INC: {0}\n".format(latinc*scale))
            outf.write("  LONG_INC: {0}\n".format(loninc*scale))
            outf.write("  GS_COUNT: {0}\n".format(ngridpt))

            nln=int(abs(lon0-lon1)/loninc+0.5)+1;
            nlt=int(abs(lat0-lat1)/latinc+0.5)+1;

            assert nln*nlt==ngridpt, "nln*nlt={0}*{1}={2} != gs_count={3}".format(nln,nlt,nln*nlt,ngridpt)

            ipt=0
            for ilat in range(nlt):
                lat=(lat0+ilat*latinc)*scale
                for ilon in range(nln):
                    ipt += 1
                    lon=(lon1-ilon*loninc)*scale

                    buffer=nt.read(16)
                    (dlt,dln,errlt,errln)=struct.unpack(gsformat,buffer)
                    dln = -dln;
                    
                    if csvfile:
                        csvfile.write('"{0}",{1},{2},{3},{4},{5},"{6}"\n'.format(
                            name.strip(),ipt,lon,lat,dln,dlt,pname.strip()))
                    if showdata:
                        outf.write("  Point {0}: {1} {2} {3} {4}\n".format(ipt+1,dlt,dln,errlt,errln))

if csvfile:
    csvfile.close()
