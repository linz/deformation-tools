#!/usr/bin/python

from LINZ.Geodetic.Ellipsoid import GRS80
import sys
import os.path
import argparse
import math

parser=argparse.ArgumentParser(description='Basic comparison of lat/lon coordinate files')
parser.add_argument('infile1',help="First input file")
parser.add_argument('infile2',help="Second input file")
parser.add_argument('outfile',help="Output file")
parser.add_argument('-s','--summary-id',help="Output summary difference with specified id")
args=parser.parse_args()

clists=[]
for infile in (args.infile1, args.infile2):
    if not os.path.exists(infile):
        print("Input file {0} missing".format(infile))
    coords=[]
    with open(infile) as ifh:
        for l in ifh:
            p=l.split()
            if len(p) >= 2:
                coords.append([float(x) for x in p[:2]])
    clists.append(coords)

maxdiff=0
npt=0
with open(args.outfile,'w') as of:
    for c1,c2 in zip(clists[0],clists[1]):
        dedln,dndlt=GRS80.metres_per_degree(*c1)
        de=(c2[0]-c1[0])*dedln
        dn=(c2[1]-c1[1])*dndlt
        diff=math.sqrt(de*de+dn*dn)
        if diff > maxdiff:
            maxdiff=diff
        npt += 1
        of.write("{0:.8f} {1:.8f} {2:.8f} {3:.8f} {4:7.3f} {5:7.3f}\n"
                 .format(c1[0],c1[1],c2[0],c2[1],de,dn))

if args.summary_id:
    print("{0} {1} test points max diff {2:.3f}m".format(args.summary_id,npt,maxdiff))




