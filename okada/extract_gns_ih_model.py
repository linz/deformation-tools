#!/usr/bin/python
import sys
import os
import os.path
import re
from math import radians, cos, sin, tan
from collections import namedtuple

syntax='''

extract_gns_model:  Generates a model file suitable for the calc_okada program
from the GNS text file format (Ian Hamling) and metadata file.  The text file 
is a simple space delimited file with columns defined by the metadata file.
The metadata file defines the header information to be used with the model, 
and the columns in the data file.

Syntax:
python extract_gns_model.py metadata_file [root_filename]

metadata_file is the metadata file for the GNS model
root_file_name is the output file base name

The script will be generated
   root_file_name.model      Containing the model definition

Example metadata file:
Event: Mw 6.6 Cook Strait earthquake, 21 July 2013
Source model: Geodetic source model, based on GPS; elastic half-space assumption
Version: 24 January 2014
Author: Ian Hamling, GNS Science
Source file: Cook_strait_slip_utm60.txt
Projection: UTM60
Columns: east_surface_m north_surface_m strike_deg dip_deg rake_deg slip_m depth_top_km depth_bottom_km
Additional data:  length_km=1.0

The "Source file", and "columns" entries are required for this program.
the Projection entry is required by the calc_okada program.

The current implementation requires all the above columns or additional data
to be defined.

'''

def extract_model( meta_file, outfile=None ):
    sourcefile=''
    metadata=[]
    columns=[]
    extracols=[]
    extradata=[]
    with open(meta_file) as mf:
        for l in mf:
            m=re.match(r'([^\:]+)\:\s*(.*?)\s*$',l)
            if not m:
                continue
            command=m.group(1)
            data=m.group(2)

            if command.lower() == 'source file':
                sourcefile=data
            elif command.lower() == 'columns':
                columns.extend(data.split())
                continue
            elif command.lower() == 'additional data':
                for cv in data.split():
                    if '=' in cv:
                        c,v=cv.split('=')
                        extracols.append(c)
                        extradata.append(float(v))
                continue
            metadata.append(l)

    if not sourcefile:
        raise RuntimeError("\"Source file\" not defined in metadata")
    ncols = len(columns)
    if not columns:
        raise RuntimeError("\"Columns\" not defined in metadata")
    columns.extend(extracols)
    FaultData=namedtuple('FaultData',columns)
    faults=[]

    datadir=os.path.dirname(meta_file)
    sfn=sourcefile
    if datadir:
        sfn=os.path.join(datadir,sourcefile)

    with open(sfn) as sf:
        for l in sf:
            if l.strip() == '':
                continue
            data=[float(x) for x in l.split()]
            data.extend(extradata)
            faults.append(FaultData(*data))

    if not outfile:
        outfile=os.path.splitext(sourcefile)[0]+'.model'

    with open(outfile,'w') as mf:
        for l in metadata:
            mf.write(l)
        mf.write("fault_num strike_deg dip_deg rake_deg length_km width_km slip_m depth_km east_m north_m\n")
                
        faultno=0
        for f in faults:
            faultno += 1
            depth_km=(f.depth_top_km+f.depth_bottom_km)/2.0
            dip_deg = f.dip_deg
            offset_m=1000.0*depth_km/tan(radians(dip_deg))
            cstr=cos(radians(f.strike_deg))
            sstr=sin(radians(f.strike_deg))
            east_mid_m=round(f.east_surface_m + cstr*offset_m,2)
            north_mid_m=round(f.north_surface_m - sstr*offset_m,2)
            width_km=round((f.depth_bottom_km-f.depth_top_km)/sin(radians(dip_deg)),5)
            mf.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n'.format(
                faultno,
                f.strike_deg,
                dip_deg,
                f.rake_deg,
                f.length_km,
                width_km,
                # Note - Ian Hamlings files seem to use opposite convention for
                # slip vector..
                f.slip_m,
                depth_km,
                east_mid_m,
                north_mid_m
            ))
    return outfile

if __name__ == "__main__":
    import sys
    import glob
    if len(sys.argv) < 2:
        print syntax
        sys.exit()

    metafile = sys.argv[1]

    try:
        print "Extracting model from",metafile
        outfile=extract_model( metafile, sys.argv[2] if len(sys.argv) > 2 else '')
        print "Created",outfile
    except:
        raise
        print sys.exc_info()[1]




