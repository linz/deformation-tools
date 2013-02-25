#!/usr/bin/python
import sys
import xlrd
import os
import os.path
import re

help='''

extract_gns_model:  Generates a model file suitable for the calc_okada program
from the spreadsheet format provided by GNS

Syntax:
python extract_gns_model.py gns_spreadsheet [root_filename]

gns_spreadsheet is the .xsl file provided by GNS
root_file_name is the output file base name

Two files will be generated
   root_file_name.model      Containing the model definition
   root_file_name.test.csv   CSV file containing the test data from the 
                             spreadsheet that can be used with calc_okada
                             to verify the data
'''

def extract_model( xls_file, outfile=None ):
    book = xls_file
    if not outfile:
        outfile= os.path.splitext(book)[0]
    
    filename = os.path.basename(book)
    if not os.path.exists(book):
        print "Spreadsheet ",book," does not exist"
        sys.exit()
    wb = None
    ms = None
    ds = None
    try:
        wb = xlrd.open_workbook( book, encoding_override='cp1252' )
        ms = wb.sheet_by_name('source_model')
        ds = wb.sheet_by_name('displacements')
    except:
        pass
    if not wb:
        raise RuntimeError( book+": is not a spreadsheet" )
    if not ms or not ds:
        raise RuntimeError(book+": does not contain source_model and displacements tabs")

    header = []
    model = []
    lat0 = None
    lon0 = None
    testdata = []
    
    in_data = False
    re_data = re.compile(r'strike_deg\s+dip_deg\s+rake_deg\s+length_km\s+width_km\s+')
    
    for i in xrange(ms.nrows):
        r = ms.row(i)
        rdata = [unicode(c.value) for c in r]
        record = ' '.join(rdata).strip()
        if not in_data and re_data.search(record):
            in_data = True
            model.append('\t'.join(rdata))
            continue
        if in_data:
            model.append('\t'.join([re.sub(r'\.0$','',r) for r in rdata]))
            continue
        if re.match(r'\s*$',record):
            continue
        header.append(record)
    
    in_data = False
    cols = []
    
    for i in xrange(ds.nrows):
        r = ds.row(i)
        rdata = [unicode(c.value) for c in r]
        record = '\t'.join(rdata).strip()
        match = re.match(r'\s*(lat0|lon0)\s*\=\s*(-?\d+\.?\d*)\s*$',record)
        if match:
            if match.group(1) == 'lat0':
                lat0 = match.group(2)
            else:
                lon0 = match.group(2)
            continue
        if re.match(r'\s*site_code\s+lat_deg\s+lon_deg',record):
            in_data = True
            cols = [rdata.index(i) for i in 'site_code lon_deg lat_deg east_mm south_mm up_mm'.split(' ')]
            testdata.append('code\tlon\tlat\tux\tuy\tuz')
            continue
        if in_data:
            tdata = [rdata[i] for i in cols]
            tdata[0] = re.sub('\.0$','',tdata[0])
            tdata[3] = str(float(tdata[3])/1000)
            tdata[4] = str(-float(tdata[4])/1000)
            tdata[5] = str(float(tdata[5])/1000)
            testdata.append('\t'.join(tdata))
    
    if len(header) < 4:
        raise RuntimeError(book + ": model header information missing")
    
    try:
        modfile = open(outfile+'.model','w')
        modfile.write("Event: "+header[0]+"\n")
        modfile.write("Model: "+header[1]+"\n")
        modfile.write("Version: "+header[2]+"\n")
        modfile.write("Author: "+header[3]+"\n")
        modfile.write("SourceFile: "+filename+"\n")
        modfile.write("Origin: "+str(lon0)+' '+str(lat0)+"\n")
        for line in model:
            modfile.write(line+"\n")
        modfile.close()
        
        testfile = open(outfile+'.test.csv','w')
        for line in testdata:
            testfile.write(line+"\n")
        testfile.close()
    except:
        if os.path.exists(outfile+'.model'):
            os.unlink(outfile+'.model')
        if os.path.exists(outfile+'.test.csv'):
            os.unlink(outfile+'.test.csv')
        raise
    

if __name__ == "__main__":
    import sys
    import glob
    if len(sys.argv) < 2:
        print help
        sys.exit()

    xls = sys.argv[1]
    if not xls.lower().endswith('.xls'):
        xls += ".xls"
    if len(sys.argv) == 2:
        for g in glob.glob(xls):
            try:
                print "Extracting model from",g
                extract_model(g)
            except RuntimeError as e:
                print e.args[0]

    else:
        try:
            print "Extracting model from",xls
            extract_model( xls, sys.argv[2] )
        except:
            print e.args[0]




