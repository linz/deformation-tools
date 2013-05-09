import os
import os.path
import re
import sys
import csv

print '''
build_deformation_model_csv

Interrogates the component.csv files in a deformation format model directory
to build model.csv file.

Takes one parameter which is the name of the model root directory. 

'''

if len(sys.argv) < 2:
    print "Need model root directory as parameter"
    sys.exit()

rootdir = sys.argv[1]
if not os.path.isdir(rootdir):
    print rootdir,"is not a directory"
    sys.exit()

modeldir = os.path.join(rootdir,'model')
if not os.path.isdir(modeldir):
    print rootdir,'does not have a "model" subdirectory'
    sys.exit()

components=sorted([name for name in os.listdir(modeldir) 
            if re.match(r'(ndm|patch_)',name) and
               os.path.isdir(os.path.join(modeldir,name))])

complist=[['component','version_added','version_revoked','reverse_patch','description']]

for c in components:
    print "Processing component:",c
    compcsv = os.path.join(modeldir,c,'component.csv')
    if not os.path.exists(compcsv):
        print "Component csv file",compcsv," is missing"
    with file(compcsv) as ccsvf:
        version_added = '0'
        version_revoked = '0'
        reverse_patch="N"
        description=""
        ccsv = csv.DictReader(ccsvf)
        for item in ccsv:
            if version_added == '0' or item['version_added'] < version_added:
                version_added = item['version_added']
            if item['version_revoked'] > version_revoked:
                version_revoked = item['version_revoked']
            if item['reverse_patch'] == 'Y':
                reverse_patch='Y'
            if description == '':
                description = item['description']
        complist.append([c,version_added,version_revoked,reverse_patch,description])

modelcsv=os.path.join(modeldir,'model.csv')
with file(modelcsv,'w') as mcsvf:
    mcsv = csv.writer(mcsvf)
    for row in complist:
        mcsv.writerow(row)
            

