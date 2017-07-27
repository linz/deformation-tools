#!/usr/bin/python
# Script to build a patch definition file from an NZGD2000 deformation model.


import sys
import argparse
import os
import os.path
import re
import subprocess
from datetime import datetime

sys.path.insert(0,'/home/ccrook/projects_git/python-linz-deformationmodel')
from LINZ.DeformationModel.Model import Model

gridtool=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'gridtool','gridtool')
if not os.path.exists(gridtool):
    raise RuntimeError('Cannot find gridtool program at '+gridtool)

makeshiftpl=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),'perl','makelinzshiftmodel.pl')
if not os.path.exists(makeshiftpl):
    raise RuntimeError('Cannot find makelinzshiftmodel.pl program at '+makeshiftpl)

class _template( object ):

    def __init__( self, definition ):
        definition=re.sub(r'^\s+','',definition,flags=re.M)
        self._definition=definition


    def __call__( self, **args ):
        v=self._definition
        for k in args:
            r=args[k]
            v=v.replace('{'+k+'}',unicode(r))
        return v

patch_header=_template('''
    SHIFT_MODEL_MODEL {patch_name}
    FORMAT LINZSHIFT2B
    VERSION_NUMBER {patch_version}
    COORDSYS NZGD2000
    DESCRIPTION
    {patch_description}
    END_DESCRIPTION
    ''')

patch_component=_template('''

    SHIFT_MODEL_COMPONENT {grid_file}
    MODEL_TYPE grid
    GROUP_ID {group_id}
    FACTOR {factor}
    DESCRIPTION
    {comp_description}
    END_DESCRIPTION
    ''')


output_ordinates={
    ('3d','3d'): ('de','dn','du'),
    ('3d','horizontal'): ('de','dn'),
    ('3d','vertical'): ('du',),
    ('horizontal','3d'): ('de','dn'),
    ('horizontal','horizontal'): ('de','dn'),
    ('vertical','3d'): ('du',),
    ('vertical','vertical'): ('du',),
    }


def main():
    import argparse

    parser=argparse.ArgumentParser("Build reverse patch definition file")
    parser.add_argument('model_dir',help="Model directory")
    parser.add_argument('build_dir',help="Patch build directory")
    parser.add_argument('patch_name',nargs="?",help="Patch file name")
    parser.add_argument('--format',choices=('linzdef','ntv2'),default='linzdef',help="Patch model format")
    parser.add_argument('--version',help="Deformation model for which patch applies, default is current version")
    parser.add_argument('--ordinates',choices=('3d','horizontal','vertical'),default='3d',help="Ordinates required")

    args=parser.parse_args()
    modeldir=args.model_dir
    builddir=args.build_dir
    format=args.format
    if format != 'linzdef':
        raise RuntimeError('Currently only linzdef format is supported')

    if not os.path.isdir(builddir):
        raise RuntimeError('Build directory {0} does not exist or is not a directory'
                           .format(builddir))

    model=Model(modeldir)
    version=args.version
    if version is None:
        version=model.version()
    datumcode=model.metadata('datum_code')
    modelname=model.metadata('model_name')

    patchname=args.patch_name
    if patchname is None:
        patchname='{0}_patch{2}_{1}'.format(datumcode,version,
                        '' if args.ordinates=='3d' else '_'+args.ordinates[:1].upper())

    deffile=os.path.join(builddir,patchname+'.def')
    binfile=os.path.join(builddir,patchname+'.bin')
    gdfname=os.path.join(builddir,patchname+'_g{0:03d}.gdf')

    ngrid=0
    ncomp=0

    with open(deffile,'w') as deff:
        deff.write(patch_header(
            patch_name=modelname,
            patch_version='1.0',
            patch_description=(
                'Model version: '+version + '\n' +
                '\nReverse patch calculated '+
                datetime.today().strftime("%Y-%m-%d"))
            ))

        revcomps=model.reversePatchComponents(args.version)
        for factor,c in revcomps:
            ncomp += 1
            spatial=c.spatialModel
            ordinates=output_ordinates.get((args.ordinates,spatial.displacement_type))
            if ordinates is None:
                continue
            gridfiles=[]
            for gm in spatial.models():
                grid=gm.model()
                if type(grid).__name__ != 'Grid':
                    raise RuntimeError('Only grid models handled by build_patch.py')
                gridfiles.append(model.getFileName(grid.gridFile()))

            # Reverse to ensure highest priority is listed first
            for g in reversed(gridfiles):
                ngrid += 1
                gdf=gdfname.format(ngrid)
                commands=[
                    gridtool,
                    'read','csv',g,
                    'write_linzgrid',datumcode,
                    modelname,
                    'Reverse patch for '+c.name,
                    'Grid file '+os.path.basename(g),
                    'resolution','0.0001',
                    'columns','+'.join(ordinates),
                    gdf]
                subprocess.call(commands)
                deff.write(patch_component(
                    grid_file=os.path.basename(gdf),
                    group_id=ncomp,
                    factor=factor,
                    comp_description="{0}: {1}".format(
                        c.name,os.path.basename(g))
                    ))

    subprocess.call([makeshiftpl,'-f','LINZSHIFT2B',deffile,binfile])


if __name__ == "__main__":
    main()
