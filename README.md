Tools for building the NZGD2000 deformation model
=================================================

This contains a number of programs and tools used to build the 
[NZGD2000 deformation model] (https://github.com/linz/nzgd2000-deformation-model).  These are used to build the deformation model, see
the [deformation-files](https://github.com/linz/deformation-files) repository
where the data and scripts used for building the model are held.

The main tools are:

* gridtool

  a C++ program for processing grids, including trimming
  smoothing, etc (see [gridtool.help](https://github.com/linz/deformation-tools/blob/master/gridtool/gridtool.help))

* okada

  a program (calc_okada) for calculating surface deformation from a rectangular
  fault model (see [calc_okada.help](https://github.com/linz/deformation-tools/blob/master/okada/calc_okada.help))

* gns_velocity

  a Fortran program based on GNS code for evaluating the GNS binary velocity model format

Also some python scripts for compiling components.  The most significant are:

* python/build_grids.py

  Complicated script to generate a set of nested grids to represent deformation caused by
  a geophysical dislocation model.  Automatically handles nesting based on accuracy tolerances.
  (see [build_grids.md](https://github.com/linz/deformation-tools/blob/master/python/build_grids.md))

* python/build_model_csv.py

  Builds the model.csv file in the deformation model based directory by interrogating the component.csv
  files in the component subdirectories

* python/buffer_land_extents.py

  Build a multipolygon defining the land extents of New Zealand, over which the deformation
  model patches are required to match the geophysical model (used by build_grids.py).  Input is 
  a shapefile of New Zealand island polygons

* python/buffer_wkt.py

  Similar to buffer_land_extents - used to buffer the plate boundary polygons used to construct the
  offshore national velocity model component.

* python/dump_shift.py

  Converts a reverse patch to ntv2 format using make_ntv2.py, and generates a set of test data
  at semi-randomly selected points to go with it.

