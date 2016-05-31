=========================
build_grids documentation
=========================

The build_grids script develops a gridded deformation model from a geophysical
fault dislocation model.  The fault model defines a set of rectangular fault planes on which
a slip vector is defined.  These are used to determine the resulting surface 
displacement.

The build_grids script determines a set of gridded displacement models that 
match the displacement field.  The grid spacing is calculated to meet accuracy
requirements.  Accuracy here refers to the difference between the displacements
calculated from the fault model and that interpolated from the grid.

The script attempts to build an efficient set of nested grids that meet the 
accuracy specification.  Note that the accuracy specification take account of whether 
the grid is over land or sea - the specification is less for offshore grid cells.  
The model is effectively truncated at the coast (after buffering). The land area
is defined by a WKT (well known text) multipolygon specification as the single item
in a text file.

The software attempts to build objects with 


Syntax:

        usage: build_grids.py [-h] [--shift-model-path SHIFT_MODEL_PATH]
                              [--submodel-path SUBMODEL_PATH]
                              [--submodel-version SUBMODEL_VERSION] [--subgrids-nest]
                              [--parcel-shift] [--apply-ramp-scale] [--apply-sf-conv]
                              [--test-settings] [--max-level MAX_LEVEL] [--reverse]
                              [--base-tolerance BASE_TOLERANCE] [--split-base]
                              [--no-trim-subgrids] [--clean-dir]
                              [--land-area LAND_AREA] [--precision PRECISION]
                              patch_file model_file [model_file ...]

        Build set of grids for deformation patch

        positional arguments:
          patch_file            Base name used for output files
          model_file            Model file(s) used to calculate deformation, passed to
                                calc_okada

        optional arguments:
          -h, --help            show this help message and exit
          --shift-model-path SHIFT_MODEL_PATH
                                Create a linzshiftmodel in the specified directory
          --submodel-path SUBMODEL_PATH
                                Create publishable component in the specified
                                directory
          --submodel-version SUBMODEL_VERSION
                                Deformation model version for which submodel first
                                applies
          --subgrids-nest       Grid CSV files calculated to replace each other rather
                                than total to deformation
          --parcel-shift        Configure for calculating parcel_shift rather than
                                rigorous deformation patch
          --apply-ramp-scale    Scale the grid by the ramp final value
          --apply-sf-conv       Apply the projection scale factor and convergence to
                                the model calculated displacements
          --test-settings       Configure for testing - generate lower accuracy grids
          --max-level MAX_LEVEL
                                Maximum number of split levels to generate (each level
                                increases resolution by 4)
          --reverse             Published model will be a reverse patch
          --base-tolerance BASE_TOLERANCE
                                Base level tolerance - depends on base column
          --split-base          Base deformation will be split into separate patches
                                if possible
          --no-trim-subgrids    Subgrids will not have redundant rows/columns trimmed
          --clean-dir           Clean publishable component subdirectory
          --land-area LAND_AREA
                                WKT file containing area land area over which model
                                must be defined
          --precision PRECISION

Fault model file format
=======================

The fault model(s) are input files to the [calc_okada](https://github.com/linz/deformation-tools/blob/master/okada/calc_okada.help) program.
The script expects that the first three records in the file are formatted as 

```
   Event: .... d mmm yyyy
   Source model: ....
   Version: ...
```

These metadata will be associated with the published model, and in particular the date is used to build the time model, which
by default is a step from 0 to 1 at the specified time.  Note that the version is the version of the fault
model (generally based on the calculation date).  This is not the same as the version of the
NZGD2000 deformation model in which it is published.

The fault model is expected to have extension .model.  It may be associated with a file .ramp which
defines one or more ramp time models.  Each line in the ramp file is formatted as

```
   yyyy-mm-dd #.###
```

which specifies the factor by which the fault model is multiplied at a specific
date.  The scale factor is assumed to be 0 before the first date.  (If a reverse
patch is being generated this will be automatically shifted to set the factor to
zero after the last date.


Operation of the build_grids program
====================================

The script operates at a high level by first calculating the fault model over a test area
at the coarsest resolution.  It then attempts to work out where that might not meet the
accuracy requirements, and based on that create one or more subgrids to cover those areas.
It repeats this recursively for the subgrids until either the accuracy requirements are
met or the limiting accuracy has been been achieved.

Note here that accuracy requirements are in terms of the accuracy with which the grid
represents the source dislocation model.  This does not relate in any way to the 
accuracy of the source dislocation model itself.  The purpose of this script is 
to calculate grids which represent the source model without introducing additional
significant error by converting the model to a set of grids.  In doing so it 
attempts to keep the total size of the grids small, so that storage, transmission,
and evaluation is efficient.

The mathematical approach used to determine the extents and cell size of the grids
is described in the paper "The Application of a Localised Deformation Model after an Earthquake", 
[Winefield, Crook, and Beavan, 2010](http://www.linz.govt.nz/system/files_force/media/file-attachments/winefield-crook-beavan-application-localised-deformation-model-after-earthquake.pdf)).

The initial origin, size, and spacing are defined in the script.  These are fixed 
so that grids calculated for different models share a common set of points.  Similarly
the subdivision of the grids by factor of 4 is fixed.

Each grid or subgrid is evaluated by the `calc_okada` program using the grid generation
options.  The resulting gridfile is loaded into python DeformationGrid class, defined in
the `degrid.py` script.  This is used to extract information from the grid such 
testing the extents of the grid within which grid parameters exceed a tolerance
(regionsExceedingTolerance), and estimating the grid resolution required to 
meet tolerance requirements (calcResolution).  The calculations also make 
extensive use of the `gridtool` program.

The DeformationGrid.regionsExceedingTolerance function uses the matplotlib contourf 
function to define the region, and then the shapely geometry library to convert this to a 
Polygon that supports spatial operations.

The DeformationGrid.calcResolution function calculates a new grid based on the source
grid which determines the offset at a grid node from the values calculated at 
adjacent nodes, and the grid spacing required to bring this below a tolerance. The
grid spacing is calculated based on a 2nd differential estimate of the error (ie 
on the curvature of the dislocation component fields).  This
function accounts for the precision of the grid data to avoid rounding errors being
misinterpreted as significant grid differences (it is susceptible to rounding
errors as it doing something like numeric differentiation.)

The accuracy of grids is tested against the LINZ network accuracy standards, which
are encoded into the `accuracy_standards.py` script.  

The calculation of patches involves building a number of trial grids of various
resolutions and extents. These are all saved along with the metadata defining 
their source parameters. If the script is rerun then it will check these files
to see if it can be reuse them.  For complex fault models this can save 
a significant amount of time in calculating the gridded surface dislocation values.
The script also builds a file of WKT (well known text) definitions of intermediate
products which can be viewed in a spatial analysis tool such as QGIS
to review the calculation of the deformation.

The grids are built by the build_deformation_grids function, which calculates the
initial grid extents required to represent the deformation, then calls the 
create_grids function recursively to build subgrids within them to 
achieve the required accuracy.  The initial building of the grids is controlled 
by the script configuration parameters listed below.  

Once these grids have been built they are modified by the create_patch_csv script to:

* trim the grids to the extents over which the are required

* add a buffer across which the difference between the subgrid and its parent grid
  is smoothed to zero, to ensure the model defines a continous dislocation function

* if the grids being constructed are to be separated deformation components which
  are added, rather than nested subgrid definitions, then subtract the parent grids
  deformation from the subgrid

* write the final .CSV definition of each grid


Grid calculation configuration
==============================

The main script configuration is defined by settings at the top of the script file.  These
are largely intended to remain unchanged in order to generate consistent, compatible grids.

These parameters include:

* origin, base_size

  define the coarsest grid used.  Keeping these fixed means that grid nodes for different
  patches will match (over the common area of the grids), which makes it easier to combine
  them into products such as reverse shift patches, or epoch dislocation grids

* base_extents

  The maximum extent over which patches will be compiled.  These are curently based
  on the NZ Exclusive Economic Zone

* scale_factor

  Approximate conversion of longitude, latitude to metres, used for accuracy 
  and resolution and buffering calculations..  

* tolerance_factor

  The accuracy requirements of the grid are based on the LINZ network accuracy
  standards.  These are multiplied by a tolerance factor (currently 0.4) to 
  ensure that this does use the full error budget for coordinates of the 
  (see [Winefield et al., 2010](http://www.linz.govt.nz/system/files_force/media/file-attachments/winefield-crook-beavan-application-localised-deformation-model-after-earthquake.pdf)).

The accuracy is assessed by two sets of parameters.  The _base_ parameters apply to the 
initial (coarsest) grid, and the _subcell_ parameters apply when grids are subdivided 
to assess where subdivision is required.  

Two parameters calculated by the calc_okada 
program are used to define the extents over which a dislocation model is required. 
These are the size of the dislocation vector, `ds`, and the part per million change
in an observed vector `err`.

The _base_ parameters include:

* base_limit_test_column
* base_limit_tolerance

  Two parameters define parameters used to define the extents of the grid. The
  test column is either `ds` or `err`.  Whichever is selected is compared with the
  tolerance to define the extents of the patch.  The default settings use the 
  `err` column.  However for calculating shift grids for applying the reverse
  patch in Landonline this was replaced with the `ds` column (ie the size of 
  the dislocation vector, bounded at 0.05metre).

* base_limit_bounds_column
* base_limit_bounds_tolerance
* base_limit_sea_bounds_tolerance

  The bounds column and tolerance are used in a similar way to the test column.
  The intent is that the base_limit_test_column test column is based on
  the part per million component of the network accuracy, whereas the 
  base_limit_bounds_column is based on the constant component of the
  accuracy.  The area of interest is that in which both apply.

* base_ramp_tolerance
* base_ramp_sea_tolerance

  Defines the extents of the ramp applied at the edge of the patch to ensure 
  that there is not a significant discontinuity at the edge of gridded dislocation 
  field.

Applying the base grid parameters will ususally find one or more areas within
which the dislocation model is required to meet the accuracy specificiations.
Each of these is then subdivided as necessary by recursively calling the 
_create_grids_ function.  

The subcell definition in _create_grids_ is controlled by parameters:

* cell_split_factor

  Defines how subgrids are defined in terms of the parent grid.
  Each cell is divided into this number of rows and columns.

* subcell_resolution_tolerance
* subcell_resolution_sea_tolerance

  The subcell resolution tolerances are used in the DeformationGrid.calcResolution
  function.  This determines the size of grid required to model the grid with the
  subcell_resolution_tolerance (the maximum difference between the grid and the
  model dislocations). Where this size is less than the size of the grid being 
  tested the grid is subdivided.

Two factors that control splitting of cells into subcells are:

* max_subcell_ratio

  When a grid is being split into finer subgrids the script has the options
  of replacing the entire grid with a finer grid, rather than replacing some
  parts of it with subgrids.  The subgrids are defined based the variation 
  (curvature) of the dislocation surface.  If the total area of the selected 
  subgrids exceeds the max_subcell_ratio of the parent grid area then the
  parent grid is replaced.

* max_split_level

  This defines the maximum number times the grid will be divided (ie it
  defines the finest grid that will be calculated.

Outputs
=======

In addition to building the set of grids that represent the dislocation the script
can output two other outputs. 

Firstly it can build a patch directory in a [LINZ dislocation model](https://github.com/linz/nzgd2000-deformation-model). 
This copies the final grid .csv files to the patch directory, and creates a component.csv 
file as required by the deformation model format.

Secondly it can create a _LINZ shiftgrid_ file.  This is a file which along with
the grid definitions can be compiled to a binary file representing the reverse patch.  This 
format is used by routines in Landonline to apply the reverse patch to coordinates and
geometry objects.

