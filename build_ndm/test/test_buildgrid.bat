@echo off
setlocal
cd data
perl ..\..\build_velocity_grid.pl test_buildgrid.gdf ..\out\test_buildgrid.grd > ..\out\test_buildgrid.log

