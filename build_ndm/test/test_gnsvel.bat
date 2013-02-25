@echo off
cd data
echo Y > test_gnsvel.in
echo  -30.243 221.785 0.09735 >> test_gnsvel.in
echo Failed > velocity.out
..\..\gns_velocity\gns_velocity < test_gnsvel.in > ..\out\test_gnsvel.log
move /Y velocity.out ..\out\test_gnsvel.out
del test_gnsvel.in
