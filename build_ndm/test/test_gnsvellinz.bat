rem - Currently building with g77 from mingw
rem
rem Generates warnings about inconsistent sizes of common block /work/

del test_gnsvellinz.out
..\gns_velocity\gns_velocity_linz.exe data\solution.gns data\lat_long_linz.dat out\test_gnsvellinz.out -30.243 221.785 0.09735 > out\test_gnsvellinz.log
