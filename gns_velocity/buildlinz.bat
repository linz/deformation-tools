rem - Currently building with g77 from mingw
rem
rem Generates warnings about inconsistent sizes of common block /work/

del gns_velocity_linz.exe
g77 gns_velocity_linz.f -o gns_velocity_linz.exe > gns_velocity.err 2>&1
