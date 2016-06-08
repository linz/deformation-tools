setlocal
cd /d %~dp0
del /q *.obj
cl /EHsc /TP calc_okada.cpp okada.cpp get_image_path.cpp tmproj.c

del /q *.obj
