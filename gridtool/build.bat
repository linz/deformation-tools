setlocal
cd /d %~dp0
del /q *.obj
cl /EHsc -Dstrncasecmp=strnicmp gridtool.cpp grid.cpp gridutil.cpp smoothgrid.cpp bltmatrx.cpp progress.cpp get_image_path.cpp
del /q *.obj
