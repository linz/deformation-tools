SOURCE=gridtool.cpp bltmatrx.cpp grid.cpp gridutil.cpp progress.cpp smoothgrid.cpp get_image_path.cpp
OBJ=$(SOURCE:.cpp=.o)

gridtool: $(OBJ)
	g++ $(OBJ) -lm -o gridtool
	
clean: 
	rm $(OBJ)
