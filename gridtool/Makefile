SOURCE=gridtool.cpp bltmatrx.cpp grid.cpp gridutil.cpp progress.cpp smoothgrid.cpp get_image_path.cpp
OBJ=$(SOURCE:.cpp=.o)

CXXFLAGS=-g

gridtool: $(OBJ)
	g++ -g $(OBJ) -lm -o gridtool
	
clean: 
	rm $(OBJ)
