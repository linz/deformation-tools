
CPPFLAGS=-std=c++0x -DUSE_BOOST_REGEX -D_stricmp=strcasecmp
LDLIBS=-lboost_regex -lm

all: calc_okada 
#test_okada

calc_okada: calc_okada.o okada.o get_image_path.o

test_okada: test_okada.o okada.o

clean: 
	rm calc_okada.o okada.o test_okada.o get_image_path.o
