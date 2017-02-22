
CPPFLAGS=-std=c++0x -DUSE_BOOST_REGEX -DUSE_BOOST_STRINGS -D_stricmp=strcasecmp
LDLIBS=-lboost_regex -lm 
CC=g++

all: calc_okada 
#test_okada

calc_okada: calc_okada.o okada.o get_image_path.o tmproj.o

test_okada: test_okada.o okada.o

clean: 
	rm -f calc_okada.o okada.o test_okada.o get_image_path.o tmproj.o
