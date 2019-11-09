# Comment lines
# General makefile for c - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CC= g++ -std=c++11 -Wall -Wextra -Wpedantic -larmadillo -lblas -llapack
CFLAGS=-I.
DEPS = metroMain.h metroAlg.h





metroMain: metroMain.cpp metroAlg.cpp
	${CC} -o metroMain.exe metroMain.cpp metroAlg.cpp
