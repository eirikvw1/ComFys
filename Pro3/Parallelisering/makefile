# Comment lines
# General makefile for c - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CC= mpiCC -std=c++11 -Wall -Wextra -Wpedantic -larmadillo -lblas -llapack
CFLAGS=-I.
DEPS = func.h







all:pro3_ce pro3_de
 pro3_ce: pro3_ce.cpp
	${CC} -o  pro3_ce.exe  pro3_ce.cpp

pro3_de: pro3_de.cpp
	${CC} -o  pro3_de.exe  pro3_de.cpp
