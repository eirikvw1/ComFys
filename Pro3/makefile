# Comment lines
# General makefile for c - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CC= g++ -std=c++11 -Wall -Wextra -Wpedantic -larmadillo -lblas -llapack -fopenmp
CFLAGS=-I.
DEPS = func.h




all:pro3_ab pro3_cd

pro3_ab: pro3_ab.cpp
	${CC} -o pro3_ab.exe pro3_ab.cpp

pro3_cd: pro3_cd.cpp
	${CC} -o pro3_cd.exe pro3_cd.cpp
