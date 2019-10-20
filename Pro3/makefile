# Comment lines
# General makefile for c - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CC= g++ -std=c++11 -Wall -Wextra -Wpedantic -larmadillo -lblas -llapack -fopenmp
CFLAGS=-I.





all:pro3_c2 pro3_d

pro3_c2: pro3_c2.cpp
	${CC} -o pro3_c2.exe pro3_c2.cpp
pro3_d: pro3_d.cpp
	${CC} -o pro3_d.exe pro3_d.cpp
