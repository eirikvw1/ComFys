# Comment lines
# General makefile for c - choose PROG =   name of given program
# Here we define compiler option, libraries and the  target
CC= g++ -std=c++11 -Wall -Wextra -Wpedantic -larmadillo -lblas -llapack -fopenmp
CFLAGS=-I.





all:pro3_ce par_test
pro3_ce: pro3_ce.cpp
	${CC} -o pro3_ce.exe pro3_ce.cpp
par_test: par_test.cpp
	${CC} -o par_test.exe par_test.cpp
