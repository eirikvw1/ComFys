#ifndef metroAlg
#define metroAlg


#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <chrono>
#include <string>
#include <chrono>
#include <random>

using namespace  std;
using namespace arma;
// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add){
  return (i+limit+add) % (limit);
}


void initialize(int, double, int **spin_matrix , double&, double&);
// The metropolis algorithm including the loop over Monte Carlo cycles
void  Metropolis(int, long&, int **spin_matrix , double&, double&, double *, int&);
// prints to file the results of the calculations
void initialize(int n_spins ,double temp, int **spin_matrix, double& E, double& M, long& idum);

void Write_File(char *,double,double);

double ran1(long *);

#endif
