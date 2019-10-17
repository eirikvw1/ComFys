#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <tuple>
#include "time.h"
#define _USE_MATH_DEFINES
#define EPS 3.0e-14
#define MAXIT 30
#define ZERO 1.0E-10
#include <iomanip>
#include <stdio.h>
//#include <lib.h>
using namespace arma;
using namespace std;

double int_function(double, double, double, double, double, double);

int main(int argc, char** argv){


  int n;
   double x[6], y, fx;
   double z[6];
   double int_mc = 0.;  double variance = 0.;
   double sum_sigma= 0. ; long idum=-1 ;
   double length = 5.; // we fix the max size of the box to L=5
   double jacobidet = pow((2*length),6);
   cout << "Read in the number of Monte-Carlo samples" << endl;
   cin >> n;
//   evaluate the integral with importance sampling
   for ( int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
     for (int j = 0; j< 6; j++) {
         x[j]=-length+2*length*ran0(&idum);
         z[j]=-2.5 +(2.5+2.5)*x[j];
     }
     fx=int_function(x[0],x[1],x[2],x[3],x[4],x[5]);
     int_mc += fx;
     sum_sigma += fx*fx;
   }
   int_mc = int_mc/((double) n );
   sum_sigma = sum_sigma/((double) n );
   variance=sum_sigma-int_mc*int_mc;
//   final output
    cout << setiosflags(ios::showpoint | ios::uppercase);
    cout << " Monte carlo result= " << setw(10) << setprecision(8) << jacobidet*int_mc;
    cout << " Sigma= " << setw(10) << setprecision(8) << volume*sqrt(variance/((double) n )) << endl;


  return 0
}



double int_function(double x1, double y1, double z1, double x2, double y2, double z2){
   double alpha = 2.;
// evaluate the different terms of the exponential
   double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
   double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
   double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
   double f=exp(exp1+exp2)/deno;
   if(deno>ZERO){
     return f;
   }
   else{
     return 0;
   }
} // end of function to evaluate
