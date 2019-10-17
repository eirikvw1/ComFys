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
using namespace std;

//     Here we define various functions called by the main program
//     this function defines the function to integrate

double int_function(double x1, double y1, double z1, double x2, double y2, double z2);

//     Main function begins here
int main(int argc, char** argv){
  if((argc <= 1)||(atoi(argv[1])<=0)){
    cout << "Need number argument larger than 0"<< endl;
    return 1;
  }

     int n=atoi(argv[1]);


     double MCint, MCintsqr2, Variance;
     double sum= 0;
     double a,b;
     a= -2.5;b=2.5;
     MCint = MCintsqr2=0.;
     double invers_period = 1./RAND_MAX; // initialise the random number generator
     srand(time(NULL));  // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method
     for ( int i = 1;  i <= n; i++){
  // obtain a floating number x in [0,1]
           double x1 = double(rand())*invers_period;
           double y1 = double(rand())*invers_period;
           double z1 = double(rand())*invers_period;
           double x2 = double(rand())*invers_period;
           double y2 = double(rand())*invers_period;
           double z2 = double(rand())*invers_period;
           x1=a +(b-a)*x1;
           y1=a +(b-a)*y1;
           z1=a +(b-a)*z1;
           x2=a +(b-a)*x2;
           y2=a +(b-a)*y2;
           z2=a +(b-a)*z2;
           sum+=int_function(x1,y1,z1,x2,y2,z2);
           MCintsqr2=sum*sum;

     }
     //MCint = sum*(b-a)*(b-a)*(b-a)*(b-a)*(b-a)*(b-a);
     MCint=sum*pow((b-a),6);


     MCint = MCint/((double) n );
     MCintsqr2 = MCintsqr2/((double) n );
     double variance=(MCintsqr2-sum*sum);
     double ana=(5*(M_PI*M_PI))/(16*16);
//   final output
     cout << " variance= " << variance << " Integral = " << MCint << " Exact= " << ana << endl;
}  // end of main program
// this function defines the function to integrate

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
