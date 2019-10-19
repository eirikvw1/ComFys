#include <omp.h>
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
#define CHUNKSIZE 1000000

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

      double *mc_int = new double [n];
     double var = 0;
     double MCint, MCintsqr2, Variance;
     MCint=0;
     MCintsqr2=0;
     double funk= 0;
     double a,b;
     a= -2.5;b=2.5;
     MCint = MCintsqr2=0.;
     double invers_period = 1./RAND_MAX; // initialise the random number generator
     srand(time(NULL));  // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method

    int teller =0;
    int chunk = CHUNKSIZE;
    #pragma omp parallel shared(chunk)
    {
    #pragma omp for schedule(dynamic,chunk)
     for ( int i = 0;  i <= n; i++){
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
           funk=int_function(x1,y1,z1,x2,y2,z2);
           MCint +=funk;
           mc_int[i]=funk;
           //cout<<funk<<endl;
           MCintsqr2+=funk*funk;
           teller +=1;

     }
     cout<<"teller "<<teller<<" "<<4000000<<endl;
   }

   
     double e_value = MCint/((double)n);
     for(int i = 0;i<n;i++){
       var+=pow(mc_int[i]-e_value,2);

     }
     var=(var*pow(b-a,6))/n;

     //MCint = sum*(b-a)*(b-a)*(b-a)*(b-a)*(b-a)*(b-a);
     MCint=(MCint*pow((b-a),6))/((double)n); //integral summen multiplisert med jakobi verdien
     MCintsqr2 = (MCintsqr2*pow((b-a),12))/((double) n );

     //Variance=MCintsqr2-MCint*MCint;
     double ana=(5*(M_PI*M_PI))/(16*16);
//   final output
     cout << " variance= " << var << " Integral = " << MCint << " Exact= " << ana << endl;
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
