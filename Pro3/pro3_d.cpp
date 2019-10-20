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

double int_GauseLA_function(double r1, double r2, double theta1,double theta2,double phi1,double phi2);

//     Main function begins here
int main(int argc, char** argv){
  if((argc <= 1)||(atoi(argv[1])<=0)){
    cout << "Need number argument larger than 0"<< endl;
    return 1;
  }

     int n=atoi(argv[1]);

     double *mc_int = new double [n];
     double MCint, MCintsqr2, Variance,funk;
     double sum= 0;
     double a,b;
     a= -2.5;b=2.5;
     MCint = MCintsqr2=0.;
     double invers_period = 1./RAND_MAX; // initialise the random number generator
     srand(time(NULL));  // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method
     for ( int i = 1;  i <= n; i++){
  // obtain a floating number x in [0,1]
           double r1 = double(rand())*invers_period;
           double r2 = double(rand())*invers_period;
           double theta1 = double(rand())*invers_period;
           double theta2 = double(rand())*invers_period;
           double phi1 = double(rand())*invers_period;
           double phi2 = double(rand())*invers_period;
           r1=-log(1-r1)/4;
           r2=-log(1-r2)/4;
           theta1=0 +(M_PI-0)*theta1;
           theta2=0 +(M_PI-0)*theta2;
           phi1=0 +(2*M_PI-0)*phi1;
           phi2=0 +(2*M_PI-0)*phi2;
           funk=int_GauseLA_function(r1,r2,theta1, theta2, phi1, phi2);
           sum+=funk;
           mc_int[i]=funk;
           //cout<<funk<<endl;
           MCintsqr2+=funk*funk;


     }
     double jac=(4*pow(M_PI,4))/16;
     double e_value = sum/((double)n);
     //MCint = sum*(b-a)*(b-a)*(b-a)*(b-a)*(b-a)*(b-a);
     for(int i = 0;i<n;i++){
      Variance+=pow(mc_int[i]-e_value,2);

     }
     Variance=(Variance*jac)/((double)n);


     MCint=(sum*jac)/((double)n); //integral summen multiplisert med jakobi verdien




     double ana=(5*(M_PI*M_PI))/(16*16);
//   final output
     cout << " variance= " << Variance<< " Integral = " << MCint << " Exact= " << ana << endl;
}  // end of main program
// this function defines the function to integrate

double int_GauseLA_function(double r1, double r2, double theta1,double theta2,double phi1,double phi2){

  double cosbeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
  double r12=r1*r1+r2*r2-2*r1*r2*cosbeta;
        double f = pow(r1,2)*pow(r2,2)*sin(theta1)*sin(theta2)/sqrt(r12);
        if(r12 > ZERO)
                return f;
        else
                return 0;


}// end of function to evaluate
