#include <armadillo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <tuple>
#include "time.h"
#define _USE_MATH_DEFINES
#define EPS 3.0e-14
#define MAXIT 30
#define ZERO 1.0E-10
using namespace arma;
using namespace std;
#ifndef func
#define func

void gauleg(double, double, double *, double *, int);
void WriteVec(double*,double*, int);
void Write_To_File(char *, double*, double*, int);
double Kart_integral(double *W_le,double *x_le,int N);
//double int_function(double, double, double, double, double, double);
void gauss_laguerre(double *, double *, int, double);
double Polar_integral(double *r_la,double *Wr_la,double *W_t,double *theta_la,double *W_p,double *phi_la,int N);
double int_GauseLA_function(double r1, double r2, double theta1,double theta2,double phi1,double phi2);
double gammln(double);
void Write_File(char *fil_navn,int,double, double);
void Write_File(char *fil_navn, int,double, double, double, double);
double polar_MC_function(double r1, double r2, double theta1,double theta2,double phi1,double phi2);
double kart_function(double x1, double y1, double z1, double x2, double y2, double z2);
void MC_kart_integral(double *, double *, double *, double *,int);
void MC_Polar_integral(double *MC_int, double *var, int n);



#endif




  void gauleg(double x1, double x2, double x[], double w[], int n){
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */

	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()


void gauss_laguerre(double *x, double *w, int n, double alf){
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}

double gammln( double xx){
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}


void WriteVec(double* vec1,double* vec2, int n){
  for(int i=0;i < n;i++){
    cout<<vec1[i]<<"  "<<vec2[i]<<endl;
  }

}
void Write_To_File(char *fil_navn, double* x, double* w, int n){
  ofstream myfile;
  myfile.open(fil_navn);
  myfile.precision(14);

  for(int j = 0; j<n; j++){
    myfile<<x[j]<<" "<<w[j]<<endl;
  }
  myfile.close();
  return;
}
//  this function defines the function to integrate
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


double Kart_integral(double *w_le,double *x_le,int N){


  double int_gauss = 0.0;
//   six-double loops
  for (int i=0;i<N;i++){
    for (int j = 0;j<N;j++){
    for (int k = 0;k<N;k++){
    for (int l = 0;l<N;l++){
    for (int m = 0;m<N;m++){
    for (int n = 0;n<N;n++){
     int_gauss+=w_le[i]*w_le[j]*w_le[k]*w_le[l]*w_le[m]*w_le[n]
    *int_function(x_le[i],x_le[j],x_le[k],x_le[l],x_le[m],x_le[n]);
  }}}}}}
  return int_gauss;

}

double int_GauseLA_function(double r1, double r2, double theta1,double theta2,double phi1,double phi2){

  double cosbeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
  double r12=r1*r1+r2*r2-2*r1*r2*cosbeta;
        double f = exp(-3*(r1+r2))*sin(theta1)*sin(theta2)/sqrt(r12);
        if(r12 > ZERO)
                return f;
        else
                return 0;


}
double Polar_integral(double *r_la,double *Wr_la,double *W_t,double *theta_la,double *W_p,double *phi_la,int N){


  double int_gauss_polar =0.0;


  for (int i=1;i<N+1;i++){
    for (int j = 1;j<N+1;j++){
    for (int k = 0;k<N;k++){
    for (int l = 0;l<N;l++){
    for (int m = 0;m<N;m++){
    for (int n = 0;n<N;n++){
     int_gauss_polar+=Wr_la[i]*Wr_la[j]*W_t[k]*W_t[l]*W_p[m]*W_p[n]
    *int_GauseLA_function(r_la[i],r_la[j],theta_la[k],theta_la[l],phi_la[m],phi_la[n]);
  }}}}}}
  return int_gauss_polar;

}

void Write_File(char *fil_navn, int N,double gauss, double polar){
  ofstream myfile(fil_navn, ios_base::app);
  myfile.precision(8);
  myfile<<N<<" "<< gauss<<" "<<polar<<endl;
  myfile.close();

}
void Write_File(char *fil_navn, int N,double kart, double polar, double var_kart, double var_polar){
  ofstream myfile(fil_navn, ios_base::app);
  myfile.precision(8);
  myfile<<N<<" "<< kart<<" "<<polar<< " "<<var_kart <<" "<<var_polar <<endl;
  myfile.close();

}

double kart_function(double x1, double y1, double z1, double x2, double y2, double z2){
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


double polar_MC_function(double r1, double r2, double theta1,double theta2,double phi1,double phi2){

  double cosbeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
  double r12=r1*r1+r2*r2-2*r1*r2*cosbeta;
        double f = pow(r1,2)*pow(r2,2)*sin(theta1)*sin(theta2)/sqrt(r12);
        if(r12 > ZERO)
                return f;
        else
                return 0;


}// end of function to evaluate

void MC_kart_integral(double *MCint, double *var, double *a, double *b, int n){
    double invers_period = 1./RAND_MAX; // initialise the random number generator
    srand(time(NULL));
    // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method
      double e_value =0;double funk= 0;
      double *mc_int = new double [n];
      double sum_MC =0;
      //start for loop
      for ( int i = 0;  i <= n; i++){

    // obtain a floating number x in [0,1]
             double x1 = double(rand())*invers_period;
             double y1 = double(rand())*invers_period;
             double z1 = double(rand())*invers_period;
             double x2 = double(rand())*invers_period;
             double y2 = double(rand())*invers_period;
             double z2 = double(rand())*invers_period;
             x1=*a +(*b-*a)*x1;
             y1=*a +(*b-*a)*y1;
             z1=*a +(*b-*a)*z1;
             x2=*a +(*b-*a)*x2;
             y2=*a +(*b-*a)*y2;
              z2=*a +(*b-*a)*z2;
             funk=kart_function(x1,y1,z1,x2,y2,z2);
             sum_MC +=funk;
             mc_int[i]=funk;
       }// end of For loop

       e_value = (*MCint)/((double)n);
       double sum_var=0;
       //start for loop
       for(int i = 0;i<n;i++){
         sum_var+=pow(mc_int[i]-e_value,2);

       }// end of for loop
       *MCint=sum_MC;
       *var=sum_var;
       delete mc_int;

}//end of fuktion


void MC_Polar_integral(double *MC_int, double *var,int n){
  double invers_period = 1./RAND_MAX; // initialise the random number generator
  srand(time(NULL));
  // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method
    double e_value =0;double funk= 0;
    double *mc_int = new double [n];
    double sum_MC =0;
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
          funk=polar_MC_function(r1,r2,theta1, theta2, phi1, phi2);
          sum_MC+=funk;
          mc_int[i]=funk;
          //cout<<funk<<endl;



    }// end of loop


  e_value = (sum_MC)/((double)n);
  double sum_var=0;
  //start for loop
  for(int i = 0;i<n;i++){
    sum_var+=pow(mc_int[i]-e_value,2);

  }// end of for loop
  *MC_int=sum_MC;
  *var=sum_var;
  delete mc_int;
}// end of function MC_Polar_integral
