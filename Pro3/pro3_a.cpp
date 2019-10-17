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

using namespace arma;
using namespace std;

void gauleg(double, double, double *, double *, int);
void WriteVec(double*,double*, int);
void Write_To_File(char *, double*, double*, int);
double Kart_integral(double *W_le,double *x_le,int N);
double int_function(double, double, double, double, double, double);
void gauss_laguerre(double *, double *, int, double);
double Polar_integral(double *r_la,double *Wr_la,double *W_t,double *theta_la,double *W_p,double *phi_la,int N);
double int_GauseLA_function(double r1, double r2, double theta1,double theta2,double phi1,double phi2);
double gammln(double);
void Write_File(char *fil_navn,int,double, double);

int main(int argc, char** argv){
  if((argc <= 1)||(atoi(argv[1])<=0)){
    cout << "Need number argument larger than 0"<< endl;
    return 1;
  }



  int N=atoi(argv[1]);


  //   reserve space in memory for vectors containing the mesh points
//   weights and function values for the use of the gauss-legendre
//   method

  double *x_le = new double [N];
  double *w_le = new double [N];

  double *Wr_la = new double [N+1];
  double *r_la = new double [N+1];
  double *W_t_la = new double [N];
  double *W_p_la = new double [N];
  double *t_la = new double [N];
  double *p_la = new double [N];


  clock_t start, finnish;
  double time_Polar = 0;
  double time_kart = 0;

  start = clock();
  gauleg(-2,2,x_le,w_le,N);
  double int_kart=Kart_integral(w_le,x_le, N);
  finnish = clock();
  time_kart += (finnish-start)/(CLOCKS_PER_SEC/1000000);

  //   set up the mesh points and weights and the power of x^alf
  double alf = 2.0;




  start = clock();
  gauleg(0.0,M_PI,t_la,W_t_la,N);
  gauleg(0.0,2*M_PI,p_la,W_p_la,N);
  gauss_laguerre(r_la,Wr_la, N, alf);
  double int_polar=Polar_integral(r_la,Wr_la,W_t_la,t_la,W_p_la,p_la,N);
  finnish = clock();
  time_Polar += (finnish-start)/(CLOCKS_PER_SEC/1000000);






  double ana=(5*(M_PI*M_PI))/(16*16);
  cout<<"Analytisk: "<<ana<<endl;
  char *int_file="int_data_2";
  Write_File(int_file, N ,int_kart,int_polar);
  char *time_file="time_data_2";
  Write_File(time_file, N ,time_kart,time_Polar);

  return 0;
  }
  // end of main function



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
  myfile.precision(14);
  myfile<<N<<" "<< gauss<<" "<<polar<<endl;
  myfile.close();

}
