#include <mpi.h>
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


using namespace std;
using namespace arma;

//     Here we define various functions called by the main program
//     this function defines the function to integrate
tuple<double ,double *>  MCintegrering(double local_a,double local_b, int local_n, int my_rank,double (*)(double, double, double, double, double, double));
double int_GauseLA_function(double r1, double r2, double theta1,double theta2,double phi1,double phi2);

int main(int nargs, char* args[]){

  int local_n, numprocs, my_rank;
  double a, b, h, local_a, local_b, total_sum, local_sum, local_var, variance,rel_error, integral;
  double  time_start, time_end, total_time;

  MPI_Init (&nargs, &args);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  int n=250000;
  double * local_vek= new double[n];
  int teller = 0;

  time_start = MPI_Wtime();
  local_n = n/numprocs;



  a= -2.5;b=2.5;
  //h = (b-a)/n;
  //local_a = a + my_rank*local_n*h;
  //local_b = local_a + local_n*h;


   // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method
  total_sum = 0.0;


  tie(local_sum, local_vek) = MCintegrering(a, b,n, my_rank,&int_GauseLA_function);


  MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE,
         MPI_SUM, 0, MPI_COMM_WORLD);
  time_end = MPI_Wtime();
  total_time = time_end-time_start;

  //double integral = (total_sum*pow(b-a,6))/((double)n);

  MPI_Barrier (MPI_COMM_WORLD);
  //beregner variansen
  double e_value =0;
  e_value =total_sum/((double)n);
  for(int i =0;i<n;i++){
    local_var+=pow(local_vek[i]-e_value,2);
  }
  MPI_Reduce(&local_var, &variance, 1, MPI_DOUBLE,
         MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  cout<< "rank"<<my_rank<<endl;
  if (my_rank==0 ){
    double jac = 4*pow(M_PI,4)/16;
    double ana=(5*(M_PI*M_PI))/(16*16);
    variance=variance*jac/((double)n*4);
    integral = total_sum*jac/((double)n*4);
    rel_error = abs(ana-integral)/ana;

    cout << " variance= " << variance << " Integral = " << integral << " Exact= " << ana <<"relativ error" <<rel_error<<endl;


    cout << "Time = " <<  total_time
        << " on number of processors: "  << numprocs  << endl;
  }
 // End MPI
 MPI_Finalize ();

 return 0;
}  // end of main program
// this function defines the function to integrate











tuple<double, double *> MCintegrering(double a,double b, int n,int my_rank, double (*func)(double, double, double , double, double, double)){


  cout<<"My rank ="<<my_rank<<"L a "<<a<<endl;
  cout<<"My rank ="<<my_rank<<"L b "<<b<<endl;
  double MCint=0.0;double mean =0.0;double var= 0.0;double funk = 0.0;
  double *MCint_vek= new double [n];

  double invers_period = 1./RAND_MAX; // initialise the random number generator
  srand(time(NULL)*((double)my_rank));


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
          MCint+=funk;
          MCint_vek[i]=funk;
          //cout<<funk<<endl;



    }




  return make_tuple(MCint, MCint_vek);
}





double int_GauseLA_function(double r1, double r2, double theta1,double theta2,double phi1,double phi2){

  double cosbeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
  double r12=r1*r1+r2*r2-2*r1*r2*cosbeta;
        double f = pow(r1,2)*pow(r2,2)*sin(theta1)*sin(theta2)/sqrt(r12);
        if(r12 > ZERO)
                return f;
        else
                return 0;


}//
