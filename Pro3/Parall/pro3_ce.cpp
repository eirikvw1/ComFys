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
double int_function(double x1, double y1, double z1, double x2, double y2, double z2);

//     Main function begins here
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
  // make sure n > numprocs, else integer division gives zero
  // Length of each process' interval of
  // integration = local_n*h.

  a= -2.5;b=2.5;
  //h = (b-a)/n;
  //local_a = a + my_rank*local_n*h;
  //local_b = local_a + local_n*h;


   // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method
  total_sum = 0.0;


  tie(local_sum, local_vek) = MCintegrering(a, b,n, my_rank,&int_function);



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
    double ana=(5*(M_PI*M_PI))/(16*16);
    variance=variance*pow((b-a),6)/((double)n*4);
    integral = total_sum*pow((b-a),6)/((double)n*4);
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
        funk=(*func)(x1,y1,z1,x2,y2,z2);
        MCint+= funk;
        MCint_vek[i]= funk;

        //printf("SUM= %d \n", MCint);
  }


cout<<"loksl n ="<<n<<endl;
  return make_tuple(MCint, MCint_vek);
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
