#include <mpi.h>
#include "func.h"


//     Here we define various functions called by the main program
//     this function defines the function to integrate
tuple<double ,double *>  MC_plaraisering_polar_integrering(double local_a,double local_b, int local_n, int my_rank,double (*)(double, double, double, double, double, double));


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


  local_n = n/numprocs;
  a= -2.5;b=2.5;
  //h = (b-a)/n;
  //local_a = a + my_rank*local_n*h;
  //local_b = local_a + local_n*h;


   // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method
  total_sum = 0.0;

  time_start = MPI_Wtime();
  tie(local_sum, local_vek) = MC_plaraisering_polar_integrering(a, b,n, my_rank,&polar_MC_function);
  time_end = MPI_Wtime();
  total_time = time_end-time_start;

  MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE,
         MPI_SUM, 0, MPI_COMM_WORLD);


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




  if (my_rank==0 ){
    double jac = 4*pow(M_PI,4)/16;
    double ana=(5*(M_PI*M_PI))/(16*16);
    variance=variance*jac/((double)n*4);
    integral = total_sum*jac/((double)n*4);
    rel_error = abs(ana-integral)/ana;

    char *datafile="Paralel_polar_data.text";

    Write_File(datafile,n,integral,variance, rel_error,total_time);

    cout << " variance= " << variance << " Integral = " << integral << " Exact= " << ana <<"relativ error" <<rel_error<<endl;


    cout << "Time = " <<  total_time
        << " on number of processors: "  << numprocs  << endl;
  }
 // End MPI
 MPI_Finalize ();

 return 0;
}  // end of main program
// this function defines the function to integrate
