#include "func.h"




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
  char *int_file="Gause_data";
  Write_File(int_file, N ,int_kart,int_polar);
  char *time_file="Gause_time";
  Write_File(time_file, N ,time_kart/1000000,time_Polar/1000000);

  return 0;
  }
  // end of main function
