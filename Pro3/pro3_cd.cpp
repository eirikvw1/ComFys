#include "func.h"

//     Here we define various functions called by the main program
//     this function defines the function to integrate


// I dette skriptet oppgave c og d gjort
//     Main function begins here
int main(int argc, char** argv){
  if((argc <= 1)||(atoi(argv[1])<=0)){
    cout << "Need number argument larger than 0"<< endl;
    return 1;
  }

     int n=atoi(argv[1]);
     double MCint_polar, MCint_kart,var_kart,a,b,var_polar, jacobi_kart, jacobi_polar;
     a= -2.5;b=2.5;var_kart=0;var_polar=0;jacobi_kart=0;jacobi_polar=0;
     clock_t start, finnish;
     double time_polar = 0;
     double time_kart = 0;


     start =clock();
     // beregner Montevarlo integralet for kartesiske kordinater
     MC_kart_integral(&MCint_kart, &var_kart,&a,&b,n);
     finnish = clock();
     time_kart = (finnish-start)/(CLOCKS_PER_SEC/1000000);
     // beregner ferdig integral sum og variancen for kartesiske kordinatsystem
     jacobi_kart=pow((b-a),6);
     MCint_kart=(MCint_kart*jacobi_kart)/((double)n); //integral summen
     var_kart=(var_kart*jacobi_kart)/((double)n);// variancen


      start =clock();
     //beregner Montecarlo integralet med polare kordinater
     MC_Polar_integral(&MCint_polar,&var_polar, n);
     finnish = clock();
     time_polar = (finnish-start)/(CLOCKS_PER_SEC/1000000);
     //ana variabel er analytisk løsning
     jacobi_polar=(4*pow(M_PI,4))/16;;
     MCint_polar=(MCint_polar*jacobi_polar)/((double)n); //integral summen
     var_polar=(var_polar*jacobi_polar)/((double)n);// variancen
     //analytisk losning
     double ana=(5*(M_PI*M_PI))/(16*16);
     // writing data to two files, time file and data file
     char *timefile="MC_time.text";
     char *datafile="MC_data.text";
     Write_File(timefile,n, time_kart/1000000,time_polar/1000000);
     Write_File(datafile,n,MCint_kart,MCint_polar, var_kart,var_polar);

     //printer ut verdier
     cout << " variance cart=  " << var_kart << " Integral cart = " << MCint_kart << " Exact= " << ana <<" time sek= "<<time_kart/1000000.0<< endl;
     cout << " variance polar= " << var_polar << " Integral polar = " << MCint_polar << " Exact= " << ana <<" time sek= " <<time_polar/1000000.0 <<endl;
     return 0;
}  // end of main program
