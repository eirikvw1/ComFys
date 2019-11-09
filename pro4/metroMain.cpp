/*
   Program to solve the two-dimensional Ising model
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
*/


#include "metroMain.h"






//  main program
int main(int argc, char* argv[]){
  char *outfilename;
  long idum;
  int **spin_matrix;
  int n_spins, mcs, teller;
  double w[17],  average[5], initial_temp, final_temp, E, M, temp_step;

  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }

  //    Read in initial values such as size of lattice, temp and1Z∑i=1MEie−βE cycles
  //n_spins, msc, intital_temp, final temp, temp_step;
  n_spins=atoi(argv[2]);mcs=atoi(argv[3]);initial_temp=atoi(argv[4]);
  final_temp=atoi(argv[5]);temp_step=atoi(argv[6]);

  spin_matrix= new int *[n_spins];
  for (int i = 0; i < n_spins; ++i){
    spin_matrix[i] = new int[n_spins];

  }


  idum = -1; // random starting point
  for ( double temp = initial_temp; temp <= final_temp; temp+=temp_step){
    //    initialise energy and magnetization
    E = M = 0.;teller =0;

    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
    // initialise array for expectation values
    for( int i = 0; i < 5; i++) average[i] = 0.;
    initialize(n_spins, temp, spin_matrix, E, M);
    cout<<"intialisert enrgiern"<< E <<endl;
    // start Monte Carlo computation

    for (int cycles = 1; cycles <= mcs; cycles++){

      Metropolis(n_spins, idum,spin_matrix, E, M, w, teller);
      // update expectation values

      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);

      double E_expect= ((double) average[0])/((double) cycles);
      double EE_expect= ((double) average[1])/((double) cycles);
      double Evar= (EE_expect-pow(E_expect,2))/(pow(temp,2));
      double M_expect= ((double) average[2])/((double)cycles);
      double MM_expect= ((double) average[3])/((double) cycles);
      double M_abs_expect= ((double) average[4])/((double) cycles);

      //Write_File(fil_navn,E_expect, M_abs_expect);


    }
    // print results

    double E_expect= ((double) average[0])/((double) mcs);
    double EE_expect= ((double) average[1])/((double) mcs);
    double Evar= (EE_expect-pow(E_expect,2))/(pow(temp,2));
    double M_expect= ((double) average[2])/((double) mcs);
    double MM_expect= ((double) average[3])/((double) mcs);
    double M_abs_expect= ((double) average[4])/((double) mcs);
    cout << " L= " << n_spins << "; mcs = " << mcs<<endl;
    double suscept =(MM_expect -pow(M_abs_expect,2))/(temp);
    cout <<" E= "<<average[0]<<" M= "<<average[2]<<endl;
    //free_matrix(spin_matrix); // free memory
      cout <<" E expect= "<<E_expect<<" e_var= "<<Evar<<endl;
      cout <<" M expect= "<<M_expect<<" abs(M) expexted "<< M_abs_expect<<"susceptebilety = "<< suscept<<endl;
      cout <<" E expect= "<<E_expect<<endl;
      cout << " teller sum " << teller<<endl;
  }
  //Write_File(fil_navn,teller,mcs);
  return 0;
}







/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */
