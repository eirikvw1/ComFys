#include "metroAlg.h"









void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w, int& teller){
    // loop over all spins

    for(int y  =0; y < n_spins; y++) {
      for (int x= 0; x < n_spins; x++){
        // Find random position
        int ix = (int) (ran1(&idum)*(double)n_spins);
        int iy = (int) (ran1(&idum)*(double)n_spins);
        int deltaE =  2*spin_matrix[iy][ix]*
  	     (spin_matrix[iy][periodic(ix,n_spins,-1)]+
  	      spin_matrix[periodic(iy,n_spins,-1)][ix] +
  	       spin_matrix[iy][periodic(ix,n_spins,1)] +
  	        spin_matrix[periodic(iy,n_spins,1)][ix]);


        // Here we perform the Metropolis test

        if( ran1(&idum) <= w[deltaE+8] ){
  	       spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
          // update energy and magnetization
          M += (double) 2*spin_matrix[iy][ix];
          E += (double) deltaE;
          teller += 1;


          }



    }


  }
} // end of Metropolis sampling over spins


// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temp, int **spin_matrix, double& E, double& M){
    // setup spin matrix and intial magnetization

    for(int y =0; y < n_spins; y++) {
      for (int x= 0; x < n_spins; x++){
        if (temp < 1.5) spin_matrix[y][x] = 1; // spin orientation for the ground state
        M +=  (double) spin_matrix[y][x];
      }
    }
    // setup initial energy
    for(int y =0; y < n_spins; y++) {
      for (int x= 0; x < n_spins; x++){
        E -=  (double) spin_matrix[y][x]*
  	(spin_matrix[periodic(y,n_spins,-1)][x] +
  	 spin_matrix[y][periodic(x,n_spins,-1)]);
      }
    }
  }// end function initialise



void initialize(int n_spins ,double temp, int **spin_matrix, double& E, double& M, long& idum){
      // setup spin matrix and intial magnetization
      //unsigned seed = -(std::chrono::system_clock::now().time_since_epoch().count());
      for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){

          if(ran1(&idum)<=0.5){
            spin_matrix[y][x] = 1; // spin orientation for the ground state
             M +=  (double) spin_matrix[y][x];
          }
          else{
            spin_matrix[y][x] = -1; // spin orientation for the ground state
             M +=  (double) spin_matrix[y][x];
          }

        }
      }
      // setup initial energy
      for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
          E -=  (double) spin_matrix[y][x]*
    	(spin_matrix[periodic(y,n_spins,-1)][x] +
    	 spin_matrix[y][periodic(x,n_spins,-1)]);
        }
      }
}




void Write_File(char *fil_navn,double en,double to){
  ofstream myfile(fil_navn, ios_base::app);
  myfile.precision(8);
  myfile<<to<<" "<<en<<endl;
  myfile.close();

}







#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
int             j;
long            k;
static long     iy=0;
static long     iv[NTAB];
double          temp;

if (*idum <= 0 || !iy) {
   if (-(*idum) < 1) *idum=1;
   else *idum = -(*idum);
   for(j = NTAB + 7; j >= 0; j--) {
      k     = (*idum)/IQ;
      *idum = IA*(*idum - k*IQ) - IR*k;
      if(*idum < 0) *idum += IM;
      if(j < NTAB) iv[j] = *idum;
   }
   iy = iv[0];
}
k     = (*idum)/IQ;
*idum = IA*(*idum - k*IQ) - IR*k;
if(*idum < 0) *idum += IM;
j     = iy/NDIV;
iy    = iv[j];
iv[j] = *idum;
if((temp=AM*iy) > RNMX) return RNMX;
else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
