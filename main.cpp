#include "functions.hpp"
#include "global.hpp"

#define dim 7  // Particular scenario 7 dimensional 

int ioff[10000];  // Indexing array

int main(int argc, char *argv[]){
 
  // Populate indexing array
  ioff[0] = 0;
  for(int i=1; i < 10000; i++)
    ioff[i] = ioff[i-1] + i;

  // Declare and read in nuclear repulsion energy,
  // AO-basis overlap, and one-electron integrals
  // (kinetic energy / nuclear-attraction).
  double enuc = ReadIN("enuc.dat");

  Matrix overlap = Matrix(dim);
  Matrix ke = Matrix(dim);
  Matrix nai = Matrix(dim);
  ReadIN("s.dat",overlap);
  ReadIN("t.dat",ke);
  ReadIN("v.dat",nai);

  // Declare and read in two-electron repulsion integrals.
  double tei[10000];
  ReadIN("eri.dat",tei);  
  

  // Construct core hamiltonian using kinetic energy
  // and nuclear attraction integrals.
  Matrix coreH = Matrix(dim);
  for(int i=0; i < dim; i++){
    for(int j=0; j < dim; j++){
      coreH(i,j) = ke(i,j) + nai(i,j);
    }
  }

  // Diagonalize AO-basis overlap and construct
  // orthodgoanlization matrix Sso 
  double eigenvalues[dim];
  double eigenvectors[dim*dim];
  
  Diagonalize(eigenvalues,eigenvectors,overlap);

  Matrix Sso = BuildSso(eigenvalues,eigenvectors,coreH);


  // Define and initialize density matrix,
  // SCF electroninc energy, total energy.
  Matrix Dens = Matrix(dim);
  double Eelec, Etot, Eprev;

  ///////////////////////////
  // Iteration procedure: ///
  ///////////////////////////
  // - Construct new fock matrix using core hamiltonian,
  // previous density matrix (first iteration dens = 0),
  // and two-electron integrals.
  // - Diagonalize the new fock matrix and back-transform
  // using orthoganalization matrix Sso.
  // - Construct new density matrix and compute SCF electron
  // energy.  Addition of nuclear repulsion energy gives
  // the total energy.
  // - Perform iteration procedure until current total energy
  // and previous converge.
  //

  do{

    Eprev = Etot;  // Previous iteration total energy

    Matrix Fock = BuildFock(coreH, Dens, tei);  // Construct new fock
   
    Matrix C0 = BuildC0(eigenvalues,eigenvectors,Sso,Fock); // Construct new C0
    
    Dens = BuildDensity(C0);  // Construct new density matrix
    
    Eelec = ComputeSCF(Dens,coreH,Fock);
    Etot = Eelec + enuc;
    printf("Etot: %f \n", Etot);
  
  }while(fabs(Etot - Eprev) > 0.000001);  // Test for convergence
  
  printf("Energy Value: %f \n", Etot);
  
  return 0;
}

         
         



