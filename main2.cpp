#include "functions.hpp"

#define dim 7



int main(int argc, char *argv[]){

  double enuc = ReadIN("enuc.dat");
  Matrix overlap = Matrix(dim);
  Matrix ke = Matrix(dim);
  Matrix nai = Matrix(dim);
  ReadIN("s.dat",overlap);
  ReadIN("t.dat",ke);
  ReadIN("v.dat",nai);

  Matrix coreH = Matrix(dim);
  for(int i=0; i < dim; i++){
    for(int j=0; j < dim; j++){
      coreH(i,j) = ke(i,j) + nai(i,j);
    }
  }

  double tei[BIGNUM];
  //std::vector<double> tei;
  ReadIN("eri.dat",tei);  //needs to be configured

  double eigenvalues[dim];
  double eigenvectors[dim*dim];
  
  Diagonalize(eigenvalues,eigenvectors,overlap);

  Matrix Sso = BuildSso(eigenvalues,eigenvectors,coreH);

  Matrix Dens = Matrix(dim);

  double Etot, Eprev;

  do{
    Eprev = Etot;
    Matrix Fock = BuildFock(coreH, Dens, tei);
    
    Matrix C0 = BuildC0(eigenvalues,eigenvectors,Sso,Fock);
    
    Dens = BuildDensity(C0);
    
    Eelec = ComputeSCF(Dens,coreH,Fock);
    Etot = Eelec + enuc;
  
  }while(fabs(Etot - Eprev) > 0.000001);
  
  printf("Energy Value: %f \n", Etot);
  
  return 0;
}

         
         



