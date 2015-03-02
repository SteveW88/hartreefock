#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "matrices.hpp"



void Diagonalize(double eigenval[], double eigenvec[], double origMat[]){  //Set for 7X7 matrices

  gsl_matrix_view m 
    = gsl_matrix_view_array (origMat, 7, 7);

  gsl_vector_complex *eval = gsl_vector_complex_alloc (7);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (7, 7);
  
  gsl_eigen_nonsymmv_workspace * w = 
    gsl_eigen_nonsymmv_alloc (7);
  
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);
  
  gsl_eigen_nonsymmv_free (w);
  
  gsl_eigen_nonsymmv_sort (eval, evec, 
                           GSL_EIGEN_SORT_ABS_DESC);
  
  {
    int i, j;
    
    for (i = 0; i < 7; i++)
      {
        gsl_complex eval_i 
          = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i 
          = gsl_matrix_complex_column (evec, i);
                     
         printf ("eigenvalue = %g \n",
                 GSL_REAL(eval_i));
         printf ("eigenvector = \n");   
        eigenval[i] = GSL_REAL(eval_i);
        
        for (j = 0; j < 7; ++j)
          {
            gsl_complex z = 
              gsl_vector_complex_get(&evec_i.vector, j);
            printf("%g \n", GSL_REAL(z));
            eigenvec[(i*7)+j] = GSL_REAL(z);
          }
      }
  }
  
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
}


int main(int argc, char *argv[]) {
  
  // Array Indexing
  int BIGNUM = 10000;
  int ioff[BIGNUM];
  ioff[0] = 0;
  for(int i = 1; i < BIGNUM; i++)
    ioff[i] = ioff[i-1] + i;


  int i, j, k, l, ij, kl, ijkl, ik, jl, ikjl;
  int ji;
  // Row/Column to one dimensional array indexing transformation
  // ij = (i>j) ? ioff[i] + j : ioff[j] + i;
  // ijkl = (ij > kl) ? ioff[ij] + kl : ioff[kl] + ij; //TODO: Verify

  //Step #1: Read enuc
  //
  std::ifstream enucFile;
  enucFile.open("enuc.dat",std::ifstream::in);
  double enuc;
  enucFile >> enuc;
  enucFile.close();
 
  
  //Step #2: Read overlap
  //
  FILE *input;
  double val;
  double overlap[49];
  double overlap2[7][7];
  input = fopen("s.dat", "r");
  while(fscanf(input, "%d %d %lf", &i, &j, &val) != EOF){
    i -= 1;
    j -= 1;
    overlap[(i*7)+j] = val;
    if(i != j)
      overlap[(j*7)+i] = val;
    

    overlap2[i][j] = val;
    if (i != j)
      overlap2[j][i] = val;
  }
  fclose(input);

  //Read kinetic energy
  //
  double ke[33];
  double ke2[7][7];
  input = fopen("t.dat", "r");
  while(fscanf(input, "%d %d %lf", &i, &j, &val) != EOF){
    i -= 1;
    j -= 1;
    ij = (i>j) ? ioff[i] + j : ioff[j] + i;
    ke[ij] = val;
    ke2[i][j] = val;
    if (i != j)
      ke2[j][i] = val;
  }
  fclose(input);
  

  //Read nuclear-attraction integrals
  //
  double nai[BIGNUM];
  double nai2[7][7];
  input = fopen("v.dat", "r");
  while(fscanf(input, "%d %d %lf", &i, &j, &val) != EOF){
    i -= 1;
    j -= 1;
    ij = (i>j) ? ioff[i] + j : ioff[j] + i;
    nai[ij] = val;
    nai2[i][j] = val;
    if(i != j)
      nai2[j][i] = val;
  }
  fclose(input);

  //Core Hamiltonian
  //
  double coreH[7][7];
  double coreH2[49];
  for(int i = 0; i < 7; i++){
    for(int j = 0; j < 7; j++){
      coreH[i][j] = ke2[i][j] + nai2[i][j];
      coreH2[(i*7)+j] = ke2[i][j] + nai2[i][j];
    }
  }
  
 
  //Step 3: Read two-electron integrals
  //
  double tei[BIGNUM];
  input = fopen("eri.dat", "r");
  while(fscanf(input, "%d %d %d %d %lf", &i, &j, &k, &l, &val) != EOF){
    i -= 1;
    j -= 1;
    k -= 1;
    l -= 1;
    
    ij = (i>j) ? ioff[i] + j : ioff[j] + i;
    kl = (k>l) ? ioff[k] + l : ioff[k] + l;
    ijkl = (ij > kl) ? ioff[ij] + kl : ioff[kl] + ij;
    tei[ijkl] = val;
  
  }
  fclose(input);
 
  ///////////////////////
  //Step 4: Diagonalize overlap
  ///////////////////////
  
  double eigenvalues[7];
  double eigenvectors[49];
 
  Diagonalize(eigenvalues,eigenvectors,overlap);
  Matrix Hcore = Matrix(coreH2,7);  // Correct
  DiagMat As = DiagMat(eigenvalues,7);
  Matrix Ls = Matrix(eigenvectors,7).T();
  Matrix LsT = Ls.T();
  Matrix Sso = ((Ls*As.ToPower(-0.5))*Ls.T()); // Correct
  

  //Step #5
  //
  Matrix Fockinit = (Sso.T() * Hcore * Sso); // Correct

  double * Fockinitarray = Fockinit.ToArraySTATIC();
  double Farray[49];                        // Correct
  for(int i=0; i < 49; i++){
    Farray[i] = *(Fockinitarray+i);
  }
  Diagonalize(eigenvalues,eigenvectors,Farray);



  Matrix C0prime = Matrix(eigenvectors,7).T();
  DiagMat eps0 = DiagMat(eigenvalues,7);
  Matrix C0 = (Sso * C0prime);  // C0 values match in value, not in sign



  // build density matrix
  Matrix DensInit = Matrix(7);
  for(int k=0; k < 7; k++){
    for(int j=0; j < 7; j++){
      for(int i=0; i <4; i++){ // Set to 4 works (sum from 0 to N/2)
        DensInit(k,j) += C0(k,i) * C0(j,i);   // Correct
      }
    }
  }

  // Step 6: Compute Initial SCF Energy
  //
  double Eelec = 0;
  for(int i=0; i<7; i++){
    for(int j=0; j<7; j++){
      Eelec +=  DensInit(i,j) * (Hcore(i,j) + Fockinit(i,j));
    }
  }
 
  double Etot = Eelec + enuc;  //incorret
  double Etotprev;
  

  // Step 7: Compute New Fock Matrix
  //
  do{
    Etotprev = Etot;
#define INDEX(i,j) (i>j) ? (ioff[i]+j) : (ioff[j]+i);
  double F[7][7];
  for(int i=0; i < 7; i++)
    for(int j=0; j < 7; j++) {
      F[i][j] = coreH[i][j];
      for(int k=0; k < 7; k++)
        for(int l=0; l < 7; l++) {
          ij = INDEX(i,j);
          kl = INDEX(k,l);
          ijkl = INDEX(ij,kl);
          ik = INDEX(i,k);
          jl = INDEX(j,l);
          ikjl = INDEX(ik,jl);
 
          F[i][j] += DensInit(k,l) * (2.0 * tei[ijkl] - tei[ikjl]);
        }
    }
  double F1[49];
  for(int i=0; i < 7; i++){
    for(int j=0; j<7; j++){
      F1[(i*7)+j] = F[i][j];
    }
  }
  Matrix newFock = Matrix(F1,7);
  // Step 8: Build New Density Matrix
  //
  Matrix newFockinit = (Sso.T() * newFock * Sso); 

  double * newFockinitarray = newFockinit.ToArraySTATIC();
  double newFarray[49];                       
  for(int i=0; i < 49; i++){
    newFarray[i] = *(newFockinitarray+i);
  }
  Diagonalize(eigenvalues,eigenvectors,newFarray);
     
  Matrix newC0prime = Matrix(eigenvectors,7).T();
  DiagMat neweps0 = DiagMat(eigenvalues,7);
  Matrix newC0 = (Sso * C0prime);  

  // build density matrix
  Matrix newDensInit = Matrix(7);
  for(int k=0; k < 7; k++){
    for(int j=0; j < 7; j++){
      for(int i=0; i <4; i++){ 
        newDensInit(k,j) += newC0(k,i) * newC0(j,i);    
      }
    }
  }


  // Step 9: Compute New SCF Energy
  //
  double Eelec1 = 0;
  for(int i=0; i<7; i++){
    for(int j=0; j<7; j++){
      Eelec1 +=  newDensInit(i,j) * (Hcore(i,j) + newFockinit(i,j));
    }
  }
  
  Etot = Eelec1 + enuc;  //incorret
 
  }while(fabs(Etot - Etotprev)> 1e-6);  // Step 10: Test for Convergence




 return 0;

}

