#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "matrices.hpp"
using namespace std;


int main(int argc, char *argv[]) {
  
  // Array Indexing
  int BIGNUM = 10000;
  int ioff[BIGNUM];
  ioff[0] = 0;
  for(int i = 1; i < BIGNUM; i++)
    ioff[i] = ioff[i-1] + i;


  int i, j, k, l, ij, kl, ijkl;
  // Row/Column to one dimensional array indexing transformation
  // ij = (i>j) ? ioff[i] + j : ioff[j] + i;
  // ijkl = (ij > kl) ? ioff[ij] + kl : ioff[kl] + ij; //TODO: Verify

  //Step #1: Read enuc
  //
  ifstream enucFile;
  enucFile.open("enuc.dat",ifstream::in);
  double enuc;
  enucFile >> enuc;
  enucFile.close();
 
  
  //Step #2: Read overlap
  //
  FILE *input;
  double val;
  double overlap[BIGNUM];
  double overlap2[7][7];
  input = fopen("s.dat", "r");
  while(fscanf(input, "%d %d %lf", &i, &j, &val) != EOF){
    i -= 1;
    j -= 1;
    ij = (i>j) ? ioff[i] + j : ioff[j] + i;
    overlap[ij] = val;
    overlap2[i][j] = val;
    if (i != j)
      overlap2[j][i] = val;
  }
  fclose(input);

  //Read kinetic energy
  //
  double ke[BIGNUM];
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
  
 gsl_matrix_view m 
   = gsl_matrix_view_array (coreH2, 7, 7);

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

        eigenvalues[i] = GSL_REAL(eval_i);
        for (j = 0; j < 7; ++j)
          {
            gsl_complex z = 
              gsl_vector_complex_get(&evec_i.vector, j);
            printf("%g \n", GSL_REAL(z));
            eigenvectors[(i*7)+j] = GSL_REAL(z);
          }
      }
  }
  
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
  
  R7::Matrix mat = R7::Matrix(coreH2);
  R7::DiagMat As = R7::DiagMat(eigenvalues);
  R7::Matrix Ls = R7::Matrix(eigenvectors).T();
  R7::DiagMat test = As.ToPower(-1/2);
  // R7::Matrix Snegsqrt = Ls * As.ToPower(-1/2) * Ls.T();
  //R7::Matrix Snegsqrt = Ls * Ls; //Ls.T();

  double tarray[9] = {1, 1, 0,
                      0, 2, 0,
                      0, -1,4};

gsl_matrix_view m1 
   = gsl_matrix_view_array (tarray, 3, 3);

 gsl_vector_complex *eval1 = gsl_vector_complex_alloc (3);
 gsl_matrix_complex *evec1 = gsl_matrix_complex_alloc (3, 3);

  gsl_eigen_nonsymmv_workspace * w1 = 
    gsl_eigen_nonsymmv_alloc (3);
  
  gsl_eigen_nonsymmv (&m1.matrix, eval1, evec1, w1);
 
  gsl_eigen_nonsymmv_free (w1);

  gsl_eigen_nonsymmv_sort (eval1, evec1, 
                           GSL_EIGEN_SORT_ABS_DESC);
 
  {
    int i, j;

    for (i = 0; i < 3; i++)
      {
        gsl_complex eval_i 
          = gsl_vector_complex_get (eval1, i);
        gsl_vector_complex_view evec_i 
          = gsl_matrix_complex_column (evec1, i);

        printf ("eigenvalue = %g \n",
                GSL_REAL(eval_i));
        printf ("eigenvector = \n");

        //eigenvalues[i] = GSL_REAL(eval_i);
        for (j = 0; j < 3; ++j)
          {
            gsl_complex z = 
              gsl_vector_complex_get(&evec_i.vector, j);
            printf("%g \n", GSL_REAL(z));
            //eigenvectors[(i*3)+j] = GSL_REAL(z);
          }
      }
  }
  
  gsl_vector_complex_free(eval1);
  gsl_matrix_complex_free(evec1);

  
  return 0;

}
