#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "matrices.hpp"

#define BIGNUM 10000

//////////////////////
/// Array Indexing ///
//////////////////////

int ioff[BIGNUM]
ioff[0] = 0;
for(int i=1; i < BIGNUM; i++)
  ioff[i] = ioff[i-1] + i;

int INDEX(int i,int j){
  return (i>j) ? (ioff[i]+j) : (ioff[j]+i);
}

/////////////////////
/// Read in Files ///
/////////////////////

double ReadIN(str filename){
  double val;
  std::ifstream File;
  File.open(filename,std::ifstream::in);
  File >> val;
  File.close();
  return val;
  
}

void ReadIN(str filename, double out[]){
  FILE *input;
  double i, j, val;
  double dim = out.size();
  input = fopen(filename,"r");
  while(fscanf(input,"%d %d %d %d %1f", &i,&j,&k,&l&val) != EOF){
    i -= 1;
    j -= 1;
    k -= 1;
    l -= 1;
    ij = INDEX(i,j);
    kl = INDEX(k,l);
    ijkl = INDEX(ij,kl);
    
    out[ijkl] = val;
  }
  fclose(input);
  
}

void ReadIN(str filename, Matrix & out){
  FILE *input;
  double i, j, val;
  double dim = out.getDim();
  input = fopen(filename,"r");
  while(fscanf(input,"%d %d %1f", &i,&j,&val) != EOF){
    i -= 1;
    j -= 1;
    out(i,j) = val;
    if(i != j)
      out(j,i) = val;
  }
  fclose(input);
  
}

///////////////////
/// Build Fock  ///
///////////////////

Matrix BuildFock(Matrix coreH, Matrix DensPrev, double tei[]){
  Matrix Fock = coreH;
  int ij, kl, ik, jl, ijkl, ikjl;
  for(int i=0; i < dim; i++){
    for(int j=0; j < dim; j++){
      for(int k=0; k < dim; k++){
        for(int l=0; l < dim; l++){
          ij = INDEX(i,j);
          kl = INDEX(k,l);
          ik = INDEX(i,k);
          jl = INDEX(j,l);
          ijkl = INDEX(ij,kl);
          ikjl = INDEX(ik,jl);
          Fock(i,j) += DensPrev(k,l) * (2.0 *tei[ijkl] - tei[ikjl]);
        }
      }
    }
  }
  
  return Fock;
}


/////////////////
/// Build Sso ///
/////////////////

Matrix BuildSso(double eigenval[], double eigenvec[], Matrix coreH){
  double dim = coreH.getDim();
  DiagMat As = DiagMat(eigenval,dim);
  Matrix Ls = Matrix(eigenvec,dim).T();
  return (Ls * As.ToPower(-0.5) * Ls.T());

}

////////////////
/// Build C0 ///
////////////////

Matrix BuildC0(double eigenval[], double eigenvec[], Matrix Sso, Matrix Fock){
  int dim = Fock.getDim();
  Matrix Fprime = (Sso.T() * Fock * Sso);
  Diagonalize(eigenval,eigenvec,Fprime);
  Matrix C0prime = Matrix(eigenvec,dim);
  
  return (Sso * C0prime);

}

  
/////////////////////
/// Build Density ///
/////////////////////

Matrix BuildDensity(Matrix C0){
  double dim = C0.getDim();
  Matrix Density = Matrix(dim);
  for(int k=0; k < dim; k++){
    for(int j=0; j < dim; j++}{
      for(int i=0; i < ((dim+1)*0.5); i++){
        Density(k,j) += C0(k,i) * C0(j,i);
      }
    }
  }
  return Density;
  
  
}


///////////////////
/// Compute SCF ///
///////////////////
  
double ComputeSCF(Matrix Density, Matrix coreH, Matrix Fock){
  int dim = coreH.getDim();
  double Eelec = 0;
  for(int i=0; i < dim; i++){
    for(int j=0; j < dim; j++){
      Eelec += Density(i,j) * (coreH(i,j) + Fock(i,j));
    }
  }
  return Eelec;
}

///////////////////
/// Diagonalize ///
///////////////////

void Diagonlize(double eigenval[], double eigenvec[], Matrix origMat){
  int dim = origMat.getDim();
  double * temparray = origMat.ToArray();
  double array[dim*dim];
  for(int i=0; i < (dim*dim); i++){
    array[i] = *(temparray+i);
  }

  gsl_matrix_view m
    = gsl_matrix_view_array (origMat, dim, dim);

  gsl_vector_complex *eval = gsl_vector_complex_alloc (dim);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (dim, dim);

  gsl_eigen_nonsymmv_workspace * w =
    gsl_eigen_nonsymmv_alloc (dim);

  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec,
                           GSL_EIGEN_SORT_ABS_DESC);

  {
    int i, j;

    for (i = 0; i < dim; i++)
      {
        gsl_complex eval_i
          = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i
          = gsl_matrix_complex_column (evec, i);

        eigenval[i] = GSL_REAL(eval_i);

        for (j = 0; j < dim; ++j)
          {
            gsl_complex z =
              gsl_vector_complex_get(&evec_i.vector, j);

            eigenvec[(i*dim)+j] = GSL_REAL(z);
          }
      }
  }

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
}





