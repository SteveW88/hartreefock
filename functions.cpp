#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "functions.hpp"
#include "global.hpp"


//////////////////////
/// Array Indexing ///
//////////////////////
// ioff array populated in main.cpp.
// This precomputed array reduces computation time
// of indexing.

int INDEX(int i,int j){
  return (i>j) ? (ioff[i]+j) : (ioff[j]+i);
}

/////////////////////
/// Read in Files ///
/////////////////////
// Three read in methods defined:
//
// double(filename): designed for file containing
//                   one value. (Typically nuclear repulsion energy)
//
// void(filename,array) : designed for file containing two-electron
//                        repulsion integrals.
//
// void(filename,matrix) : designed for file containing one-electron
//                         integrals.                        
// 

double ReadIN(const char filename[]){
  double val;
  std::ifstream File;
  File.open(filename,std::ifstream::in);
  File >> val;
  File.close();
  return val;
  
}

void ReadIN(const char filename[], double out[]){
  FILE *input;
  int i, j, k, l, ij, kl, ijkl;
  double val;
  input = fopen(filename,"r");
  while(fscanf(input,"%d %d %d %d %lf", &i,&j,&k,&l,&val) != EOF){
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

void ReadIN(const char filename[], Matrix & out){
  FILE *input;
  int i, j;
  double val;
  double dim = out.getDim();
  input = fopen(filename,"r");
  while(fscanf(input,"%d %d %lf", &i,&j,&val) != EOF){
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
// Constructs Fock matrix using core hamiltonian,
// density matrix, and two-electron repulsion integrals.
//

Matrix BuildFock(Matrix coreH, Matrix DensPrev, double tei[]){
  Matrix Fock = coreH;
  int dim = coreH.getDim();
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
// Constructs symmetric orthogonalization matrix. 
//

Matrix BuildSso(double eigenval[], double eigenvec[], const Matrix & coreH){
  double dim = coreH.getDim();
  DiagMat As = DiagMat(eigenval,dim);
  Matrix Ls = Matrix(eigenvec,dim).T();
  return (Ls * As.ToPower(-0.5) * Ls.T());

}

////////////////
/// Build C0 ///
////////////////
// Transforms eigenvectors of new fock matrix
// into the original AO basis.
//

Matrix BuildC0(double eigenval[], double eigenvec[], Matrix Sso, Matrix Fock){
  int dim = Fock.getDim();
  Matrix Fprime = (Sso.T() * Fock * Sso);
  Diagonalize(eigenval,eigenvec,Fprime);
  Matrix C0prime = Matrix(eigenvec,dim).T();
  Matrix C0 = Sso * C0prime;

  return C0;

}

  
/////////////////////
/// Build Density ///
/////////////////////
// Construct density matrix using occupied MOs.
//

Matrix BuildDensity(const Matrix & C0){
  double dim = C0.getDim();
  Matrix Density = Matrix(dim);
  for(int k=0; k < dim; k++){
    for(int j=0; j < dim; j++){
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
// Compute SCF electronic energy using denstiy matrix,
// core hamiltonian, and fock matrix.
//
  
double ComputeSCF(const Matrix & Density, const Matrix & coreH, const Matrix & Fock){
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
// GSL procedure to diagonalize a given matrix.
//

void Diagonalize(double eigenval[], double eigenvec[], Matrix origMat){
  int dim = origMat.getDim();
  double * temparray = origMat.ToArray();
  
  gsl_matrix_view m
    = gsl_matrix_view_array (temparray, dim, dim);

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





