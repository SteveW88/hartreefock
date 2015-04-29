#include "matrices.hpp"



//////////////////////
/// Array Indexing ///
//////////////////////

int INDEX(int i, int j);


/////////////////////
/// Read in Files ///
/////////////////////

double ReadIN(const char filename[]);

void ReadIN(const char filename[], double out[]);

void ReadIN(const char filename[], Matrix & out);


///////////////////
/// Build Fock  ///
///////////////////

Matrix BuildFock(Matrix coreH, Matrix DensPrev, double tei[]);


/////////////////
/// Build Sso ///
/////////////////

Matrix BuildSso(double eigenval[], double eigenvec[], const Matrix & coreH);


////////////////
/// Build C0 ///
////////////////

Matrix BuildC0(double eigenval[], double eigenvec[], Matrix Sso, Matrix  Fock);

  
/////////////////////
/// Build Density ///
/////////////////////

Matrix BuildDensity(const Matrix & C0);


///////////////////
/// Compute SCF ///
///////////////////

double ComputeSCF(const Matrix & Density, const Matrix & coreH, const Matrix & Fock);


///////////////////
/// Diagonalize ///
///////////////////

void Diagonalize(double eigenval[], double eigenvec[], Matrix origMat);
