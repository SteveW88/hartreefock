#include "matrices.hpp"

int INDEX(int i, int j);

double ReadIN(const char filename[]);

void ReadIN(const char filename[], double out[]);

void ReadIN(const char filename[], Matrix & out);

Matrix BuildFock(Matrix coreH, Matrix DensPrev, double tei[]);

Matrix BuildSso(double eigenval[], double eigenvec[], const Matrix & coreH);

Matrix BuildC0(double eigenval[], double eigenvec[], Matrix Sso, Matrix  Fock);


Matrix BuildDensity(const Matrix & C0);

double ComputeSCF(const Matrix & Density, const Matrix & coreH, const Matrix & Fock);

void Diagonalize(double eigenval[], double eigenvec[], Matrix origMat);
