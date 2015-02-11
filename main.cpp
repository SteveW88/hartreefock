#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;


int main(int argc, char *argv[]) {
  
  // Array Indexing
  int BIGNUM = 100;
  int ioff[BIGNUM];
  ioff[0] = 0;
  for(int i = 1; i < BIGNUM; i++)
    ioff[i] = ioff[i-1] + i;


  int i, j, k, l, ij, kl, ijkl;
  // Row/Column to one dimensional array indexing transformation
  // ij = (i>j) ? ioff[i] + j : ioff[j] + i;
  // ijkl = (ij > kl) ? ioff[ij] + kl : ioff[kl] + ij; //TODO: Verify

  //Read enuc
  ifstream enucFile;
  enucFile.open("enuc.dat",ifstream::in);
  double enuc;
  enucFile >> enuc;
  enucFile.close();
 

  //Read overlap
  FILE *input;
  double val;
  double overlap[BIGNUM];
  input = fopen("s.dat", "r");
  while(fscanf(input, "%d %d %lf", &i, &j, &val) != EOF){
    i -= 1;
    j -= 1;
    ij = (i>j) ? ioff[i] + j : ioff[j] + i;
    overlap[ij] = val;
  }
  fclose(input);

  //Read kinetic energy
  double ke[BIGNUM];
  input = fopen("t.dat", "r");
  while(fscanf(input, "%d %d %lf", &i, &j, &val) != EOF){
    i -= 1;
    j -= 1;
    ij = (i>j) ? ioff[i] + j : ioff[j] + i;
    ke[ij] = val;
  }
  fclose(input);
  

  //Read nuclear-attraction integrals
  double nai[BIGNUM];
  input = fopen("v.dat", "r");
  while(fscanf(input, "%d %d %lf", &i, &j, &val) != EOF){
    i -= 1;
    j -= 1;
    ij = (i>j) ? ioff[i] + j : ioff[j] + i;
    nai[ij] = val;
  }
  fclose(input);

  //Core Hamiltonia
  double coreH[BIGNUM];
  
  
}


