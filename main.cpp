#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <vector>
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
  double coreH[7][7];
  for(int i = 0; i < 7; i++){
    for(int j = 0; j < 7; j++){
      coreH[i][j] = ke2[i][j] + nai2[i][j];
    }
  }
  
 
  //Read two-electron integrals
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
  //Diagonalize overlap
  ///////////////////////
  //Use library for this??

  



}


