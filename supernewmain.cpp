
#define dim 7

int main(int argc, char *argv[]){

  double enuc = ReadIN("enuc.dat");
  double overlap[dim][dim], ke[dim][dim], nai[dim][dim];
  ReadIN("s.dat",overlap);
  ReadIN("t.dat",ke);
  ReadIN("v.dat",nai);

  double coreH[dim*dim];
  for(int i=0; i < dim; i++){
    for(int j=0; j < dim; j++){
      coreH[(i*dim)+j] = ke[i][j] + nai[i][j];
    }
  }

  double tei[ size ];
  ReadIN("eri.dat",tei);  //needs to be configured

  double eigenvalues[dim];
  double eigenvectors[dim*dim];
  
  Diagonalize(eigenvalues,eigenvectors,overlap);

}
