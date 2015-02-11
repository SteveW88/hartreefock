#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

class Matrix{

  std::vector<double> mElements; 

public:
  Matrix(std::vector<double> elements):
    mElements(elements){}


};





int main(int argc, char *argv[]) {
 

  ifstream enuc;
  enuc.open("enuc.dat");
  char data[100];

  enuc >> data;

  cout << data << endl;

  std::vector<double> vecta;

  vecta.push_back(1.0);
  vecta.push_back(2.0);
  vecta.push_back(3.0);
  vecta.push_back(4.0);

  Matrix a = Matrix(vecta);
  int sup = 1;
}


