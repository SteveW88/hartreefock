#include <vector>
#include <cmath>
// possibly take in a vector of values
// use vector size to determine matrix size

namespace R7{

class Matrix {
protected:
  double mX[7][7];

public:

  Matrix(std::vector<double> elements){
    for (int i = 0; i < 7; i++){
      for(int j = 0; j < 7; j++)
        {
          mX[i][j] = elements.at((i*7)+j);
        }
    }
         
  }
  Matrix(double elements[]){
    for (int i = 0; i < 7; i++){
      for(int j = 0; j < 7; j++)
        {
          mX[i][j] = elements[(i*7)+j];
        }
    }
  }
  Matrix(){
    for (int i = 0; i < 7; i++){
      for(int j = 0; j < 7; j++)
        {
          mX[i][j] = 0;
        }
    }
  }
public:

  double X(int row, int col){return mX[row][col];}
  double setX(int row, int col, double val){mX[row][col] = val;}

  Matrix T() { 
    Matrix temp;
    for(int i = 0; i < 7; i++){
      for(int j = 0; j < 7; j++){
        temp.setX(i,j,(*this).mX[j][i]);
      }
    }
    return temp;} 

  Matrix & operator*=(const Matrix & rhs){
    Matrix temp;
    double num;
    for(int i = 0; i < 7; i++){
      for(int k = 0; k < 7; k++){
        num = 0;
        for(int j = 0; j < 7; j++){
          num += mX[i][j] * rhs.mX[j][k];
        }
        temp.setX(i,k,num);
      }
    }
  }

  // Matrix operator*(Matrix lhs, const Matrix & rhs){
  //    return lhs *= rhs;
  //  }

};

  
class DiagMat : Matrix {

public:
  
  DiagMat(double elements[]) : Matrix() {
    for (int i = 0; i < 7; i++){
      setX(i,i,elements[i]);
    }
  }
  DiagMat() : Matrix(){}

  DiagMat ToPower(double pwr){
    DiagMat temp;
    for(int i = 0; i < 7; i++){
      temp.setX(i,i,pow((*this).X(i,i),pwr));
    }
    return temp;
  }

};

}               
