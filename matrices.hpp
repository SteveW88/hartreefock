#include <vector>
#include <cmath>
#include <stdio.h>
// possibly take in a vector of values
// use vector size to determine matrix size


namespace R7{

class Matrix {
protected:
  int dim;
  double **mX;// = new int*[dim];
  
public:

  Matrix(std::vector<double> elements, int size) : dim(size){
    mX = new double*[dim];
    for(int l = 0; l < dim; l++)
      mX[l] = new double[dim];
    for (int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++)
        {
          mX[i][j] = elements.at((i*dim)+j);
        }
    }
         
  }
  Matrix(double elements[], int size) : dim(size){
    mX = new double*[dim];
    for(int l = 0; l < dim; l++)
      mX[l] = new double[dim];
    for (int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++)
        {
          mX[i][j] = elements[(i*dim)+j];
        }
    }
  }
  Matrix(int size) : dim(size){
    mX = new double*[dim];
    for(int l = 0; l < dim; l++)
      mX[l] = new double[dim];
    for (int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++)
        {
          mX[i][j] = 0;
        }
    }
  }
public:

  double X(int row, int col){return mX[row][col];}
  double setX(int row, int col, double val){mX[row][col] = val;}

  void out(){
    for (int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++)
        {
          printf("%f", **(mX[i][j]);
        }
  }
  }

  Matrix T() { 
    Matrix temp = Matrix((*this).dim);
    for(int i = 0; i < dim; i++){
      for(int j = 0; j < dim; j++){
        temp.setX(i,j,(*this).mX[j][i]);
      }
    }
    return temp;} 

  Matrix & operator*=(const Matrix & rhs){
    Matrix temp = Matrix((*this).dim);
    double num;
    for(int i = 0; i < dim; i++){
      for(int k = 0; k < dim; k++){
        num = 0;
        for(int j = 0; j < dim; j++){
          num += (mX[i][j] * rhs.mX[j][k]);
        }
        temp.setX(i,k,num);
      }
    }
    return temp;
  }

  //Matrix operator*(Matrix lhs, const Matrix & rhs){
  //  return lhs *= rhs;
  //}

};

inline Matrix operator*(Matrix lhs, const Matrix & rhs){
  return lhs *= rhs;
}
  
class DiagMat : public Matrix {

public:
  
  DiagMat(double elements[], int size) : Matrix(size) {
    for (int i = 0; i < dim; i++){
      setX(i,i,elements[i]);
    }
  }
  DiagMat(int size) : Matrix(size){}

  DiagMat ToPower(double pwr){
    DiagMat temp = DiagMat((*this).dim);
    for(int i = 0; i < dim; i++){
      temp.setX(i,i,pow((*this).X(i,i),pwr));
    }
    return temp;
  }

};

}               
