#include <vector>
#include <cmath>
#include <stdio.h>
// possibly take in a vector of values
// use vector size to determine matrix size


class Matrix {
protected:
  int mDim;
  std::vector<std::vector<double> > mY;

public:

  Matrix(double elements[], int size) : mDim(size){
    mY.resize(mDim);
    for(int i=0; i < mY.size(); i++){
      mY[i].resize(mDim);
    }
    for(int i=0; i < mDim; i++){
      for(int j=0; j < mDim; j++){
        mY[i][j] = elements[(i*mDim)+j];
      }
    }
  }
  
  Matrix(int size) : mDim(size){
   mY.resize(mDim);
    for(int i=0; i < mY.size(); i++){
      mY[i].resize(mDim);
    }
    for(int i=0; i < mDim; i++){
      for(int j=0; j < mDim; j++){
        mY[i][j] = 0;
      }
    }
  }

public:

  double & operator()(const int & row, const int & col){
    return this->mY[row][col];
  }
  const double & operator()(const int & row, const int & col) const{
    return this->mY[row][col];
  }

  int getDim() const{return this->mDim;}
 
  double * ToArraySTATIC(){
    static double array[49];
      for(int i=0; i < mDim; i++){
        for(int j=0; j < mDim; j++){
          array[(i*mDim)+j] = this->mY[i][j];
        }
      }
    return array;
  }

  Matrix & operator=(const Matrix & rhs){
    if (&rhs == this)
      return *this;
    mDim = rhs.getDim();

    mY.resize(mDim);
    for(int i=0; i < mY.size(); i++){
      mY[i].resize(mDim);
    }
    
    for (int i = 0; i < mDim; i++){
      for(int j = 0; j < mDim; j++){
        mY[i][j] = rhs(i,j);
      }
    }
    return *this;
  }
  
  Matrix operator*(const Matrix & rhs){
    Matrix temp = Matrix((*this).mDim);
    // double num;
     for(int i = 0; i < mDim; i++){
      for(int k = 0; k < mDim; k++){
        //  num = 0;
        for(int j = 0; j < mDim; j++){
          temp(i,k) += (this->mY[i][j] * rhs.mY[j][k]);
        }
        //    temp.setX(i,k,num);
      }
     }
     return temp;

  }

  Matrix & operator*=(const Matrix & rhs){
    Matrix temp = (*this)*rhs;
    (*this) = temp;
    return *this;
  }

  Matrix T() { 
    Matrix temp = Matrix((*this).mDim);
    for(int i = 0; i < mDim; i++){
      for(int j = 0; j < mDim; j++){
        temp(i,j) = this->mY[j][i];
      }
    }
    return temp;
} 
         
    
};

  
class DiagMat : public Matrix {

public:
  
  DiagMat(double elements[], int size) : Matrix(size) {
    for (int i = 0; i < mDim; i++){
      mY[i][i] = elements[i];
    }
  }
  DiagMat(int size) : Matrix(size){}

  DiagMat ToPower(double pwr){
    DiagMat temp = DiagMat((*this).mDim);
    for(int i = 0; i < mDim; i++){
      temp(i,i) = pow(this->mY[i][i],pwr);
    }
    return temp;
  }

};

      
