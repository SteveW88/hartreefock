#include <vector>
#include <cmath>
#include <stdio.h>


/////////////////////
/// Class: Matrix ///
/////////////////////
// N x N matrix class

class Matrix {
protected:
  int mDim;                                 // dimension of matrix
  std::vector<std::vector<double> > mY;     // double vector of values
  std::vector<double> mArray;               // single vector of valeus
                                            // (used for ToArray function)

public:

  ////////////////////
  /// Constructors ///
  ////////////////////

  // Given array of elements
  Matrix(const double elements[], int size) : mDim(size){
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

  // Matrix of zeros
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

  ///////////////////////
  /// Class Functions ///
  ///////////////////////

  // Return matrix element i,j by reference
  double & operator()(const int & row, const int & col){
    return this->mY[row][col];
  }

  // Return matrix element i,j by constant reference
  const double & operator()(const int & row, const int & col) const{
    return this->mY[row][col];
  }

  // Return dimension of matrix
  int getDim() const{return this->mDim;}
 
  // Return pointer to mArray:
  //   (needed for GSL diagonalization procedure)
  double * ToArray(){
    for(int i=0; i < mDim; i++){
      for(int j=0; j < mDim; j++){
        mArray.push_back(this->mY[i][j]);
      }
    }
    double * a = &mArray[0];
    return a;
  }


  // Matrix assignment
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
  
  // Matrix-Matrix multiplication
  Matrix operator*(const Matrix & rhs){
    Matrix temp = Matrix((*this).mDim);
    // double num;
     for(int i = 0; i < mDim; i++){
      for(int k = 0; k < mDim; k++){
        //  num = 0;
        for(int j = 0; j < mDim; j++){
          temp(i,k) += (this->mY[i][j] * rhs.mY[j][k]);
        }
      }
     }
     return temp;

  }

  // Matrix-Matrix multiplication
  Matrix & operator*=(const Matrix & rhs){
    Matrix temp = (*this)*rhs;
    (*this) = temp;
    return *this;
  }


  // Matrix transpose
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

  
/////////////////////
/// Class DiagMat ///
/////////////////////
// Diagonal n x n matrix class

class DiagMat : public Matrix {

public:
  
  ////////////////////
  /// Constructors ///
  ////////////////////

  // Given diagonal elements as array
  DiagMat(const double elements[], int size) : Matrix(size) {
    for (int i = 0; i < mDim; i++){
      mY[i][i] = elements[i];
    }
  }

  // N x N matrix of zeros 
  DiagMat(int size) : Matrix(size){}

public:

  /////////////////
  /// Functions ///
  /////////////////

  // Raise diagonal elements to a given power
  DiagMat ToPower(double pwr){
    DiagMat temp = DiagMat((*this).mDim);
    for(int i = 0; i < mDim; i++){
      temp(i,i) = pow(this->mY[i][i],pwr);
    }
    return temp;
  }

};

      
