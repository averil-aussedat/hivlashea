#include <iostream>
#include "Vector.h"
#ifndef MATRIX_H
#define MATRIX_H
class Matrix
{
 private:
  int length;
  int width;
  Vector* Coef;

 public:
  Matrix(int=0,int=0,double** = NULL);
  Matrix(const Matrix&);
  ~Matrix();
  //Accesseurs
  int getLength()const;
  int getWidth()const;
  //Operation Elementaire
  Matrix& operator=(const Matrix&);
  Vector& operator[](int);
  Matrix operator+(const Matrix&)const;
  Matrix operator-(const Matrix&)const;
  Matrix operator*(const Matrix&)const;
  Matrix transp()const;
  double Max()const;
  friend Matrix operator*(double,const Matrix &);
  friend Vector operator*(const Matrix&,const Vector &);
  friend int Gauss(Matrix&,Vector&); //Algorithme de decomposition LU
  friend void ConjuguateGradient(Matrix&,Vector&,Vector&,Matrix&,int,double); // Algorithme Gradient Conjugue avec Preconditionnement de matrice
  friend Vector Solve(Matrix&,Vector&,Vector&);
  friend Vector TestGauss(Matrix&,Vector&);
  friend std::ostream& operator<<(std::ostream&,const Matrix&);
};
#endif
