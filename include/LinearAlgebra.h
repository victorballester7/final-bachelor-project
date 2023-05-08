#ifndef LINEARALGEBRA_H
#define LINEARALGEBRA_H

#include <iostream>

class Vector {
 public:
  // Constructors
  Vector();                              // Null vector (0-dim)
  Vector(int n);                         // Zero vector (n-dim)
  Vector(double x, double y, double z);  // 3D vector (x,y,z)
  Vector(const Vector& v_);              // Copy constructor

  // Destructors
  ~Vector();

  // Note: const member functions cannot modify the object
  void print() const;                    // Print vector to stdout
  int dim() const;                       // Return dimension of the vector
  double norm() const;                   // Return norm of the vector
  Vector normalized() const;             // Return normalized vector
  double dot(const Vector& v_) const;    // Return dot product with v_
  Vector cross(const Vector& v_) const;  // Return cross product with v_

  double operator()(int i) const;                                        // Access element i (const), i = 0, ..., n-1
  double& operator()(int i);                                             // Access element i (non-const), i = 0, ..., n-1
  Vector& operator=(const Vector v_);                                    // Assignment
  Vector operator+(const Vector& v_) const;                              // Vector addition
  Vector operator-(const Vector& v_) const;                              // Vector subtraction
  Vector operator*(double a) const;                                      // Scalar multiplication from the left
  Vector operator/(double a) const;                                      // Scalar division
  bool operator==(const Vector& v_) const;                               // Vector equality
  bool operator!=(const Vector& v_) const;                               // Vector inequality
  friend std::ostream& operator<<(std::ostream& os, const Vector& Vec);  // Print vector to output stream

 private:
  int n;      // Dimension
  double* v;  // Vector v_(n)
};

class Matrix {
 public:
  // Constructors
  Matrix();                                // Null matrix (0 x 0)
  Matrix(int m, int n);                    // Zero m x n matrix
  Matrix(const double* p, int m, int n);   // m x n matrix with elements p
  Matrix(const double** p, int m, int n);  // m x n matrix with elements p
  Matrix(const Matrix& M_);                // Copy constructor
  Matrix(int i, double angle);             // rotation matrix around axis i (i = 1=x, 2=y, 3=z)

  // Destructor
  ~Matrix();

  // Assignment
  void print() const;        // Print matrix to stdout
  int nrows() const;         // Return number of rows
  int ncols() const;         // Return number of columns
  Matrix identity(int n);    // Return n x n identity matrix
  Matrix transpose() const;  // Return transpose of the matrix
  // Matrix Inverse() const; // Return inverse of the matrix
  // double Determinant() const; // Return determinant of the matrix

  double operator()(int i, int j) const;  // Access element (i,j) (const)
  double& operator()(int i, int j);       // Access element (i,j) (non-const)

  Matrix& operator=(const Matrix& M_);  // Assignment

  Matrix operator+(const Matrix& M) const;  // Matrix addition
  Matrix operator-(const Matrix& M) const;  // Matrix subtraction
  Matrix operator*(double a) const;         // Scalar multiplication
  Matrix operator/(double a) const;         // Scalar division
  bool operator==(const Matrix& M) const;   // Matrix equality
  bool operator!=(const Matrix& M) const;   // Matrix inequality

  Matrix operator*(const Matrix& M_) const;  // Matrix multiplication
  Vector operator*(const Vector& v) const;   // Matrix-vector multiplication

 private:
  int m;       // First dimension (number of rows)
  int n;       // Second dimension (number of columns)
  double** M;  // Matrix M(m,n)
};

#endif  // LINEARALGEBRA_H