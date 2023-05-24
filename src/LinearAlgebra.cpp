#include "../include/LinearAlgebra.h"

#include <cmath>
#include <stdexcept>

Vector::Vector() : n(0) {
  v = NULL;
}

Vector::Vector(int n) : n(n) {
  v = new double[n];
  for (int i = 0; i < n; i++)
    v[i] = 0.;
}

Vector::Vector(double x, double y, double z) : n(3) {
  v = new double[3];
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

Vector::Vector(const Vector& v_) : n(v_.n) {
  v = new double[n];
  for (int i = 0; i < n; i++) {
    v[i] = v_.v[i];
  }
}

Vector::~Vector() {
  delete[] v;
}

void Vector::print() const {
  for (int i = 0; i < n; i++)
    printf("%f ", v[i]);
  printf("\n");
}

int Vector::dim() const {
  return n;
}

double Vector::norm() const {
  double norm = 0.;
  for (int i = 0; i < n; i++)
    norm += v[i] * v[i];
  return sqrt(norm);
}

Vector Vector::normalized() const {
  double norm = this->norm();
  Vector v_(n);
  for (int i = 0; i < n; i++)
    v_.v[i] = v[i] / norm;
  return v_;
}

double Vector::dot(const Vector& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("Vector::dot: dimension mismatch");
  double dot = 0.;
  for (int i = 0; i < n; i++)
    dot += v[i] * v_.v[i];
  return dot;
}

Vector Vector::cross(const Vector& v_) const {
  if (n != 3 || v_.dim() != 3)
    throw std::invalid_argument("Vector::cross: dimension mismatch");
  Vector u(3);
  u.v[0] = v[1] * v_(2) - v[2] * v_(1);
  u.v[1] = v[2] * v_(0) - v[0] * v_(2);
  u.v[2] = v[0] * v_(1) - v[1] * v_(0);
  return u;
}

double Vector::operator()(int i) const {
  if (i < 0 || i >= n)
    throw std::out_of_range("Vector::operator(): index out of range");
  return v[i];
}

double& Vector::operator()(int i) {
  if (i < 0 || i >= n)
    throw std::out_of_range("Vector::operator(): index out of range");
  return v[i];
}

Vector& Vector::operator=(const Vector v_) {
  if (n != v_.n) {
    delete[] v;
    n = v_.n;
    v = new double[n];
  }
  for (int i = 0; i < n; i++)
    v[i] = v_.v[i];
  return *this;
}

Vector Vector::operator+(const Vector& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("Vector::operator+: dimension mismatch");
  Vector u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] + v_.v[i];
  return u;
}

Vector Vector::operator+=(const Vector& v_) {
  if (n != v_.n)
    throw std::invalid_argument("Vector::operator+: dimension mismatch");
  for (int i = 0; i < n; i++)
    v[i] += v_.v[i];
  return *this;
}
Vector Vector::operator-(const Vector& v_) const {
  if (n != v_.n)
    throw std::invalid_argument("Vector::operator-: dimension mismatch");
  Vector u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] - v_.v[i];
  return u;
}

Vector Vector::operator*(double a) const {
  Vector u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] * a;
  return u;
}

Vector Vector::operator/(double a) const {
  Vector u(n);
  for (int i = 0; i < n; i++)
    u.v[i] = v[i] / a;
  return u;
}

bool Vector::operator==(const Vector& v_) const {
  if (n != v_.n)
    return false;
  for (int i = 0; i < n; i++) {
    if (v[i] != v_.v[i])
      return false;
  }
  return true;
}

bool Vector::operator!=(const Vector& v_) const {
  return !(*this == v_);
}

std::ostream& operator<<(std::ostream& os, const Vector& u) {
  os << "(" << u(0) << ", " << u(1) << ", " << u(2) << ")";
  return os;
}

Matrix::Matrix() : m(0), n(0) {
  M = NULL;
}

Matrix::Matrix(int m, int n) : m(m), n(n) {
  // Allocate memory
  M = new double*[m];
  for (int i = 0; i < m; i++)
    M[i] = new double[n];

  // Initialization
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      M[i][j] = 0.;
  }
}

Matrix::Matrix(const double* p, int m, int n) : m(m), n(n) {
  // Allocate memory
  M = new double*[m];
  for (int i = 0; i < m; i++)
    M[i] = new double[n];

  // Initialization
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      M[i][j] = p[i * n + j];
  }
}

Matrix::Matrix(const Matrix& M_) : m(M_.m), n(M_.n) {
  // Allocate memory
  M = new double*[m];
  for (int i = 0; i < m; i++)
    M[i] = new double[n];

  // Initialization
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      M[i][j] = M_.M[i][j];
  }
}

Matrix::Matrix(int i, double angle) : m(3), n(3) {
  // Allocate memory
  M = new double*[m];
  for (int i = 0; i < m; i++)
    M[i] = new double[n];

  // Initialization
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      M[i][j] = 0.;
  }

  // Rotation matrix around axis i
  if (i == 1) {
    M[0][0] = 1.;
    M[1][1] = cos(angle);
    M[1][2] = sin(angle);
    M[2][1] = -sin(angle);
    M[2][2] = cos(angle);
  } else if (i == 2) {
    M[0][0] = cos(angle);
    M[0][2] = -sin(angle);
    M[1][1] = 1.;
    M[2][0] = sin(angle);
    M[2][2] = cos(angle);
  } else if (i == 3) {
    M[0][0] = cos(angle);
    M[0][1] = sin(angle);
    M[1][0] = -sin(angle);
    M[1][1] = cos(angle);
    M[2][2] = 1.;
  } else {
    throw std::invalid_argument("Invalid axis index");
  }
}

Matrix::~Matrix() {
  for (int i = 0; i < m; i++)
    delete[] M[i];
  delete[] M;
}

void Matrix::print() const {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      printf("%f ", M[i][j]);
    printf("\n");
  }
}

int Matrix::nrows() const {
  return m;
}

int Matrix::ncols() const {
  return n;
}

Matrix Matrix::identity(int n) {
  Matrix I(n, n);
  for (int i = 0; i < n; i++)
    I.M[i][i] = 1.;
  return I;
}

Matrix Matrix::transpose() const {
  Matrix result(n, m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[j][i] = M[i][j];
  }
  return result;
}

double Matrix::operator()(int i, int j) const {
  if (i < 0 || i >= m || j < 0 || j >= n)
    throw std::out_of_range("Matrix::operator(): index out of range");
  return M[i][j];
}

double& Matrix::operator()(int i, int j) {
  if (i < 0 || i >= m || j < 0 || j >= n)
    throw std::out_of_range("Matrix::operator(): index out of range");
  return M[i][j];
}

Matrix& Matrix::operator=(const Matrix& M_) {
  if (m != M_.m || n != M_.n) {
    for (int i = 0; i < m; i++)
      delete[] M[i];
    delete[] M;
    m = M_.m;
    n = M_.n;
    M = new double*[m];
    for (int i = 0; i < m; i++)
      M[i] = new double[n];
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      M[i][j] = M_.M[i][j];
  }
  return *this;
}

Matrix Matrix::operator+(const Matrix& M_) const {
  if (m != M_.m || n != M_.n)
    throw std::invalid_argument("Matrix::operator+: incompatible dimensions");

  Matrix result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[i][j] = M[i][j] + M_.M[i][j];
  }
  return result;
}

Matrix Matrix::operator-(const Matrix& M_) const {
  if (m != M_.m || n != M_.n)
    throw std::invalid_argument("Matrix::operator-: incompatible dimensions");

  Matrix result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[i][j] = M[i][j] - M_.M[i][j];
  }
  return result;
}

Matrix Matrix::operator*(double a) const {
  Matrix result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[i][j] = M[i][j] * a;
  }
  return result;
}

Matrix Matrix::operator/(double a) const {
  Matrix result(m, n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result.M[i][j] = M[i][j] / a;
  }
  return result;
}

Matrix Matrix::operator*(const Matrix& M_) const {
  if (n != M_.m)
    throw std::invalid_argument("Matrix::operator*: incompatible dimensions");

  Matrix result(m, M_.n);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < M_.n; j++) {
      for (int k = 0; k < n; k++)
        result.M[i][j] += M[i][k] * M_.M[k][j];
    }
  }
  return result;
}

Vector Matrix::operator*(const Vector& v_) const {
  if (n != v_.dim())
    throw std::invalid_argument("Matrix::operator*: incompatible dimensions");

  Vector result(m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++)
      result(i) += M[i][j] * v_(j);
  }
  return result;
}

std::ostream& operator<<(std::ostream& os, const Matrix& M) {
  for (int i = 0; i < M.m; i++) {
    for (int j = 0; j < M.n; j++)
      os << M.M[i][j] << " ";
    os << std::endl;
  }
  return os;
}