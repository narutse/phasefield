#ifndef MATRIX_HPP_
#define MATRIX_HPP_

// My simple matrix class.
// Author: Narut Sereewattanawoot

#include <vector>

// Matrix Class
template <class T>
class Matrix {
private:
  using size_t = std::size_t;
  size_t m, n;
  std::vector<T> M;

public:
  Matrix() {}
  Matrix(size_t m, size_t n) : m(m), n(n), M(m*n) {}
  Matrix(size_t m, size_t n, const T& val) : m(m), n(n), M(m*n, val) {}
  Matrix(const Matrix<T> &x) : m(x.m), n(x.n), M(x.M) {}
  Matrix(Matrix<T>&& x) : m(x.m), n(x.n), M(std::move(x.M)) {}

  const T& operator() (size_t i, size_t j) const { return M[i*n+j]; }
  T& operator() (size_t i, size_t j) { return M[i*n + j]; }

  const T& operator[] (size_t i) const { return M[i]; }
  T& operator[] (size_t i) { return M[i]; }  

  const T* data() const { return M.data(); }
  T* data() { return M.data(); }

  size_t Sub2Ind(size_t i, size_t j) const { return i*n+j; }
  size_t GetM() const { return m; }
  size_t GetN() const { return n; }
  size_t size() const { return M.size(); }  
  void swap(Matrix<T>& x);

  Matrix<T>& operator= (const Matrix<T>& x);
  Matrix<T>& operator= (Matrix<T>&& x);
};

// Matrix Implementation:
template <class T>
void Matrix<T>::swap(Matrix<T>& x) {
  std::swap(m, x.m);
  std::swap(n, x.n);
  M.swap(x.M);
}

template <class T>
Matrix<T>& Matrix<T>::operator= (const Matrix<T>& x) {
  m = x.m;
  n = x.n;
  M = x.M;
  return *this;
}

template <class T>
Matrix<T>& Matrix<T>::operator= (Matrix<T>&& x) {
  m = x.m;
  n = x.n;
  M = std::move(x.M);
  return *this;
}

// SubMatrix Class
template <class T>
struct SubMatrix {
private:
  size_t m, n, n_orig; 
  T *M;

public:
  SubMatrix() {}  
  SubMatrix(Matrix<T> &M, size_t start_i, size_t start_j,
    size_t m, size_t n);  
  SubMatrix(SubMatrix<T> &M, size_t start_i, size_t start_j,
    size_t m, size_t n);  
  SubMatrix(T *M, size_t m, size_t n);
  SubMatrix(const SubMatrix<T> &x);

  size_t GetM() const { return m; }
  size_t GetN() const { return n; }  
  const T* data() const { return M; }  
  T* data() { return M; }
  
  SubMatrix<T>& operator= (const SubMatrix<T>& x);    

  const T& operator() (size_t i, size_t j) const { return M[i*n_orig + j]; }
  T& operator() (size_t i, size_t j) { return M[i*n_orig + j]; }

  template <class U>
  void CopyFrom(const U& src);

  template <class U>
  void CopyFromArray(const U& src);

  template <class U>
  void Fill(U& dest);

  template <class U>
  void FillArray(U&& dest);
};

// Submatrix Implementation
template <class T>
SubMatrix<T>::SubMatrix(Matrix<T> &M, size_t start_i, size_t start_j,
  size_t m, size_t n) :
  m(m), n(n), n_orig(M.GetN()), M(M.data() + (start_i * n_orig + start_j))
{}

template <class T>
SubMatrix<T>::SubMatrix(SubMatrix<T> &M, size_t start_i,
  size_t start_j, size_t m, size_t n) :
  m(m), n(n), n_orig(M.n_orig), M(M.M + (start_i * n_orig + start_j))
  {}

template <class T>
SubMatrix<T>::SubMatrix(T *M, size_t m, size_t n) :
  m(m), n(n), n_orig(n), M(M)
{}

template <class T>
SubMatrix<T>::SubMatrix(const SubMatrix<T>& x) :
  m(x.m), n(x.n), n_orig(x.n_orig), M(x.M) {}

template <class T>
SubMatrix<T>& SubMatrix<T>::operator=(const SubMatrix<T>& x) {
  this->m = x.m;
  this->n = x.n;
  this->n_orig = x.n_orig;
  this->M = x.M;
  return *this;
}

template <class T> template <class U>
void SubMatrix<T>::CopyFrom(const U& src) {
  for (size_t i = 0; i < m; ++i) {
    std::copy(&src(i, 0), &src(i, n), &M[i*n_orig]);
  }
}

template <class T> template <class U>
void SubMatrix<T>::CopyFromArray(const U& src) {
  for (size_t i = 0; i < m; ++i) {
    std::copy(&src[i*n], &src[(i+1)*n], &M[i*n_orig]);
  }
}

template <class T> template <class U>
void SubMatrix<T>::Fill(U& dest) {
  for (size_t i = 0; i < m; ++i) {
    std::copy(&M[i*n_orig], &M[i*n_orig + n], &dest(i, 0));
  }
}

template <class T> template <class U>
void SubMatrix<T>::FillArray(U&& dest) {
  for (size_t i = 0; i < m; ++i) {
    std::copy(&M[i*n_orig], &M[i*n_orig + n], &dest[i*n]);
  }
}
// template <class T> template <class U>
// void SubMatrix<T>::FillArray(U *dest) {
//   for (size_t i = 0; i < m; ++i) {
//     std::copy(&M[i*n_orig], &M[i*n_orig + n], )
//   }
// }

#endif
