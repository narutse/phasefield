#ifndef __FBESOLVER_HPP__
#define __FBESOLVER_HPP__

#include <iostream>
#include <Eigen/Sparse>
#include "MMSfunc.hpp"
#include "matrix.hpp"

// #define PRINT_ITER
namespace Solver
{
  using namespace MMSfunc;
  class FBESolver
  {
    template <class T> using vector = std::vector<T>;    
    using SpMat = Eigen::SparseMatrix<double>;
    using Vec = Eigen::VectorXd;
    using Map = Eigen::Map<Vec>;
    // Boundary locations. These never change:
    static constexpr double XL = 0.0;
    static constexpr double XR = 1.0;
    static constexpr double YB = 0.0;
    static constexpr double YT = 0.5;
    static constexpr double default_tf = 8.0;

    // member variables:
    const Param &param;
    int M;
    int N;
    double dt;
    double dx;
    double dy;

    vector<double> x;
    vector<double> y;
    Vec sol;
    Vec rhs;
    SpMat mat;

    void ConstructMatrix();
    void FillRHS(double t);
    void Integrate(double t_final);
    int sub2ind(int i, int j) { return i*N+j; };

  public:
    FBESolver(const Param &param, int M, int N, double dt)
    : param(param), M(M), N(N), dt(dt),
      dx(1.0 * (XR - XL) / (M - 1)), dy(1.0 * (YT - YB) / (N - 1)),
      x(M), y(N),
      sol((M-1)*N), rhs((M-1)*N), mat((M-1)*N,(M-1)*N)
    {
      for (int i = 0; i < M; ++i) {
        x[i] = XL + (XR-XL)*1.0*i/(M-1);
      }
      for (int j = 0; j < N; ++j) {
        y[j] = YB + (YT-YB)*1.0*j/(N-1);
      }
    }

    Matrix<double> Solve(double t_final = default_tf) 
    {
      Integrate(t_final);
      // put output in a more readable form:
      Matrix<double> result(M, N);
      for (int i = 0; i < M-1; ++i) {
        for (int j = 0; j < N; ++j) {
          result(i, j) = sol[sub2ind(i, j)];
        }
      }
      // copy grid point not used because of periodic BCs
      for (int j = 0; j < N; ++j) result(M-1, j) = result(0, j);
      return result;
    }

    Matrix<double> GetAnalytic(double t);
    double GetDx() { return dx; }
    double GetDy() { return dy; }
  };

  // implementations of non-trivial functions here:
  void FBESolver::ConstructMatrix()
  {
    int n = (M-1)*N;
    double diag = 1 + 4*dt*param.kappa/dx/dy;
    double offdiag = -dt*param.kappa/dx/dy;
    mat.reserve(Vec::Constant(n, 5));
    // These can be reordered for spatial cache locality
    // But I will edit this later since this function
    // is only called once.

    // Set up Diriclet BC:
    for (int i = 0; i < M-1; ++i) {
      int top = sub2ind(i, 0);
      int bottom = sub2ind(i, N-1);
      mat.insert(top, top) = 1.0;
      mat.insert(bottom, bottom) = 1.0;
    }
    // interior diagonal:
    for (int i = 0; i < M-1; ++i) {
      for (int j = 1; j < N-1; ++j) {
        int index = sub2ind(i, j);
        mat.insert(index, index) = diag;
      }
    }
    // interior off-diagonal:
    for (int i = 0; i < M-1; ++i) {
      for (int j = 1; j < N-1; ++j) {
        int index = sub2ind(i, j);
        int im = sub2ind((i==0) ? (M-2) : (i-1), j);
        int ip = sub2ind((i==M-2) ? 0 : (i+1), j);
        int jm = sub2ind(i, j-1);
        int jp = sub2ind(i, j+1);
        mat.insert(index, im) = offdiag;
        mat.insert(index, ip) = offdiag;
        mat.insert(index, jm) = offdiag;
        mat.insert(index, jp) = offdiag;
      }
    }
    mat.makeCompressed();
  }

  void FBESolver::FillRHS(double t)
  {
    // Dirichlet BC:
    for (int i = 0; i < M-1; ++i) {
      int top = sub2ind(i, 0);
      int bottom = sub2ind(i, N-1);
      rhs[top] = 1.0;
      rhs[bottom] = 0.0;
    }
    // interior:
    for (int i = 0; i < M-1; ++i) {
      for (int j = 1; j < N-1; ++j) {
        int index = sub2ind(i, j);
        double S = source(param, x[i], y[j], t);
        rhs[index] = sol[index] + 
          dt*(-4*pow(sol[index], 3) + 6*sol[index]*sol[index] 
              - 2*sol[index] + S);
      }
    }
  }

  void FBESolver::Integrate(double t_final)
  {
    // Initial values:
    for (int i = 0; i < M-1; ++i) {
      for (int j = 0; j < N; ++j) {
        sol[sub2ind(i, j)] = eta(param, x[i], y[j], 0.0);
      }
    }
    // Setting up solver:
    ConstructMatrix();
    Eigen::SparseLU<SpMat> solver;
    // solver.analyzePattern(mat);
    // solver.factorize(mat);
    // if (solver.info() != Eigen::Success)
    //   std::cerr << "factorize failed!!!" << std::endl;
    solver.compute(mat);

    double t = 0.0;
    double t_check = t_final - dt/2;

    #ifdef PRINT_ITER
    printf("t = %.6f", t);
    #endif
    while (t < t_check) {
      FillRHS(t);
      sol = solver.solve(rhs);
      t += dt;
      #ifdef PRINT_ITER
      printf("\rt = %.6f", t);
      #endif
    }
  }

  Matrix<double> FBESolver::GetAnalytic(double t = default_tf)
  {
    Matrix<double> sol(M, N);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < N; ++j) {
        sol(i, j) = eta(param, x[i], y[j], t);
      }
    }
    return sol;
  }

}

#endif