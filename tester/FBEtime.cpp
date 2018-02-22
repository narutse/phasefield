#include "FBEsolver.hpp"
#include <fstream>
#include <utility>

double L2 (const Matrix<double> &U, const Matrix<double> &V, 
  double dx, double dy) 
{  
  int MN = U.size();
  double sum = 0.0;
  for (int i = 0; i < MN; ++i) {
    sum += std::pow(U[i] - V[i], 2);
  }
  return std::sqrt(sum*dx*dy);
}

int main() {
  const MMSfunc::Param param;
  const int M = 1025, N = 513;
  const double DT[] = {
    ldexp(1, -2), ldexp(1, -3), 
    ldexp(1, -4), ldexp(1, -5),
    ldexp(1, -6), ldexp(1, -7), 
    ldexp(1, -8), ldexp(1, -9), 
    ldexp(1, -10), ldexp(1, -11),
    ldexp(1, -12), ldexp(1, -13),    
  };
  
  FILE *out = fopen("data/FBEtimeL2.csv", "w");
  for (auto &dt : DT) {
    printf("Running with M = %d, N = %d, dt = %e\n", M, N, dt);
    auto solver = Solver::FBESolver(param, M, N, dt);
    double dx = solver.GetDx(), dy = solver.GetDy();
    auto analytic = solver.GetAnalytic();
    auto numerical = solver.Solve();
    double l2 = L2(analytic, numerical, dx, dy);
    printf("\rdt = %e, L2 = %e\n", dt, l2);
    fprintf(out, "%e, %e\n", dt, l2);
    fflush(out);
  }
  fclose(out);
}
