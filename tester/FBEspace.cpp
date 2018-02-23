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
  const double dt = ldexp(1.0, -13);
  const std::pair<int, int> MN[] = {
    {101, 51}, {129, 65},
    {201, 101}, {257, 129},
    {401, 201}, {513, 257},
    {801, 401}, {1025, 513}
  };
  
  FILE *out = fopen("data/FBEspaceL2.csv", "w");
  fprintf(out, "grid, L2\n");
  for (auto &p : MN) {
    printf("Running with M = %d, N = %d\n", p.first, p.second);
    auto solver = Solver::FBESolver(param, p.first, p.second, dt);
    double dx = solver.GetDx(), dy = solver.GetDy();
    auto analytic = solver.GetAnalytic();
    auto numerical = solver.Solve();
    double l2 = L2(analytic, numerical, dx, dy);
    printf("\rdx = %f, dy = %f, L2 = %e\n", dx, dy, l2);
    fprintf(out, "%e, %e\n", std::sqrt(dx*dy), l2);
    fflush(out);
  }
  fclose(out);
}
