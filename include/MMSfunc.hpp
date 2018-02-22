#ifndef __MMSfunc_HPP__
#define __MMSfunc_HPP__

#include <cmath>

namespace MMSfunc {
  // default values. They are used for this benchmark problem.
  struct Param {
    double kappa = 0.0004;
    double A1 = 0.0075;
    double B1 = 8.0*M_PI;
    double A2 = 0.03;
    double B2 = 22.0*M_PI;
    double C2 = 0.0625*M_PI;
  };

  const auto alpha = [](const Param& p, auto x, auto t) 
  {
    return p.A1*t*sin(p.B1*x) + p.A2*sin(p.B2*x + p.C2*t) + 0.25;
  };

  const auto eta = [](const Param &p, auto x, auto y, auto t) {
    return -0.5*tanh((y - alpha(p, x, t))/sqrt(2.0*p.kappa)) + 0.5;
  };

  const auto dadt = [](const Param &p, auto x, auto t)
  {
    return p.A1*sin(p.B1*x) + p.A2*p.C2*cos(p.B2*x + p.C2*t);
  };

  const auto dadx = [](const Param &p, auto x, auto t)
  {
    return p.A1*p.B1*t*cos(p.B1*x) + p.A2*p.B2*cos(p.B2*x + p.C2*t);
  };

  const auto d2adx2 = [](const Param &p, auto x, auto t)
  {
    return - p.A1*pow(p.B1, 2)*t*sin(p.B1*x) 
            - p.A2*pow(p.B2, 2)*sin(p.B2*x + p.C2*t);
  };

  const auto source = [](const Param &p, auto x, auto y, auto t)
  {  
    auto ta = tanh((y-alpha(p, x, t))/sqrt(2.0*p.kappa));
    auto factor = 0.25*(1.0 - ta*ta)/sqrt(p.kappa);
    auto first_term = -2.0*sqrt(p.kappa)*ta*pow(dadx(p, x, t), 2);
    auto second_term = sqrt(2.0)*(dadt(p, x, t) - p.kappa * d2adx2(p, x, t));
    return factor*(first_term + second_term);
  };
}

#endif