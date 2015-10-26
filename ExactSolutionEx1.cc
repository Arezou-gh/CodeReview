#include "ExactSolutionEx1.hh"
#include <cmath>
// Example-1, (Taken from numerical experiments of paper Houston et al.)
// The solution is smooth on L-shape domain in [-1,1]*[-1,1]\(0,1)*(-1,0)
template <>
void ExactSolutionEx1<2>::vector_value(const Point<2> &p, Vector<double> &values) const
{
  values(0) = -1.0 * std::exp(p(0)) * (p(1) * std::cos (p(1)) + std::sin (p(1)));
  values(1) = std::exp(p(0)) * p(1) * std::sin (p(1));
  values(2) = 2.0 * std::exp(p(0)) * std::sin(p(1)) - (2.0 * (1.0 - std::exp(1.0)) *
                                                       (std::cos(1.0) - 1.0)) / 3.0;
}


template <>
void ExactSolutionEx1<2>::vector_gradient(const Point<2> &p,
                                          std::vector<Tensor<1,2>> &gradients) const
{
  gradients[0][0]= -1.0 * std::exp(p(0)) * (p(1) * std::cos(p(1)) + std::sin(p(1)));
  gradients[0][1]= -1.0 * std::exp(p(0)) * ( 2.0 * std::cos(p(1)) - p(1) * std::sin(p(1)));
  gradients[1][0]=  std::exp(p(0)) * p(1) * std::sin(p(1));
  gradients[1][1]=  std::exp(p(0)) * (std::sin(p(1)) + p(1) * std::cos(p(1)));
  gradients[2][0]=  2.0 * std::exp(p(0)) * std::sin(p(1));
  gradients[2][1]=  2.0 * std::exp(p(0)) * std::cos(p(1));

}



template <>
void ExactSolutionEx1<3>::vector_value(const Point<3> &p, Vector<double> &values) const
{
  // Silence warnings
  (void) p;
  (void) values;
}


template <>
void ExactSolutionEx1<3>::vector_gradient(const Point<3> &p,
                                          std::vector<Tensor<1,3>> &gradients) const
{
  // Silence warnings
  (void) p;
  (void) gradients;
}
