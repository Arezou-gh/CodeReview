#include "ExactSolutionEx4.hh"
#include <cmath>

using std::exp;
using std::pow;
using std::sin;
using std::cos;

template <>
void ExactSolutionEx4<2>::vector_value(const Point<2> &p, Vector<double> &values) const
{
  (void) p;
  (void) values;
}


template <>
void ExactSolutionEx4<2>::vector_gradient(const Point<2> &p,
                                          std::vector<Tensor<1,2>> &gradients) const
{
  (void) p;
  (void) gradients;
}



template <>
void ExactSolutionEx4<3>::vector_value(const Point<3> &p, Vector<double> &values) const
{
  const double mu(1.);
  const double lambda(10.);
  const double x(p[0]);
  const double y(p[1]);
  const double z(p[2]);

  values(0) = 2.*y*cos(pow(x,2)+pow(y,2));
  values(1) = -2.*x*cos(pow(x,2)+pow(y,2));
  values(2) = 0.;
  values(3) = mu*exp(-lambda*(pow(x,2)+pow(y,2)+pow(z,2)));
}


template <>
void ExactSolutionEx4<3>::vector_gradient(const Point<3> &p,
                                          std::vector<Tensor<1,3>> &gradients) const
{
  const double mu(1.);
  const double lambda(10.);
  const double x(p[0]);
  const double y(p[1]);
  const double z(p[2]);

  gradients[0][0] = -4.*x*y*sin(pow(x,2)+pow(y,2));
  gradients[0][1] = 2.*cos(pow(x,2)+pow(y,2)) - 4.*pow(y,2)*sin(pow(x,2)+pow(y,2));
  gradients[0][2] = 0.;
  gradients[1][0] = -2.*cos(pow(x,2)+pow(y,2)) + 4.*pow(x,2)*sin(pow(x,2)+pow(y,2));
  gradients[1][1] = 4.*x*y*sin(pow(x,2)+pow(y,2));
  gradients[1][2] = 0.;
  gradients[2][0] = 0.;
  gradients[2][1] = 0.;
  gradients[2][2] = 0.;
  gradients[3][0] = -2.*mu*lambda*x*exp(-lambda*(pow(x,2)+pow(y,2)+pow(z,2)));
  gradients[3][1] = -2.*mu*lambda*y*exp(-lambda*(pow(x,2)+pow(y,2)+pow(z,2)));
  gradients[3][2] = -2.*mu*lambda*z*exp(-lambda*(pow(x,2)+pow(y,2)+pow(z,2)));
}
