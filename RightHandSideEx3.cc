#include "RightHandSideEx3.hh"

#include <cmath>

using std::exp;
using std::cos;
using std::sin;
using std::pow;


template <int dim>
double RightHandSideEx3<dim>::value(const Point<dim> &p, const unsigned int component) const
{
  const double x(p[0]);
  const double y(p[1]);
  double val(0.);

  if (component==0)

    val = 8.*y*pow(x,2)*cos(pow(x,2)+pow(y,2)) + 16.*y*sin(pow(x,2)+pow(y,2)) +
          8.*pow(y,3)*cos(pow(x,2)+pow(y,2)) - 20.*exp(-10.*(pow(x,2)+pow(y,2)))*x;
  if (component==1)
    val = -8.*pow(x,3)*cos(pow(x,2)+pow(y,2)) - 16.*x*sin(pow(x,2)+pow(y,2)) -
          8.*pow(y,2)*x*cos(pow(x,2)+pow(y,2)) - 20.*exp(-10.*(pow(x,2)+pow(y,2)))*y;

  return val;
}


template <int dim>
void
RightHandSideEx3<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = this->value(p,c);
}


//Explicit initialization
template class RightHandSideEx3<2>;
template class RightHandSideEx3<3>;
