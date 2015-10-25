#include "RightHandSideEx4.hh"

#include <cmath>

using std::exp;
using std::cos;
using std::sin;
using std::pow;


template <int dim>
double RightHandSideEx4<dim>::value(const Point<dim> &p, const unsigned int component) const
{

  const double mu(1.);
  const double lambda(10.);
  const double x(p[0]);
  const double y(p[1]);
  const double z(p[2]);
  double val(0.);

  if (component==0)

    val= 16*y*sin(pow(x,2)+pow(y,2))+8*pow(x,2)*y*cos(pow(x,2)+pow(y,2))+8*pow(y,3)*cos(pow(x,2)+pow(y,2))-2*mu*lambda*x*exp(-lambda*(pow(x,2)+pow(y,2)+pow(z,2)));
  if (component==1)

    val= -16*x*sin(pow(x,2)+pow(y,2)) -8*x*pow(y,2)*cos(pow(x,2)+pow(y,2))-8*pow(x,3)*cos(pow(x,2)+pow(y,2))-2*mu*lambda*y*exp(-lambda*(pow(x,2)+pow(y,2)+pow(z,2))) ;
  return val;

  if (component==2)

    val= -2*mu*lambda*z*exp(-lambda*(pow(x,2)+pow(y,2)+pow(z,2))) ;
  return val;
}


template <int dim>
void
RightHandSideEx4<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = this->value(p,c);
}


//Explicit initialization
template class RightHandSideEx4<2>;
template class RightHandSideEx4<3>;
