#include "RightHandSideEx1.hh"


template <int dim>
double RightHandSideEx1<dim>::value(const Point<dim> &p, const unsigned int component) const
{
  // Silence warnings
  (void) p;
  (void) component;

  return 0;
}


template <int dim>
void
RightHandSideEx1<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = this->value(p,c);
}


//Explicit initialization
template class RightHandSideEx1<2>;
template class RightHandSideEx1<3>;
