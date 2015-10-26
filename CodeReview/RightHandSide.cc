#include "RightHandSide.hh"


template <int dim>
double RightHandSide <dim>::value (const Point<dim> & p, const unsigned int  component) const
{
  return 0;
}


template <int dim>
void
RightHandSide <dim>::vector_value (const Point<dim> &p, Vector<double> &values) const
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = RightHandSide::value (p, c);
}

//Explicit initialization

template class RightHandSide<2>;
