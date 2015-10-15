#include "ExactSolutionEx1.hh"
#include <cmath>

#include <deal.II/lac/vector.h> 

template <int dim>
void ExactSolutionEx1 <dim>::vector_value (const Point<dim> &p,
    Vector<double>   &values) const
{

  values(0) = -1.0 * std::exp (p(0)) * (p(1) * std::cos (p(1)) + std::sin (p(1)));
  values(1) =    std::exp(p(0)) * p(1) * std::sin (p(1));
  values(2) = 2.0 * std::exp (p(0)) * std::sin (p(1)) - (2.0 * (1.0 - std::exp(1.0)) * (std::cos (1.0) - 1.0)) / 3.0;

}


template <int dim>
void ExactSolutionEx1 <dim>::vector_gradient  (const Point<dim> &p,
    std::vector< Tensor< 1, dim > > &gradients) const
{

  gradients[0][0]= -1.0 * std::exp (p(0)) * (p(1) * std::cos (p(1)) + std::sin (p(1)));
  gradients[0][1]= -1.0 * std::exp (p(0)) * ( 2.0 * std::cos (p(1)) - p(1) * std::sin (p(1)));
  gradients[1][0]=  std::exp(p(0)) * p(1) * std::sin (p(1));
  gradients[1][1]=  std::exp(p(0)) * (std::sin (p(1)) + p(1) * std::cos(p(1)));
  gradients[2][0]=  2.0 * std::exp (p(0)) * std::sin (p(1));
  gradients[2][1]=  2.0 * std::exp (p(0)) * std::cos (p(1));

}

//Explicit initialization

template class ExactSolutionEx1<2>;
