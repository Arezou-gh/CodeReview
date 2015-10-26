#ifndef _RIGHTHANDSIDEEx1_HH_
#define _RIGHTHANDSIDEEx1_HH_

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

template <int dim>
class RightHandSideEx1 : public Function<dim>
{
public:
  RightHandSideEx1() : Function<dim>(dim+1) {}
  virtual double value (const Point<dim> &p, const unsigned int component) const;
  virtual void vector_value (const Point<dim> &p, Vector<double> &value) const;
};

#endif
