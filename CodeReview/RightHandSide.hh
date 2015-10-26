#ifndef _RIGHTHANDSIDE_HH_
#define _RIGHTHANDSIDE_HH_

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim>(dim+1) {}
    virtual double value (const Point<dim> &p, const unsigned int component = 0) const;
    virtual void vector_value (const Point<dim> &p, Vector<double> &value) const;
};

#endif
