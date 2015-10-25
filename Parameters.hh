#ifndef _PARAMETERS_HH_
#define _PARAMETERS_HH_

#include <deal.II/base/parameter_handler.h>

using namespace dealii;

// Example 1 is the first example of "hp-Adaptive Discontinuous Galerkin Finite
// Element Methods for the Stokes Problem."
// Example 2 is the second example of "An adaptive Uzawa FEM for the Stokes
// Problem: Convergence without the Inf-Sup Condition"
// Example 3 is the first example of "An adaptive Uzawa FEM for the Stokes
// Problem: Convergence without the Inf-Sup Condition"
// Example 4 is the third example of "An adaptive Uzawa FEM for the Stokes
// Problem: Convergence without the Inf-Sup Condition"
enum EXAMPLE {example_1,example_2,example_3,example_4};
enum QUADRATURE {gauss_legendre,gauss_lobatto};
enum REFINEMENT {h_refine,p_refine,hp_refine};

class Parameters
{
public :
  Parameters(const std::string &input_filename);

  bool get_verbose() const;

  EXAMPLE get_example() const;

  QUADRATURE get_quadrature() const;

  REFINEMENT get_refinement() const;

  unsigned int get_dim() const;

  unsigned int get_max_degree() const;

  unsigned int get_max_n_cycles() const;

  double get_theta() const;

  double get_tolerance() const;

private :
  void declare_parameters(ParameterHandler &prm);

  void parse_parameters(ParameterHandler &prm);

  ParameterHandler prm;
  bool verbose;
  EXAMPLE example;
  QUADRATURE quadrature;
  REFINEMENT refinement;
  unsigned int dim;
  unsigned int max_degree;
  unsigned int max_n_cycles;
  double theta;
  double tolerance;
};


inline
bool Parameters::get_verbose() const
{
  return verbose;
}


inline
EXAMPLE Parameters::get_example() const
{
  return example;
}


inline
QUADRATURE Parameters::get_quadrature() const
{
  return quadrature;
}


inline
REFINEMENT Parameters::get_refinement() const
{
  return refinement;
}


inline
unsigned int Parameters::get_dim() const
{
  return dim;
}


inline
unsigned int Parameters::get_max_degree() const
{
  return max_degree;
}


inline
unsigned int Parameters::get_max_n_cycles() const
{
  return max_n_cycles;
}


inline
double Parameters::get_theta() const
{
  return theta;
}


inline
double Parameters::get_tolerance() const
{
  return tolerance;
}

#endif
