#include "ExactSolutionEx3.hh"
#include <cmath>

#include <deal.II/lac/vector.h> 
#include <stdio.h> 

template <int dim>
void ExactSolutionEx3 <dim>::vector_value (const Point<dim> &p,
		Vector<double>   &values) const
		{
	
	values(0) = 2*p(1)*std::cos(std::pow(p(0),2)+std::pow(p(1),2));
	values(1) =  -2*p(0)*std::cos(std::pow(p(0),2)+std::pow(p(1),2));
	values(2) = std::exp(-10*(std::pow(p(0),2)+std::pow(p(1),2)))- 0.3142;	
// 0.3142 : is the pressure mean value
		}

template <int dim>
void ExactSolutionEx3 <dim>::vector_gradient  (const Point<dim> &p,
		std::vector< Tensor< 1, dim > > &gradients) const
		{	
	//L-shape/ actual problem / solution from Houston, et.al's paper
	gradients[0][0]= -4*p(0)*p(1)*std::sin(std::pow(p(0),2)+std::pow(p(1),2));		
	gradients[0][1]= 2*std::cos(std::pow(p(0),2)+std::pow(p(1),2))-4*std::pow(p(1),2)*std::sin(std::pow(p(0),2)+std::pow(p(1),2));
	gradients[1][0]=  -2*std::cos(std::pow(p(0),2)+std::pow(p(1),2) ) + 4*std::pow(p(0),2)*std::sin(std::pow(p(0),2)+std::pow(p(1),2));
	gradients[1][1]=  4*p(0)*p(1)*std::sin(std::pow(p(0),2)+std::pow(p(1),2)) ;
	gradients[2][0]= std::exp(-10*(std::pow(p(0),2)+std::pow(p(1),2)) ) * -20*p(0)  ;
	gradients[2][1]= std::exp(-10*(std::pow(p(0),2)+std::pow(p(1),2)) ) * -20*p(1)  ;
		}

//Explicit initialization

template class ExactSolutionEx3 <2>;
