#include "ExactSolutionEx2.hh"
#include <cmath>

// Example-2 (Taken from numerical experiments of paper Houston et al.)
// In this example pressure and gradient of velocity both are singular on re-entrant corner (0,0)
// We consider the L-shape [-1,1]*[-1,1]\(0,1)*(-1,0) as problem domain

using std::cos;
using std::sin;
using std::pow;
using std::exp;
using std::sqrt;

template <>
void ExactSolutionEx2<2>::vector_value(const Point<2> &p, Vector<double> &values) const
{

  const double x(p[0]);
  const double y(p[1]);

  const double PI = std::atan(1.0)*4;

  const double Lambda= 0.544;
  const double Omega=(3*(PI))/2;
  double r= sqrt(pow(x,2)+pow(y,2));
  double Phi;

  // in order to get rid of NAN(Not A Number) case
  if (std::fabs(x-0.0)<1e-12 &&  std::fabs(y-0.0)<1e-12 )
    Phi=PI/4;
  else if (x-0.0>1e-12 && std::fabs(y-0.0)<1e-12)
    Phi=0.0;
  else if (x-0.0<1e-12 && std::fabs(y-0.0)<1e-12)
    Phi=PI;
  else if (std::fabs(x-0.0)<1e-12 && y-0.0>1e-12)
    Phi=PI/2;
  else if (std::fabs(x-0.0)<1e-12 && y-0.0<1e-12)
    Phi=3*PI/2;
  else if ((x-0.0<1e-12 && y-0.0>1e-12) || (x-0.0<1e-12 && y-0.0<1e-12))
    Phi=std::atan(p(1)/p(0))+PI;
  else
    Phi=std::atan(p(1)/p(0));


  double Psi= (sin((1+Lambda)*(Phi))* cos(Lambda*Omega))/(1+Lambda)-cos((1+Lambda)*Phi)
              -(sin((1-Lambda)*Phi) * cos(Lambda*Omega))/(1-Lambda) + cos((1-Lambda)*Phi);

  double Deriv_Psi= (cos(Lambda*Omega) * cos((1+Lambda)*Phi))+ (1+Lambda)*sin((1+Lambda)*Phi)-
                    cos(Lambda*Omega)*cos((1-Lambda)*Phi)-(1-Lambda)*sin((1-Lambda)*Phi);

  double third_Deriv_Psi= - (pow((1+Lambda),2)) * cos(Lambda*Omega)* cos((1+Lambda)*Phi)- pow((1+Lambda),3)* sin((1+Lambda)*Phi)+ pow((1-Lambda),2)* cos(Lambda*Omega)* cos((1-Lambda)*Phi)+ pow((1-Lambda),3)* sin((1-Lambda)* Phi);

  //L-shape/ actual problem / solution from Houston, et.al's paper

  values(0) = pow(r, Lambda) *((1+Lambda)*sin(Phi) * Psi + cos(Phi) *
                               Deriv_Psi ) ;

  values(1) =  pow(r, Lambda) * (sin(Phi)* Deriv_Psi - (1+Lambda) * cos(Phi) * Psi);

  if (std::fabs(x-0.0)<1e-12 &&  std::fabs(y-0.0)<1e-12 )
    // since in this case, presssure would be infinity...Hence, we have to take care of infinity value of pressure by setting it with a large number.
    values(2) = pow (10,100);
  else
    values(2) = -pow(r,(Lambda-1)) * (pow((1+Lambda),2)*Deriv_Psi + third_Deriv_Psi )/ (1-Lambda);
}


template <>
void ExactSolutionEx2<2>::vector_gradient(const Point<2> &p,
                                          std::vector<Tensor<1,2>> &gradients) const
{

  const double x(p[0]);
  const double y(p[1]);
  const double PI = std::atan(1.0)*4;

  const double Lambda= 0.544;
  const double Omega=(3*(PI))/2;
  double r= sqrt(pow(x,2)+pow(y,2));
  double Phi;

  // in order to get rid of NAN(Not A Number) case
  if (std::fabs(x-0.0)<=1e-12 &&  std::fabs(y-0.0)<=1e-12 )
    Phi=PI/4;
  else if (x-0.>1e-12 && std::fabs(y-0.0)<=1e-12)
    Phi=0.0;
  else if (x-0.<1e-12 && std::fabs(y-0.0)<=1e-12)
    Phi=PI;
  else if (std::fabs(x-0.0)<=1e-12 && y-0.0>1e-12)
    Phi=PI/2;
  else if (std::fabs(x-0.0)<=1e-12 && y-0.0<1e-12)
    Phi=3*PI/2;
  else if ((x-0.0<1e-12 && y-0.0>1e-12) || (x-0.0<1e-12 && y-0.0<1e-12))
    Phi=std::atan(p(1)/p(0))+PI;
  else
    Phi=std::atan(p(1)/p(0));


  //L-shape/ actual problem / solution from Houston, et.al's paper


  if (std::fabs(x-0.0)<1e-12 &&  std::fabs(y-0.0)<1e-12 )
    gradients[0][0]= pow (10,100);
  else
    gradients[0][0]= Lambda* pow(r, (Lambda-1))*( (1+Lambda)*sin(Phi) *( sin((1+Lambda)*Phi)*cos(Lambda*Omega)/(1+Lambda) -cos((1+Lambda)*Phi) - cos(Lambda*Omega) *sin((1-Lambda)*Phi)/(1-Lambda)+
                                                  cos((1-Lambda)*Phi) ) + cos(Phi)* ( cos(Lambda*Omega) *cos((1+Lambda)*Phi) +(1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega)*cos((1-Lambda)*Phi) - (1-Lambda)*sin((1-Lambda)*Phi)  ) )* x/r
                     +( pow(r,Lambda)*( (1+Lambda)*cos(Phi)* (  sin((1+Lambda)*Phi) *cos(Lambda*Omega)/(1+Lambda) - cos((1+Lambda)*Phi) - cos(Lambda*Omega) *sin((1-Lambda)*Phi) /(1-Lambda) + cos((1-Lambda)*Phi) ) + (1+Lambda)*sin(Phi)*( cos((1+Lambda)*Phi) *cos(Lambda*Omega) + (1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega)*cos((1-Lambda)*Phi) - (1-Lambda)*sin((1-Lambda)*Phi) )-
                                        sin(Phi) *( cos(Lambda*Omega)*cos((1+Lambda)*Phi) +(1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega)*cos((1-Lambda)*Phi) -(1-Lambda)*sin((1-Lambda)*Phi) ) + cos(Phi)*
                                        (-(1+Lambda)*cos(Lambda*Omega)*sin((1+Lambda)*Phi) + pow((1+Lambda),2)*cos((1+Lambda)*Phi) + (1-Lambda) *cos(Lambda*Omega)*sin((1-Lambda)*Phi) - pow((1-Lambda),2) *cos((1-Lambda)*Phi))))
                     * (-y)/pow(r,2);

  if (std::fabs(x-0.0)<1e-12 &&  std::fabs(y-0.0)<1e-12 )
    gradients[0][1]= pow (10,100);
  else
    gradients[0][1]= Lambda* pow(r, (Lambda-1))*( (1+Lambda)*sin(Phi) *( sin((1+Lambda)*Phi)*cos(Lambda*Omega)/(1+Lambda) -cos((1+Lambda)*Phi) - cos(Lambda*Omega) *sin((1-Lambda)*Phi)/(1-Lambda)+
                                                  cos((1-Lambda)*Phi) ) + cos(Phi)* ( cos(Lambda*Omega) *cos((1+Lambda)*Phi) +(1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega)*cos((1-Lambda)*Phi) - (1-Lambda)*sin((1-Lambda)*Phi)  ) )* y/r
                     + ( pow(r,Lambda)*( (1+Lambda)*cos(Phi)* (  sin((1+Lambda)*Phi) *cos(Lambda*Omega)/(1+Lambda) - cos((1+Lambda)*Phi) - cos(Lambda*Omega) *sin((1-Lambda)*Phi) /(1-Lambda) + cos((1-Lambda)*Phi) ) + (1+Lambda)*sin(Phi)*( cos((1+Lambda)*Phi) *cos(Lambda*Omega) + (1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega)*cos((1-Lambda)*Phi) - (1-Lambda)*sin((1-Lambda)*Phi) )-
                                         sin(Phi) *( cos(Lambda*Omega)*cos((1+Lambda)*Phi) +(1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega)*cos((1-Lambda)*Phi) -(1-Lambda)*sin((1-Lambda)*Phi) ) + cos(Phi)*
                                         (-(1+Lambda)*cos(Lambda*Omega)*sin((1+Lambda)*Phi) + pow((1+Lambda),2)*cos((1+Lambda)*Phi) + (1-Lambda) *cos(Lambda*Omega)*sin((1-Lambda)*Phi) - pow((1-Lambda),2) *cos((1-Lambda)*Phi))))
                     * x/pow(r,2);



  if (std::fabs(x-0.0)<1e-12 &&  std::fabs(y-0.0)<1e-12 )
    gradients[1][0]= pow (10,100);
  else
    gradients[1][0]= Lambda*pow(r, (Lambda-1)) *( sin(Phi)*( cos(Lambda*Omega) *cos((1+Lambda)*Phi) + (1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega)*cos((1-Lambda)*Phi) - (1-Lambda)*sin((1-Lambda)*Phi)  ) - (1+Lambda)*cos(Phi)*( sin((1+Lambda)*Phi) *cos(Lambda*Omega)/(1+Lambda) - cos((1+Lambda)*Phi) - cos(Lambda*Omega)*sin((1-Lambda)*Phi)/(1-Lambda) + cos((1-Lambda)*Phi))  ) * x/r
                     + pow(r,Lambda)* ( cos(Phi)* (cos(Lambda*Omega)*cos((1+Lambda)*Phi) +(1+Lambda) *sin((1+Lambda)*Phi) - cos(Lambda*Omega) *cos((1-Lambda)*Phi) - (1-Lambda) *sin((1-Lambda)*Phi))
                                        + sin(Phi) *( -cos(Lambda*Omega)*(1+Lambda)*sin((1+Lambda)*Phi) + pow((1+Lambda),2) *cos((1+Lambda)*Phi) + cos(Lambda*Omega) *(1-Lambda)* sin((1-Lambda)*Phi) - pow((1-Lambda),2) *cos((1-Lambda)*Phi))
                                        + (1+Lambda)*sin(Phi) *( sin((1+Lambda)*Phi) *cos(Lambda*Omega)/ (1+Lambda) - cos((1+Lambda)*Phi) - cos(Lambda*Omega)*sin((1-Lambda)*Phi)/(1-Lambda) + cos((1-Lambda)*Phi))
                                        -(1+Lambda)*cos(Phi) *( cos(Lambda*Omega)* cos((1+Lambda)*Phi) + (1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega) *cos((1-Lambda)*Phi) -(1-Lambda) *sin((1-Lambda)*Phi)) ) *(-y/pow(r,2));


  if (std::fabs(x-0.0)<1e-12 &&  std::fabs(y-0.0)<1e-12 )
    gradients[1][1]= pow (10,100);
  else
    gradients[1][1]=  Lambda*pow(r, (Lambda-1)) *( sin(Phi)*( cos(Lambda*Omega) *cos((1+Lambda)*Phi) + (1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega)*cos((1-Lambda)*Phi) - (1-Lambda)*sin((1-Lambda)*Phi)  ) - (1+Lambda)*cos(Phi)*( sin((1+Lambda)*Phi) *cos(Lambda*Omega)/(1+Lambda) -cos((1+Lambda)*Phi) - cos(Lambda*Omega)*sin((1-Lambda)*Phi)/(1-Lambda) + cos((1-Lambda)*Phi))  ) * y/r
                      + pow(r,Lambda)* ( cos(Phi)* (cos(Lambda*Omega)*cos((1+Lambda)*Phi) +(1+Lambda) *sin((1+Lambda)*Phi) - cos(Lambda*Omega) *cos((1-Lambda)*Phi) - (1-Lambda) *sin((1-Lambda)*Phi))
                                         + sin(Phi) *( -cos(Lambda*Omega)*(1+Lambda)*sin((1+Lambda)*Phi) + pow((1+Lambda),2) *cos((1+Lambda)*Phi) + cos(Lambda*Omega) *(1-Lambda)* sin((1-Lambda)*Phi) - pow((1-Lambda),2) *cos((1-Lambda)*Phi))
                                         + (1+Lambda)*sin(Phi) *( sin((1+Lambda)*Phi) *cos(Lambda*Omega)/ (1+Lambda) -cos((1+Lambda)*Phi) - cos(Lambda*Omega)*sin((1-Lambda)*Phi)/(1-Lambda) + cos((1-Lambda)*Phi))
                                         -(1+Lambda)*cos(Phi) *( cos(Lambda*Omega)* cos((1+Lambda)*Phi) + (1+Lambda)*sin((1+Lambda)*Phi) - cos(Lambda*Omega) *cos((1-Lambda)*Phi) -(1-Lambda) *sin((1-Lambda)*Phi)) ) *(x/pow(r,2));


  gradients[2][0]= 0 ;
  gradients[2][1]=  0 ;
}



template <>
void ExactSolutionEx2<3>::vector_value(const Point<3> &p, Vector<double> &values) const
{
  // Silence warnings
  (void) p;
  (void) values;
}


template <>
void ExactSolutionEx2<3>::vector_gradient(const Point<3> &p,
                                          std::vector<Tensor<1,3>> &gradients) const
{
  // Silence warnings
  (void) p;
  (void) gradients;
}
