#include "ApertureStatistics.h"

#include <complex>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

std::tuple<double, double> g3lcong::ApertureStatistics::a_Theta8(const double& theta1, const double& theta2, const double& theta3)
{
  double theta1sq=theta1*theta1;
  double theta2sq=theta2*theta2;
  double theta3sq=theta3*theta3;
  double Theta8=(theta1sq*theta2sq+theta1sq*theta3sq+theta2sq*theta3sq)/3.;
  double a=2./3.*theta1sq*theta2sq*theta3sq/Theta8;

  return std::make_tuple(a, Theta8);

}

std::complex<double> g3lcong::ApertureStatistics::A_NNM(const double& vartheta1, const double& vartheta2, const double& phi,
							const double& theta1, const double& theta2,
							const double& a2, const double& bigT)
{

  // See Schneider & Watts(2005) for these formulae
  
  std::complex <double> eiphi3(cos(phi),sin(phi));

  std::complex<double>  b0  = 
    vartheta1*vartheta1/2./theta1/theta1
    +vartheta2*vartheta2/2./theta2/theta2
    -a2/4./theta1/theta1/theta2/theta2
    *(vartheta1*vartheta1*theta2*theta2/theta1/theta1
      +vartheta2*vartheta2*theta1*theta1/theta2/theta2
      +2.*vartheta1*vartheta2*cos(phi));

  std::complex<double> g1   = vartheta1-.5*a2*(vartheta1/theta1/theta1+vartheta2/eiphi3/theta2/theta2);
  std::complex<double> g2   = vartheta2-.5*a2*(vartheta2/theta2/theta2+vartheta1*eiphi3/theta1/theta1);
  std::complex<double> F1   = 2.*theta1*theta1-g1*conj(g1);
  std::complex<double> F2   = 2.*theta2*theta2-g2*conj(g2);

  std::complex<double> ANNM = ((g1-vartheta1)
			       *(g2-vartheta2)
			       *(F1*F2/a2 -(F1+F2) +2.*a2+conj(g1)*g2/eiphi3
				 +g1*conj(g2)*eiphi3));
  
  ANNM-=(((g2-vartheta2)+(g1-vartheta1)*eiphi3)
	 *(g1*(F2-2.*a2)+g2*(F1-2.*a2)/eiphi3));
  ANNM+= 2.*g1*g2*a2;
  ANNM*= exp(-real(b0))/(72.*3.1416*bigT*bigT);

  return ANNM;
}

void g3lcong::ApertureStatistics::readGtilde(std::string filename, bool tesselated)
{
  // Reading in
  std::ifstream input(filename);

  if(!input)
    {
      std::cerr<<"ApertureStatistics::readGtilde: Couldn't open file "<<filename<<" .Exiting\n";
      exit(1);
    };
  
  if(tesselated)
    {
      double val1,val2,val3,val4,val5,val6, val7;
      while(input>>val1>>val2>>val3>>val4>>val5>>val6>>val7)
	{
	  vartheta1_.push_back(val1);
	  vartheta2_.push_back(val2);
	  phi_.push_back(val3);
	  V_.push_back(val4);
	  std::complex<double> G(val5,val6);
	  Gtilde_.push_back(G);
	};
    }
  else
    {
      double val1,val2,val3,val4,val5,val6, val7, val8, val9;
      while(input>>val1>>val2>>val3>>val4>>val5>>val6>>val7>>val8>>val9)
	{
	  vartheta1_.push_back(val1);
	  vartheta2_.push_back(val2);
	  phi_.push_back(val3);
	  V_.push_back(val4*val5*val6);
	  std::complex<double> G(val7,val8);
	  Gtilde_.push_back(G);
	};
    };
}


g3lcong::ApertureStatistics::ApertureStatistics(std::string filenameGtilde, bool tesselated)
{
  readGtilde(filenameGtilde, tesselated);
}



std::complex<double> g3lcong::ApertureStatistics::NNM(const double& theta1, const double& theta2, const double& theta3)
{
  //Precompute a and Theta8
  double a, Theta8;
  std::tie(a, Theta8) = a_Theta8(theta1, theta2, theta3);

  // Calculate N2Map
  std::complex<double> result(0.,0.);
  
  for(unsigned int i=0; i<Gtilde_.size(); i++)
    {
    
      result += A_NNM(vartheta1_.at(i), vartheta2_.at(i), phi_.at(i), theta1, theta2, a, Theta8) * Gtilde_.at(i)
   	*V_.at(i)*vartheta1_.at(i)*vartheta2_.at(i);
    }

  return result;
}
