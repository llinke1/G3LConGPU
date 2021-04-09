#ifndef G3LCONG_CONSTANTS_H
#define G3LCONG_CONSTANTS_H

#include <complex>

//For GPU Parallelisation, match this to maximum of computing GPU
#define THREADS 256 //Maximum Threads per Block
#define BLOCKS 184 //Maximum blocks for all SMs in GPU

namespace g3lcong
{
 /*****************************************/
  /* USEFUL NUMBERS ************************/
  /*****************************************/

  ///Pi
  const double pi=3.14159;

  ///Imaginary unit
  const std::complex<double> i(0.,1.);
  const std::complex<float> i_fl(0.,1.);

  
  ///1 Solarmass in kg
  const double solarmass_in_kg=1.98845e30;

  ///1 Mpc in m
  const double Mpc_in_m=3.086e22;

  ///1 rad in arcmin
  const double rad_in_arcmin=pi/180.0/60.0;
  
  ///Speed of light in m/s
  const double c=2.99792e8;

  ///Newtonion Gravitational Constant in m^3/s^2/kg
  const double G=6.673e-11; 


}

#endif //GGGL_CONSTANTS_H
