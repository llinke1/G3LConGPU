#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <complex>

#include "ApertureStatistics.h"
#include "helpers.h"

/**
 * Program computing <N^2 M_ap> for predefined thetas
 *
 * Output is in units of Gtilde (unitless or Msun/Mpc^2)
 * Usage: calculateApertureStatistics.x filename_Gtilde filename_Thetas Tesselated?
 *
 * Example: calculateApertureStatistics.x 
 * ../products/gtilde_allTiles/all.gtilde.dat ../products/NNMap/thetas.dat 0
 * 
 * @author Laila Linke, llinke@astro.uni-bonn.de
 */


int main(int argc, char* argv[])
{
 
  int n_params = 3; //Expected number of params
  std::string usage="calculateNNM.x filename_Gtilde filename_Thetas Tesselated?"; //Usage description

  std::string example="calculateNNM.x gtilde_allTiles/all.gtilde.dat thetas.dat 0"; //Example usage

  // Check Number of CMD Line arguments
  g3lcong::checkCmdLine(argc, n_params, usage, example);

  // Reading in from commandline
  std::string filename_Gtilde = argv[1]; //Name of File containing Gtilde
  std::string filename_Thetas = argv[2]; //Name of File with thetas for which NNMap shall be calculated
  bool tesselation = std::stoi(argv[3]); //Bool if Gtilde was tesselated or not (0 if not tesselated, 1 else)

  //Initializing Aperture Statistics and read in Gtilde
  g3lcong::ApertureStatistics aperture_statistics(filename_Gtilde, tesselation);
  
  // Reading in Thetas for which N2Map should be calculated and do calculation
  std::ifstream input(filename_Thetas);
  double theta1,theta2,theta3;
  while(input>>theta1>>theta2>>theta3)
    {
      // Calculate Aperture Statistics
      std::complex<double> NNMap = aperture_statistics.NNM(theta1, theta2, theta3);

      // Print out Aperture Statistics
      std::cout<<theta1<<" "<<theta2<<" "<<theta3<<" "
	       <<std::real(NNMap)<<" "<<std::imag(NNMap)<<std::endl;
    };
  

  return 0;

  
}
