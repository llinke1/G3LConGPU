#include "TwoPtStat.h"
#include "Kdtree.h"
#include "Function.h"

#include <cmath>

#include <iostream>
#include <string>

/**
 * Program computing the paircount 
 * Needs as input stddev of redshiftweighting
 * If stddev of redshiftweighting is set to 0, no redshift weighting is done
 *
 * Usage: calculatePaircount_kdtrees.x catalog1 catalog2 weight-z-file theta_min theta_max number_bins
 *
 * Example: calculatePaircount_kdtrees.x ../../Data/KIDS_129.0_-0.5.objects.dat
 * ../../Data/KIDS_129.0_-0.5.objects.dat none
 *
 * @author Laila Linke, llinke@astro.uni-bonn.de
 */

int main(int argc, char* argv[])
{
 

  int n_params=6; // Expected number of params
  std::string usage="calculatePaircount_kdtrees.x filename_galaxies1 filename_galaxies2 sigmaZ theta_min theta_max num_bins"; //Usage description

  std::string example="calculatePaircount_kdtrees.x ../../Data/KIDS_129.0_-0.5.objects.dat ../../Data/KIDS_129.0_-0.5.objects.dat 0 0.15 0.79 100"; //Example usage
  
  // Check Number of CMD Line arguments
  g3lcong::checkCmdLine(argc, n_params, usage, example);
  
  
  
  // Reading in commandline
  std::string filename1 = argv[1]; // File with first galaxy catalog
  std::string filename2 = argv[2]; // File with second galaxy catalog
  double sigmaZ = std::stod(argv[3]); //Width of redshift weighting Gaussian, if 0: no weighting
  double theta_min = std::stod(argv[4]); //[arcmin]
  double theta_max = std::stod(argv[5]); //[arcmin]
  int num_bins = std::stoi(argv[6]); 

  // Setting Binsize
  double theta_binsize = (log(theta_max) - log(theta_min))/num_bins;
  g3lcong::TwoPtStat two_pt_stat(theta_min, theta_max, num_bins);

  // Building Trees from catalog
  g3lcong::Kdtree tree1(filename1, 6, 1, 2, 5, 3, 4, 6);
  g3lcong::Kdtree tree2(filename2, 6, 1, 2, 5, 3, 4, 6);

   // Calculating Paircount
  two_pt_stat.dualTreeCount(tree1.getRoot(), tree2.getRoot(), sigmaZ);


  // Calculate Normalization
  double normalization;
  if (filename1 == filename2)
    {
      normalization = tree1.getSize()*(tree2.getSize()-1);
    }
  else
    {
      normalization = tree1.getSize()*tree2.getSize();
    };

  // Print out Paircount
  for(unsigned int i=0; i<two_pt_stat.counts_total_.size(); i++)
    {
      std::cout<<theta_min*exp(i*theta_binsize)<<" "<<two_pt_stat.counts_total_.at(i)<<" "<<normalization<<std::endl;
    }

  return 0;
  
}
