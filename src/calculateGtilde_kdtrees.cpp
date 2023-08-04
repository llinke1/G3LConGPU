#include "Kdtree.h"
#include "ThreePtStat.h"
#include "constants.h"
#include "Function.h"
#include <string>
#include <iostream>

/**
 * Program computing Gtilde for one tile, using kdTrees
 *
 * Needs as input file with omega and stddev of redshiftweighting
 * If stddev of redshiftweighting is set to 0, no redshift weighting is done
 *
 *
 * Output has logarithmic binning in theta1 and theta2 and linear in phi
 *
 * Usage: calculateGtilde_kdtrees.x
 * filename_sources filename_lenses1 filename_lenses2
 * filename_omega theta_min theta_max num_bins
 * filename_weight_redshift filename_sigma_crit Physical?
 *
 *
 * Example: calculateGtilde_kdtrees.x  ../../Data/KIDS_129.0_-0.5.sources.dat
 * ../../Data/KIDS_129.0_-0.5.objects.dat ../../Data/KIDS_129.0_-0.5.objects.dat
 * ../products/omega_allTiles/all.omega.dat  0.15 79.9 100 "none" "none" 0
 *
 * @author Laila Linke, llinke@astro.uni-bonn.de
 */

int main(int argc, char *argv[])
{
  // Checking Command Line

  int n_params = 15;                                                                                                                                                                               // Expected number of params
  std::string usage = "calculateGtilde_kdtrees.x filename_sources filename_lenses1 filename_lenses2 filename_omega theta_min theta_max r_min r_max num_bins sigmaZ filename_sigma_crit Physical?"; // Usage description

  std::string example = "calculateGtilde_kdtrees.x  ../../Data/KIDS_129.0_-0.5.sources.dat ../../Data/KIDS_129.0_-0.5.objects.dat ../../Data/KIDS_129.0_-0.5.objects.dat ../products/omega_allTiles/all.omega.dat  0.15 79.9 0.1 300 100 0 none 0"; // Example usage

  // Check Number of CMD Line arguments
  g3lcong::checkCmdLine(argc, n_params, usage, example);

  // Read in command line arguments
  std::string filename_sources = argv[1]; // File with Source Galaxies
  std::string filename_lenses1 = argv[2]; // File with Lens1 Galaxies
  std::string filename_lenses2 = argv[3]; // File with Lens2 Galaxies

  std::string filename_omega = argv[4]; // File with Angular Correlation Func

  double theta_min = std::stod(argv[5]); // Min Theta for calc of Gtilde [arcmin]
  double theta_max = std::stod(argv[6]); // Max Theta for calc of Gtilde [arcmin]
  double r_min = std::stod(argv[7]);     // Min R for calc of Gtilde [Mpc] (ignored if not physical)
  double r_max = std::stod(argv[8]);     // Max R for calc of Gtilde [Mpc] (ignored if not physical)
  int num_bins = std::stoi(argv[9]);     // Number of Bins for Gtilde on ONE axis

  // Width of Redshift Weighting Gaussian, if 0:  no weighting
  double sigmaZ = std::stod(argv[10]);
  // File with source-averaged Sigma Crit, if "none": SigCrit is 1
  std::string filename_sigma_crit = argv[11];
  // File with angular diameter distance, if "none": distance is 1
  std::string filename_angular_distance = argv[12];

  bool physical = std::stoi(argv[13]); // Is 1 if physical Gtilde, 0 if angular Gtilde
  bool flipE1 = false;
  flipE1 = std::stoi(argv[14]); // If 1: Sign of ellipticity component 1 is flipped
  bool flipE2 = false;
  flipE2 = std::stoi(argv[15]);

  double phi_min = 0.0;             // Min Phi for calc of Gtilde [radians]
  double phi_max = 2 * g3lcong::pi; // Max Phi for calc of Gtilde [radians]

  // Building trees for sources and lenses
  g3lcong::Kdtree sources(filename_sources, 6, 1, 2, 5, 3, 4, 6, flipE1, flipE2);
  g3lcong::Kdtree lenses1(filename_lenses1, 6, 1, 2, 5, 3, 4, 6);
  g3lcong::Kdtree lenses2(filename_lenses2, 6, 1, 2, 5, 3, 4, 6);

  // Reading in Functions
  g3lcong::Function omega(filename_omega, 0.);
  g3lcong::Function sigma_crit(filename_sigma_crit, 1.);
  g3lcong::Function angular_distance(filename_angular_distance, 1.);

  // Initialize ThreePtStat
  g3lcong::ThreePtStat three_pt_stat(theta_min, theta_max, r_min, r_max, phi_min, phi_max, num_bins);

  // Calculate Gtilde
  if (physical)
  {
    three_pt_stat.tripleTreeCount(sources, lenses1, lenses2, &omega,
                                  sigmaZ);
  }
  else
  {
    three_pt_stat.tripleTreeCount(sources, lenses1, lenses2, &omega,
                                  sigmaZ, &sigma_crit, &angular_distance);
  }

  // Print out vartheta1, vartheta2, psi, binsize, Gtilde
  double phi_binsize = three_pt_stat.phi_binsize_;
  double theta_binsize = three_pt_stat.theta_binsize_;

  for (int i = 0; i < num_bins; i++)
  {
    // Theta1 Center of this bin
    double theta1 = 0.5 * (exp(log(theta_min) + theta_binsize * (i + 1)) + exp(log(theta_min) + theta_binsize * i));
    // Theta1 Binsize of this bin
    double deltaTheta1 = exp(log(theta_min) + theta_binsize * (i + 1)) - exp(log(theta_min) + theta_binsize * i);

    for (int j = 0; j < num_bins; j++)
    {
      // Theta2 Center of this bin
      double theta2 = 0.5 * (exp(log(theta_min) + theta_binsize * (j + 1)) + exp(log(theta_min) + theta_binsize * j));
      // Theta2 Binsize of this bin
      double deltaTheta2 = exp(log(theta_min) + theta_binsize * (j + 1)) - exp(log(theta_min) + theta_binsize * j);

      for (int k = 0; k < num_bins; k++)
      {
        // Phi Center of this bin
        double phi = (k + 0.5) * phi_binsize + phi_min;

        // Output
        std::cout
            << theta1 << " "                                                                     // bin center theta 1 [arcmin or Mpc]
            << theta2 << " "                                                                     // bin center theta 2 [arcmin or Mpc]
            << phi << " "                                                                        // phi center [radians]
            << deltaTheta1 << " "                                                                // bin size theta 1[arcmin or Mpc]
            << deltaTheta2 << " "                                                                // bin size theta 2[arcmin or Mpc]
            << phi_binsize << " "                                                                // phi bin size [radians]
            << three_pt_stat.G_real_total_.at(i * num_bins * num_bins + j * num_bins + k) << " " // Real part of Gtilde [dimensionless or Msun/Mpc²]
            << three_pt_stat.G_imag_total_.at(i * num_bins * num_bins + j * num_bins + k) << " " // Imaginary part of Gtilde [dimensionless or Msun/Mpc²]
            << three_pt_stat.weight_total_.at(i * num_bins * num_bins + j * num_bins + k)        // Weight of Gtilde [dimensionless or Msun/Mpc²]
            << std::endl;
      };
    };
  };
}
