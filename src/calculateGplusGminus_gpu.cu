#include "helpers.h"
#include "constants.h"
#include "Function.h"
#include "kernelfunctions.cuh"
#include "cudaHelpers.cuh"
#include <string>
#include <iostream>

/**
 * Program computing Gtilde for one tile in angular units, using GPU parallelization
 * Needs as input file with omega and stddev of redshiftweighting
 * If stddev of redshiftweighting is set to 0, no redshift weighting is done
 *
 * Computation is done brute-force on submatrices of lens1-lens2-source space
 *
 * Output has logarithmic binning in theta1 and theta2 and linear in phi
 *
 *
 * @warning Number of parallel threads, blocks and size of submatrices is hard
 * coded in Precompiler flags and needs to be fit to GPU
 *
 * Usage: calculateGtilde_gpu.x filename_lenses1 filename_lenses2
 * filename_objects filename_omega theta_min [arcmin], theta_max [arcmin]
 * sigma
 *
 * Example: calculateGtilde_gpu.x KIDS_129.0_-0.5.objects.dat
 * KIDS_129.0_-0.5.objects.dat KIDS_129.0_-0.5.sources.dat
 * ../products/omega_allTiles/all.omega.dat 0.15 79.9 0
 *
 * @author Laila Linke, llinke@astro.uni-bonn.de
 */

int main(int argc, char *argv[])
{
  // Checking Command Line

  int n_params = 8;                                                                                                                                                                           // Expected number of params
  std::string usage = "calculateGplusGminus_gpu.x filename_sources filename_lenses1 filename_lenses2 theta_min theta_max num_bins flipE1? flipE2?"; // Usage description

  std::string example = "calculateGplusGminus_gpu.x  ../../Data/KIDS_129.0_-0.5.sources.dat ../../Data/KIDS_129.0_-0.5.objects.dat ../../Data/KIDS_129.0_-0.5.objects.dat 0.15 79.9 300 0 0"; // Example usage

  // Check Number of CMD Line arguments
  g3lcong::checkCmdLine(argc, n_params, usage, example);

  // Read in command line arguments
  std::string filename_lenses = argv[1];   // File with Source Galaxies
  std::string filename_sources1 = argv[2]; // File with Lens1 Galaxies
  std::string filename_sources2 = argv[3]; // File with Lens2 Galaxies

  double theta_min = std::stod(argv[4]); // Min Theta for calc of Gtilde [arcmin]
  double theta_max = std::stod(argv[5]); // Max Theta for calc of Gtilde [arcmin]
  int num_bins = std::stoi(argv[6]);     // Number of Bins for Gtilde on ONE axis

  bool flipE1 = false;
  flipE1 = std::stoi(argv[7]); // If 1: Sign of ellipticity component 1 is flipped
  bool flipE2 = false;
  flipE2 = std::stoi(argv[8]);

  double phi_min = 0.0;             // Min Phi for calc of Gtilde [radians]
  double phi_max = 2 * g3lcong::pi; // Max Phi for calc of Gtilde [radians]

  // Reading in galaxies and copying to device

  // x,y, z vectors
  std::vector<double> x1, y1, x2, y2, xL, yL, zL, e11, e21, e12, e22, w1, w2;

  if (g3lcong::readSources2Dev(filename_sources1, 6, 1, 2, 3, 4, 6, x1,
                               y1, e11, e21, w1, flipE1, flipE2))
    return 1;

  if (g3lcong::readSources2Dev(filename_sources2, 6, 1, 2, 3, 4, 6, x2,
                               y2, e12, e22, w2, flipE1, flipE2))
    return 1;

  if (g3lcong::readLenses2Dev(filename_lenses, 6, 1, 2, 5, xL, yL, zL))
    return 1;

  // Declare arrays for coordinates of galaxies on device
  double *dev_x1, *dev_y1, *dev_x2, *dev_y2;
  double *dev_xL, *dev_yL, *dev_zL;
  double *dev_e11, *dev_e21, *dev_e12, *dev_e22, *dev_w1, *dev_w2;

  // Numbers of sources and lenses
  int N1, N2, NL;
  N1 = x1.size(); // Number of galaxies
  N2 = x2.size();
  NL = xL.size();

  // Allocate memory on device
  CUDA_SAFE_CALL(cudaMalloc(&dev_x1, N1 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_y1, N1 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_e11, N1 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_e21, N1 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_w1, N1 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_x2, N2 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_y2, N2 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_e12, N2 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_e22, N2 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_w2, N2 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_xL, NL * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_yL, NL * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_zL, NL * sizeof(double)));

  // Copy values
  CUDA_SAFE_CALL(cudaMemcpy(dev_x1, x1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_y1, y1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_e11, e11.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_e21, e21.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_w1, w1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_x2, x2.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_y2, y2.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_e12, e12.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_e22, e22.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_w2, w2.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_xL, xL.data(), NL * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_yL, yL.data(), NL * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_zL, zL.data(), NL * sizeof(double), cudaMemcpyHostToDevice));

  // Calculate Gtilde

  // Declare container for Gplus, Gminus and weight on device
  double *dev_Gplus_real, *dev_Gplus_imag, *dev_Gminus_real, *dev_Gminus_imag, *dev_weight;

  // Allocate memory for Gplus, Gminus, and weight on device
  CUDA_SAFE_CALL(cudaMalloc(&dev_Gplus_real, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_Gplus_imag, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_Gminus_real, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_Gminus_imag, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_weight, num_bins * num_bins * num_bins * sizeof(double)));

  // Set Gtilde to 0 on device
  CUDA_SAFE_CALL(cudaMemset(dev_Gplus_real, 0, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMemset(dev_Gplus_imag, 0, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMemset(dev_Gminus_real, 0, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMemset(dev_Gminus_imag, 0, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMemset(dev_weight, 0, num_bins * num_bins * num_bins * sizeof(double)));

  g3lcong::addToGplusGminus<<<BLOCKS, THREADS>>>(dev_xL, dev_yL, dev_x1, dev_y1, dev_x2, dev_y2, dev_e11, dev_e21, dev_e12,
                                                 dev_e22, dev_w1, dev_w2, num_bins, NL, N1, N2, theta_min, theta_max,
                                                 dev_Gplus_real, dev_Gplus_imag, dev_Gminus_real, dev_Gminus_imag, dev_weight);

  // Declare arrays for Gplus, Gminus and weight on host
  double *Gplus_real, *Gplus_imag, *Gminus_real, *Gminus_imag, *weight;

  // Allocate memory for Gplus, Gminus and weight on host
  Gplus_real = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
  Gplus_imag = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
  Gminus_real = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
  Gminus_imag = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
  weight = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));

  // Copy Gplus and Gminus from device to host
  cudaMemcpy(Gplus_real, dev_Gplus_real, num_bins * num_bins * num_bins * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(Gplus_imag, dev_Gplus_imag, num_bins * num_bins * num_bins * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(Gminus_real, dev_Gminus_real, num_bins * num_bins * num_bins * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(Gminus_imag, dev_Gminus_imag, num_bins * num_bins * num_bins * sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy(weight, dev_weight, num_bins * num_bins * num_bins * sizeof(double), cudaMemcpyDeviceToHost);

  // Output

  // Print out vartheta1, vartheta2, psi, binsize, Gplus Gminus
  double phi_binsize = (phi_max - phi_min) / num_bins;
  double theta_binsize = log(theta_max / theta_min) / num_bins;

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

        int index = i * num_bins * num_bins + j * num_bins + k;

        // Weight
        double weight_ = weight[index];
        // Gplus real
        double Gplus_real_ = Gplus_real[index];
        // Gplus imag
        double Gplus_imag_ = Gplus_imag[index];
        // Gminus real
        double Gminus_real_ = Gminus_real[index];
        // Gminus imag
        double Gminus_imag_ = Gminus_imag[index];

        if (weight_ != 0)
        {
          Gplus_real_ /= weight_;
          Gplus_imag_ /= weight_;
          Gminus_real_ /= weight_;
          Gminus_imag_ /= weight_;
        }
        else
        {
          Gplus_real_ = 0;
          Gplus_imag_ = 0;
          Gminus_real_ = 0;
          Gminus_imag_ = 0;
        }

        // Output
        std::cout
            << theta1 << " "      // bin center theta 1 [arcmin]
            << theta2 << " "      // bin center theta 2 [arcmin]
            << phi << " "         // phi center [radians]
            << deltaTheta1 << " " // bin size theta 1[arcmin]
            << deltaTheta2 << " " // bin size theta 2[arcmin]
            << phi_binsize << " " // phi bin size [radians]
            << Gplus_real << " "  // Real part of Gplus [dimensionless]
            << Gplus_imag << " "  // Imaginary part of Gtilde [dimensionless]
            << Gminus_real << " " // Real part of Gplus [dimensionless]
            << Gminus_imag << " " // Imaginary part of Gtilde [dimensionless]
            << weight             // Weight of Gtilde [dimensionless]
            << std::endl;
      };
    };
  };

  return 0;
}
