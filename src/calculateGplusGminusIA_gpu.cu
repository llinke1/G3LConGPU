#include "helpers.h"
#include "constants.h"
#include "Function.h"
#include "kernelfunctions.cuh"
#include "cudaHelpers.cuh"
#include <string>
#include <iostream>

/**
 * Program computing GplusGminus_IA for one tile using GPU parallelization
 * The code requires a maximum distance pi between lenses and sources.
 * If this distance is too large, the estimator will give a sum of the Gplus and Gminus due to IA and the Gtilde due to cosmic shear
 * The x,y and z positions of the galaxies and pi need to be given in the same (physical!) units, e.g. Mpc/h!
 * The input files must be ASCII and formatted as
 * x y e1 e2 z
 *
 * Computation is done brute-force on submatrices of lens-source1-source2 space
 *
 * Output has logarithmic binning in r1 and r2 and linear in phi.
 * Units of r1 and r2 are the same as the x,y,z positions (e.g. Mpc/h)
 * phi is binned linearly
 *
 * @warning Number of parallel threads, blocks and size of submatrices is hard coded in Precompiler flags and needs to be fit to GPU
 *
 * Usage: calculateGplusGminusIA_gpu.x filename_lenses filename_shapes1 filename_shapes2 rMin rMax Nbins Pi
 *
 * Example: calculateGplusGminusIA_gpu.x tile1.lenses.dat tile1.shapes.dat tile1.shapes.dat 0.1 400 128 10
 *
 * @author Laila Linke, laila.linke@uibk.ac.at
 */

int main(int argc, char *argv[])
{
  // Checking command line
  int n_params = 7;
  std::string usage = "calculateGplusGminusIA_gpu.x  filename_lenses filename_shapes1 filename_shapes2 rMin rMax Nbins Pi";

  std::string example = "calculateGplusGminusIA_gpu.x tile1.lenses.dat  tile1.shapes.dat tile1.shapes.dat 0.1 400 128 10";

  //   Check Number of CMD Line arguments
  g3lcong::checkCmdLine(argc, n_params, usage, example);

  // Read in command line arguments
  std::string filename_lenses = argv[1];   // File with Lens Galaxies
  std::string filename_sources1 = argv[2]; // File with Shape1 Galaxies
  std::string filename_sources2 = argv[3]; // File with Shape2 Galaxies

  double r_min = std::stod(argv[4]); // Min R for calc of Gtilde (same units as X and Y positions of galaxies)
  double r_max = std::stod(argv[5]); // Max R for calc of Gtilde (same untis as X and Y positions of galaxies)
  int num_bins = std::stoi(argv[6]); // Number of Bins for Gtilde on ONE axis

  double Pi = std::stod(argv[7]); // Max line-of-sight distance between lenses and sources

  double phi_min = 0.0;      // Min Phi for calc of Gtilde [radians]
  double phi_max = 2 * M_PI; // Max Phi for calc of Gtilde [radians]

  // Reading in galaxies and copying to device

  // x,y, z vectors
  std::vector<double> x1, y1, z1, x2, y2, z2, xL, yL, zL, e11, e21, e12, e22, tmp;

  if (g3lcong::readSources2Dev(filename_sources1, 5, 1, 2, 3, 4, 5, x1,
                               y1, e11, e21, z1))
    return 1;


  if (g3lcong::readSources2Dev(filename_sources2, 5, 1, 2, 3, 4, 5, x2,
                               y2, e12, e22, z2))
    
    return 1;

  if (g3lcong::readLenses2Dev(filename_lenses, 5, 1, 2, 5, xL, yL, zL))
    return 1;

  
  // Declare arrays for coordinates of galaxies on device
  double *dev_x1, *dev_y1, *dev_z1, *dev_x2, *dev_y2, *dev_z2;
  double *dev_xL, *dev_yL, *dev_zL;
  double *dev_e11, *dev_e21, *dev_e12, *dev_e22;

  // Numbers of sources and lenses
  int N1, N2, NL;
  N1 = x1.size(); // Number of galaxies
  N2 = x2.size();
  NL = xL.size();

  // Allocate memory on device
  CUDA_SAFE_CALL(cudaMalloc(&dev_x1, N1 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_y1, N1 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_z1, N1 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_e11, N1 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_e21, N1 * sizeof(double)));

  CUDA_SAFE_CALL(cudaMalloc(&dev_x2, N2 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_y2, N2 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_z2, N2 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_e12, N2 * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_e22, N2 * sizeof(double)));

  CUDA_SAFE_CALL(cudaMalloc(&dev_xL, NL * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_yL, NL * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_zL, NL * sizeof(double)));

  // Copy values
  CUDA_SAFE_CALL(cudaMemcpy(dev_x1, x1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_y1, y1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_e11, e11.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_e21, e21.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_z1, z1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_x2, x2.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_y2, y2.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_e12, e12.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_e22, e22.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_z2, z2.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_xL, xL.data(), NL * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_yL, yL.data(), NL * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(dev_zL, zL.data(), NL * sizeof(double), cudaMemcpyHostToDevice));

  // Declare container for Gplus, Gminus and weight on device
  double *dev_Gplus_real, *dev_Gplus_imag, *dev_Gminus_real, *dev_Gminus_imag, *dev_weight;

  // Allocate memory for Gplus, Gminus, and weight on device
  CUDA_SAFE_CALL(cudaMalloc(&dev_Gplus_imag, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_Gplus_real, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_Gminus_real, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_Gminus_imag, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMalloc(&dev_weight, num_bins * num_bins * num_bins * sizeof(double)));

  // Set Gtilde to 0 on device
  CUDA_SAFE_CALL(cudaMemset(dev_Gplus_real, 0, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMemset(dev_Gplus_imag, 0, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMemset(dev_Gminus_real, 0, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMemset(dev_Gminus_imag, 0, num_bins * num_bins * num_bins * sizeof(double)));
  CUDA_SAFE_CALL(cudaMemset(dev_weight, 0, num_bins * num_bins * num_bins * sizeof(double)));

  g3lcong::addToGplusGminusIA<<<BLOCKS, THREADS>>>(dev_xL, dev_yL, dev_zL,
                                                   dev_x1, dev_y1, dev_z1,
                                                   dev_x2, dev_y2, dev_z2,
                                                   dev_e11, dev_e21, dev_e12, dev_e22, 
                                                   Pi, NL, N1, N2, r_min, r_max, num_bins,
                                                   dev_Gplus_real, dev_Gplus_imag, dev_Gminus_real, dev_Gminus_imag, dev_weight);


  CUDA_SAFE_CALL(cudaPeekAtLastError());
  CUDA_SAFE_CALL(cudaDeviceSynchronize());

  // Declare arrays for Gplus, Gminus and weight on host
  double *Gplus_real, *Gplus_imag, *Gminus_real, *Gminus_imag, *weight;

  // Allocate memory for Gplus, Gminus and weight on host
  Gplus_real = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
  Gplus_imag = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
  Gminus_real = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
  Gminus_imag = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
  weight = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));

  // Copy Gplus and Gminus from device to host
  CUDA_SAFE_CALL(cudaMemcpy(Gplus_real, dev_Gplus_real, num_bins * num_bins * num_bins * sizeof(double),
             cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Gplus_imag, dev_Gplus_imag, num_bins * num_bins * num_bins * sizeof(double),
             cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Gminus_real, dev_Gminus_real, num_bins * num_bins * num_bins * sizeof(double),
             cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Gminus_imag, dev_Gminus_imag, num_bins * num_bins * num_bins * sizeof(double),
             cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(weight, dev_weight, num_bins * num_bins * num_bins * sizeof(double),
             cudaMemcpyDeviceToHost));
		

  // Output

  // Print out vartheta1, vartheta2, psi, binsize, Gplus Gminus
  double phi_binsize = (phi_max - phi_min) / num_bins;
  double r_binsize = log(r_max / r_min) / num_bins;

  std::vector<double> bin_centers_theta;
  std::vector<double> newbin_edges_theta;
  for (int i = 0; i < num_bins; i++)
  {
    bin_centers_theta.push_back(0.5 * (exp(log(r_min) + r_binsize * (i + 1)) + exp(log(r_min) + r_binsize * i)));
  }

  newbin_edges_theta.push_back(r_min);
  for (int i = 0; i < num_bins - 1; i++)
  {
    newbin_edges_theta.push_back(0.5 * (bin_centers_theta[i + 1] + bin_centers_theta[i]));
  }
  newbin_edges_theta.push_back(r_max);

  for (int i = 0; i < num_bins; i++)
  {
    // Theta1 Center of this bin
    double theta1 = bin_centers_theta[i];
    // Theta1 Binsize of this bin
    double deltaTheta1 = newbin_edges_theta[i + 1] - newbin_edges_theta[i];

    for (int j = 0; j < num_bins; j++)
    {
      // Theta2 Center of this bin
      double theta2 = bin_centers_theta[j];

      // Theta2 Binsize of this bin
      double deltaTheta2 = newbin_edges_theta[j + 1] - newbin_edges_theta[j];

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
            << theta1 << " "       // bin center theta 1 [arcmin]
            << theta2 << " "       // bin center theta 2 [arcmin]
            << phi << " "          // phi center [radians]
            << deltaTheta1 << " "  // bin size theta 1[arcmin]
            << deltaTheta2 << " "  // bin size theta 2[arcmin]
            << phi_binsize << " "  // phi bin size [radians]
            << Gplus_real_ << " "  // Real part of Gplus [dimensionless]
            << Gplus_imag_ << " "  // Imaginary part of Gtilde [dimensionless]
            << Gminus_real_ << " " // Real part of Gplus [dimensionless]
            << Gminus_imag_ << " " // Imaginary part of Gtilde [dimensionless]
            << weight_             // Weight of Gtilde [dimensionless]
            << std::endl;
      };
    };
  };

  return 0;
}
