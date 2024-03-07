#include "helpers.h"
#include "constants.h"
#include "Function.h"
#include "kernelfunctions.cuh"
#include "cudaHelpers.cuh"
#include <string>
#include <iostream>

/**
 * Program computing Gtilde_IA_biased for one tile using GPU parallelization
 * The output is not corrected for the second-order correlation between lens galaxies
 * The code requires a maximum distance pi between lenses and sources.
 * If this distance is too large, the estimator will give a sum of the Gtilde due to IA and the Gtilde due to cosmic shear
 * The x,y and z positions of the galaxies and pi need to be given in the same (physical!) units, e.g. Mpc/h!
 * The input files must be ASCII and formatted as
 * x y e1 e2 z
 *
 * Computation is done brute-force on submatrices of lens1-lens2-source space
 *
 * Output has logarithmic binning in r1 and r2 and linear in phi.
 * Units of r1 and r2 are the same as the x,y,z positions (e.g. Mpc/h)
 * phi is binned linearly
 *
 * @warning Number of parallel threads, blocks and size of submatrices is hard coded in Precompiler flags and needs to be fit to GPU
 *
 * Usage: calculateGtildeIABiased_gpu.x filename_lenses1 filename_lenses2 filename_shapes rMin rMax Nbins Pi
 *
 * Example: calculateGtildeIABiased_gpu.x tile1.lenses.dat tile1.lenses.dat tile1.shapes.dat 0.1 400 128 10
 *
 * @author Laila Linke, laila.linke@uibk.ac.at
 */

int main(int argc, char *argv[])
{
    // Checking command line
    int n_params = 7;
    std::string usage = "calculateGtildeIABiased_gpu.x  filename_shapes filename_lenses1 filename_lenses2 rMin rMax Nbins Pi";

    std::string example = "calculateGtildeIABiased_gpu.x tile1.shapes.dat  tile1.lenses.dat tile1.lenses.dat 0.1 400 128 10";

    //   Check Number of CMD Line arguments
    g3lcong::checkCmdLine(argc, n_params, usage, example);

    // Read in command line arguments
    std::string filename_sources = argv[1]; // File with Source Galaxies
    std::string filename_lenses1 = argv[2]; // File with Lens1 Galaxies
    std::string filename_lenses2 = argv[3]; // File with Lens2 Galaxies

    double r_min = std::stod(argv[4]); // Min R for calc of Gtilde (same units as X and Y positions of galaxies)
    double r_max = std::stod(argv[5]); // Max R for calc of Gtilde (same untis as X and Y positions of galaxies)
    int num_bins = std::stoi(argv[6]); // Number of Bins for Gtilde on ONE axis

    double Pi = std::stod(argv[7]); // Max line-of-sight distance between lenses and sources

    double phi_min = 0.0;      // Min Phi for calc of Gtilde [radians]
    double phi_max = 2 * M_PI; // Max Phi for calc of Gtilde [radians]

    // Reading in galaxies and copying to device

    std::vector<double> x1, y1, z1, x2, y2, z2, xS, yS, zS, e1, e2, tmp;

    // Reading in sources
    if (g3lcong::readSources2Dev(filename_sources, 5, 1, 2, 3, 4, 5, xS, yS, e1, e2, zS))
        return 1;

    // Reading in lenses
    if (g3lcong::readLenses2Dev(filename_lenses1, 5, 1, 2, 5, x1, y1, z1))
        return 1;

    if (g3lcong::readLenses2Dev(filename_lenses2, 5, 1, 2, 5, x2, y2, z2))
        return 1;

    // Declare arrays for coordinates of galaxies on device
    double *dev_x1, *dev_y1, *dev_x2, *dev_y2, *dev_z1, *dev_z2, *dev_zS;
    double *dev_xS, *dev_yS, *dev_e1, *dev_e2;

    // Numbers of sources and lenses
    int N1 = x1.size(); // Number of lens 1galaxies
    int N2 = x2.size(); // Number of lens2 galaxies
    int NS = xS.size(); // Number of source galaxies

    // Allocate memory on device
    CUDA_SAFE_CALL(cudaMalloc(&dev_x1, N1 * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_y1, N1 * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_z1, N1 * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_x2, N2 * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_y2, N2 * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_z2, N2 * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_xS, NS * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_yS, NS * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_zS, NS * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_e1, NS * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_e2, NS * sizeof(double)));

    // Copy values
    CUDA_SAFE_CALL(cudaMemcpy(dev_x1, x1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_y1, y1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_z1, z1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_x2, x2.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_y2, y2.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_z2, z2.data(), N2 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_xS, xS.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_yS, yS.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_e1, e1.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_e2, e2.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_zS, zS.data(), NS * sizeof(double), cudaMemcpyHostToDevice));

    // Declare container for Greal, Gimag and weight on device
    double *dev_Greal, *dev_Gimag, *dev_weight;

    // Allocate memory for Greal, Gimag and weight on device
    CUDA_SAFE_CALL(cudaMalloc(&dev_Greal, num_bins * num_bins * num_bins * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_Gimag, num_bins * num_bins * num_bins * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_weight, num_bins * num_bins * num_bins * sizeof(double)));

    // Set Gtilde to 0 on device
    CUDA_SAFE_CALL(cudaMemset(dev_Greal, 0, num_bins * num_bins * num_bins * sizeof(double)));
    CUDA_SAFE_CALL(cudaMemset(dev_Gimag, 0, num_bins * num_bins * num_bins * sizeof(double)));
    CUDA_SAFE_CALL(cudaMemset(dev_weight, 0, num_bins * num_bins * num_bins * sizeof(double)));

    // Calculate Gtilde
    g3lcong::addToGtildeIABiased<<<BLOCKS, THREADS>>>(dev_x1, dev_y1, dev_z1, dev_x2, dev_y2, dev_z2,
                                                      dev_xS, dev_yS, dev_zS, dev_e1, dev_e2,
                                                      Pi, N1, N2, NS, r_min, r_max, num_bins, dev_Greal, dev_Gimag, dev_weight);

    // Declare arrays for Greal, Gimag and weight on host
    double *Greal, *Gimag, *weight;

    // Allocate memory for total Greal, Gimag and weight on host
    Greal = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
    Gimag = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));
    weight = (double *)malloc(num_bins * num_bins * num_bins * sizeof(double));

    // Copy Gtilde from device to host
    CUDA_SAFE_CALL(cudaMemcpy(Greal, dev_Greal, num_bins * num_bins * num_bins * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(Gimag, dev_Gimag, num_bins * num_bins * num_bins * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(weight, dev_weight, num_bins * num_bins * num_bins * sizeof(double), cudaMemcpyDeviceToHost));

    // Output

    double phi_binsize = (phi_max - phi_min) / num_bins; // size of phi bins
    double r_binsize = log(r_max / r_min) / num_bins;    // log size of r bins

    for (int i = 0; i < num_bins; i++)
    {
        // Theta1 Center of this bin
        double theta1 = 0.5 * (exp(log(r_min) + r_binsize * (i + 1)) + exp(log(r_min) + r_binsize * i));
        // Theta1 Binsize of this bin
        double deltaTheta1 = exp(log(r_min) + r_binsize * (i + 1)) - exp(log(r_min) + r_binsize * i);

        for (int j = 0; j < num_bins; j++)
        {
            // Theta2 Center of this bin
            double theta2 = 0.5 * (exp(log(r_min) + r_binsize * (j + 1)) + exp(log(r_min) + r_binsize * j));
            // Theta2 Binsize of this bin
            double deltaTheta2 = exp(log(r_min) + r_binsize * (j + 1)) - exp(log(r_min) + r_binsize * j);

            for (int k = 0; k < num_bins; k++)
            {
                // Phi Center of this bin
                double phi = (k + 0.5) * phi_binsize + phi_min;

                int index = i * num_bins * num_bins + j * num_bins + k;

                // Weight
                double w = weight[index];
                // Greal
                double Gr = Greal[index];
                // Gimag
                double Gi = Gimag[index];

                if (w != 0)
                {
                    Gr /= w;
                    Gi /= w;
                }
                else
                {
                    Gr = 0;
                    Gi = 0;
                }

                // Output
                std::cout
                    << theta1 << " "      // bin center theta 1 [arcmin or Mpc]
                    << theta2 << " "      // bin center theta 2 [arcmin or Mpc]
                    << phi << " "         // phi center [radians]
                    << deltaTheta1 << " " // bin size theta 1[arcmin or Mpc]
                    << deltaTheta2 << " " // bin size theta 2[arcmin or Mpc]
                    << phi_binsize << " " // phi bin size [radians]
                    << Gr << " "          // Real part of Gtilde [dimensionless or Msun/Mpc²]
                    << Gi << " "          // Imaginary part of Gtilde [dimensionless or Msun/Mpc²]
                    << w                  // Weight of Gtilde [dimensionless or Msun/Mpc²]
                    << std::endl;
            };
        };
    };

    return 0;
}