#include "helpers.h"
#include "constants.h"
#include "cuda_helpers.cuh"
#include "kernelfunctions.cuh"
#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
    // Checking Command Line

    int n_params = 6;       // Expected number of params
    std::string usage = ""; // Usage description

    std::string example = ""; // Example usage

    // Check Number of CMD Line arguments
    g3lcong::checkCmdLine(argc, n_params, usage, example);

    // Read in command line arguments
    std::string filename_sources = argv[1]; // File with Source Galaxies
    std::string filename_lenses = argv[2];  // File with Lens Galaxies

    double theta_min = std::stod(argv[3]); // Min Theta for calc of Gtilde [arcmin]
    double theta_max = std::stod(argv[4]); // Max Theta for calc of Gtilde [arcmin]
    int num_bins = std::stoi(argv[5]);     // Number of Bins for Gtilde on ONE axis

    double deltaChi = std::stod(argv[6]);

    // Reading in galaxies and copying to device

    // x,y, z vectors
    std::vector<double> x1, y1, z1, chi1, xS, yS, e1, e2, w, chiS;

    if (g3lcong::readPhysSources(filename_sources, 7, 1, 2, 3, 4, 6, 7, xS, yS, e1, e2, w, chiS))
        return 1;
    if (g3lcong::readPhysLenses(filename_lenses, 7, 1, 2, 7, x1, y1, chi1))
        return 1;

    // Declare arrays for coordinates of galaxies on device
    double *dev_x1, *dev_y1, *dev_chi1;
    double *dev_xS, *dev_yS, *dev_e1, *dev_e2, *dev_w, *dev_chiS;

    // Numbers of sources and lenses
    int N1, NS;
    N1 = x1.size(); // Number of galaxies

    NS = xS.size();

    // Allocate memory on device
    CUDA_SAFE_CALL(cudaMalloc(&dev_x1, N1 * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_y1, N1 * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_chi1, N1 * sizeof(double)));

    CUDA_SAFE_CALL(cudaMalloc(&dev_xS, NS * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_yS, NS * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_e1, NS * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_e2, NS * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_w, NS * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_chiS, NS * sizeof(double)));

    // Copy values
    CUDA_SAFE_CALL(cudaMemcpy(dev_x1, x1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_y1, y1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_chi1, chi1.data(), N1 * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_xS, xS.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_yS, yS.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_e1, e1.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_e2, e2.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_w, w.data(), NS * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(dev_chiS, chiS.data(), NS * sizeof(double), cudaMemcpyHostToDevice));

    // Calculate Gtilde

    // Declare container for Greal, Gimag and weight on device
    double *dev_Greal, *dev_Gimag, *dev_weight;

    // Allocate memory for Greal, Gimag and weight on device
    CUDA_SAFE_CALL(cudaMalloc(&dev_Greal, num_bins * num_bins * num_bins * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_Gimag, num_bins * num_bins * num_bins * sizeof(double)));
    CUDA_SAFE_CALL(cudaMalloc(&dev_weight, num_bins * num_bins * num_bins * sizeof(double)));

    // Check if memory for Gtilde could be allocated on device
    if (0 == dev_Greal || 0 == dev_Gimag || 0 == dev_weight)
    {
        std::cout << "Couldnt allocate Space for Gtilde" << std::endl;
        return 1;
    };

    // Set Gtilde to 0 on device
    CUDA_SAFE_CALL(cudaMemset(dev_Greal, 0, num_bins * num_bins * num_bins * sizeof(double)));
    CUDA_SAFE_CALL(cudaMemset(dev_Gimag, 0, num_bins * num_bins * num_bins * sizeof(double)));
    CUDA_SAFE_CALL(cudaMemset(dev_weight, 0, num_bins * num_bins * num_bins * sizeof(double)));

    g3lcong::addToGE<<<BLOCKS, THREADS>>>(dev_x1, dev_y1,
                                          dev_xS, dev_yS,
                                          dev_chi1, dev_chiS,
                                          dev_e1, dev_e2,
                                          dev_w, deltaChi,
                                          num_bins, N1, NS,
                                          theta_min, theta_max,
                                          dev_Greal,
                                          dev_Gimag, dev_weight);

    // Declare arrays for total Greal, Gimag and weight on host
    double *Greal_tot, *Gimag_tot, *weight_tot;

    // Allocate memory for total Greal, Gimag and weight on host
    Greal_tot = (double *)malloc(num_bins * sizeof(double));
    Gimag_tot = (double *)malloc(num_bins * sizeof(double));
    weight_tot = (double *)malloc(num_bins * sizeof(double));

    if (Greal_tot == NULL || Gimag_tot == NULL || weight_tot == NULL)
    {
        std::cerr << "calculateGtilde_gpu: Couldn't allocate memory" << std::endl;
        exit(1);
    };

    // Copy Gtilde from device to host
    CUDA_SAFE_CALL(cudaMemcpy(Greal_tot, dev_Greal, num_bins * sizeof(double),
               cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(Gimag_tot, dev_Gimag, num_bins * sizeof(double),
               cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(weight_tot, dev_weight, num_bins * sizeof(double),
               cudaMemcpyDeviceToHost));

    // Output

    // Print out vartheta1, vartheta2, psi, binsize, Gtilde
    double theta_binsize = log(theta_max / theta_min) / num_bins;

    for (int i = 0; i < num_bins; i++)
    {
        // Theta1 Center of this bin
        double theta = 0.5 * (exp(log(theta_min) + theta_binsize * (i + 1)) + exp(log(theta_min) + theta_binsize * i));
        // Theta1 Binsize of this bin
        double deltaTheta = exp(log(theta_min) + theta_binsize * (i + 1)) - exp(log(theta_min) + theta_binsize * i);

        // Weight
        double weight = weight_tot[i];
        // Greal
        double Greal = Greal_tot[i];
        // Gimag
        double Gimag = Gimag_tot[i];

        if (weight != 0)
        {
            Greal /= weight;
            Gimag /= weight;
        }
        else
        {
            Greal = 0;
            Gimag = 0;
        }

        // Output
        std::cout
            << theta << " "      // bin center theta 1 [Mpc/h]
            << deltaTheta << " " // bin size theta 1[arcmin or Mpc]
            << Greal << " "      // Real part of Gtilde [dimensionless or Msun/Mpc²]
            << Gimag << " "      // Imaginary part of Gtilde [dimensionless or Msun/Mpc²]
            << weight            // Weight of Gtilde [dimensionless or Msun/Mpc²]
            << std::endl;
    };

    return 0;
}