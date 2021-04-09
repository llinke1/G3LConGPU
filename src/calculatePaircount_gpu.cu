#include "kernelfunctions.cuh"
#include "helpers.h"
#include "constants.h"
#include <string>
#include <iostream>


/**
 * Program computing the paircount using GPU parallelization
 * Needs as input stddev of redshiftweighting
 * If stddev of redshiftweighting is set to 0, no redshift weighting is done
 *
 * Usage: calculatePaircount_gpu.x catalog1 catalog2 weight-z-file theta_min theta_max number_bins
 *
 * Example: calculatePaircount_gpi.x ../../Data/KIDS_129.0_-0.5.objects.dat
 * ../../Data/KIDS_129.0_-0.5.objects.dat none
 *
 * @author Laila Linke, llinke@astro.uni-bonn.de
 */

int main(int argc, char* argv[])
{
  int n_params=6; // Expected number of params
  std::string usage="calculatePaircount_gpu.x filename_galaxies1 filename_galaxies2 sigmaZ theta_min theta_max num_bins"; //Usage description

  std::string example="calculatePaircount_gpu.x ../../Data/KIDS_129.0_-0.5.objects.dat ../../Data/KIDS_129.0_-0.5.objects.dat 0 0.15 0.79 100"; //Example usage
  
  // Check Number of CMD Line arguments
  g3lcong::checkCmdLine(argc, n_params, usage, example);
  
   
  // Reading in commandline
  std::string filename1 = argv[1]; // File with first galaxy catalog
  std::string filename2 = argv[2]; // File with second galaxy catalog
  double sigmaZ = std::stod(argv[3]); //Width of redshift weighting Gaussian, if 0: no weighting
  double theta_min = std::stod(argv[4]); //[arcmin]
  double theta_max = std::stod(argv[5]); //[arcmin]
  int num_bins = std::stoi(argv[6]); 

  double theta_binsize = log(theta_max / theta_min)/num_bins;

   // Reading in galaxies and copying to device

  // x,y, z vectors
  std::vector<double> x1, y1, z1, x2, y2, z2;
  
  
  //Numbers of galaxies
  int N1, N2;

  if(g3lcong::readLenses2Dev(filename1, 6, 1, 2, 5, x1, y1, z1)) return 1;
  if(g3lcong::readLenses2Dev(filename2, 6, 1, 2, 5, x2, y2, z2)) return 1;
  

 //Declare arrays for coordinates of galaxies on device
  double *dev_x1, *dev_y1, *dev_x2, *dev_y2, *dev_z1, *dev_z2;
  N1=x1.size(); //Number of galaxies
  N2=x2.size(); 

  // Allocate memory on device
  cudaError_t err1 = cudaMalloc(&dev_x1, N1*sizeof(double));
  cudaError_t err2 = cudaMalloc(&dev_y1, N1*sizeof(double));
  cudaError_t err3 = cudaMalloc(&dev_z1, N1*sizeof(double));

  if(err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess)
    {
      std::cout<<"Could not allocate memory"<<std::endl;
      return 1;
    };

  err1 = cudaMalloc(&dev_x2, N2*sizeof(double));
  err2 = cudaMalloc(&dev_y2, N2*sizeof(double));
  err3 = cudaMalloc(&dev_z2, N2*sizeof(double));

  if(err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess)
    {
      std::cout<<"Could not allocate memory"<<std::endl;
      return 1;
    };

  // Copy values
  cudaMemcpy(dev_x1, x1.data(), N1*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_y1, y1.data(), N1*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_z1, z1.data(), N1*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_x2, x2.data(), N2*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_y2, y2.data(), N2*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_z2, z2.data(), N2*sizeof(double), cudaMemcpyHostToDevice);
  
  // Calculate Paircount

  // Declare container for paircount on device
  double * dev_paircount;

  // Allocate memory for paircount
  cudaMalloc(&dev_paircount, num_bins*sizeof(double));

  // Check if memory for paircount could be allocated
  if(0 == dev_paircount)
    {
      std::cout<<"calculatePaircount_gpu:: Couldn't allocate space for paircount"<<std::endl;
      return 1;
    };
  
  //Set values of paircount to 0 (on device)
  cudaMemset(dev_paircount, 0, num_bins*sizeof(double));
  
  double sigma2=sigmaZ*sigmaZ; //Square of stddev of z weighting


  g3lcong::addToPaircount<<< BLOCKS, THREADS>>>(dev_x1,dev_y1,dev_x2,dev_y2,
						dev_z1, dev_z2, N1, N2,
						num_bins,
						theta_min, theta_binsize,
						sigma2, dev_paircount);

  //Declare array for total paircount on host
  double *paircount;
      
  //Allocate memory for total paircount on host
  paircount=(double *) malloc(num_bins*sizeof(double));

  //Copy result for paircount from device to host
  cudaMemcpy(paircount, dev_paircount, num_bins*sizeof(double), cudaMemcpyDeviceToHost);

  // Calculate Normalization
  double normalization=0;
  for(int i=0; i<num_bins; i++)
    {
      normalization+=paircount[i];
    };

  // Print out Paircount
  for(unsigned int i=0; i<num_bins; i++)
    {
      std::cout<<theta_min*exp(i*theta_binsize)<<" "<<paircount[i]<<" "<<normalization<<std::endl;
    }

  
  return 0;
}
