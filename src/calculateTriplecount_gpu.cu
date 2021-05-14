#include "helpers.h"
#include "constants.h"
#include "Function.h"
#include "kernelfunctions.cuh"
#include <string>
#include <iostream>


int main(int argc, char* argv[])
{
  // Checking Command line
 int n_params=7; // Expected number of params
  std::string usage="calculateTriplecount_gpu.x filename_gal1 filename_gal2 filename_gal3 r_min r_max num_bins filename_Dcomov"; //Usage description

  std::string example=""; //Example usage
  
  // Check Number of CMD Line arguments
  g3lcong::checkCmdLine(argc, n_params, usage, example);


  // Read in command line arguments
  std::string filename_gal1 = argv[1]; // File with Galaxies 1
  std::string filename_gal2 = argv[2]; // File with Galaxies 2
  std::string filename_gal3 = argv[3]; // File with Galaxies 3

  double r_min=std::stod(argv[4]); // Min R for calc of triplecount [Mpc]
  double r_max=std::stod(argv[5]); // Max R for calc of triplecount [Mpc]
  int num_bins=std::stoi(argv[6]); // Number of Bins for triplecount on ONE axis

  // File with angular diameter distance, if "none": distance is 1 [Mpc]
  std::string filename_com_distance = argv[7];

  std::cerr<<"Finished reading CMD Line"<<std::endl;
  
  // Reading in galaxies and copying to device

  // x,y, z vectors
  std::vector<double> x1, y1, z1, x2, y2, z2, x3, y3, z3;
  

  if(g3lcong::readLenses2Dev(filename_gal1, 3, 1, 2, 3, x1, y1, z1)) return 1;
  if(g3lcong::readLenses2Dev(filename_gal2, 3, 1, 2, 3, x2, y2, z2)) return 1;
  if(g3lcong::readLenses2Dev(filename_gal3, 3, 1, 2, 3, x3, y3, z3)) return 1;
  

 //Declare arrays for coordinates of galaxies on device
  double *dev_x1, *dev_y1, *dev_z1, *dev_x2, *dev_y2, *dev_z2, *dev_x3, *dev_y3, *dev_z3;
  
  //Numbers of galaxies
  int N1, N2, N3;
  N1=x1.size(); //Number of galaxies
  N2=x2.size(); 
  N3=x3.size();
  
  // Allocate memory on device
  cudaError_t err1 = cudaMalloc(&dev_x1, N1*sizeof(double));
  cudaError_t err2 = cudaMalloc(&dev_y1, N1*sizeof(double));
  cudaError_t err3 = cudaMalloc(&dev_z1, N1*sizeof(double));

  if(err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess)
    {
      std::cerr<<"Could not allocate memory"<<std::endl;
      exit(1);
    };

  err1 = cudaMalloc(&dev_x2, N2*sizeof(double));
  err2 = cudaMalloc(&dev_y2, N2*sizeof(double));
  err3 = cudaMalloc(&dev_z2, N2*sizeof(double));

  if(err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess)
    {
      std::cerr<<"Could not allocate memory"<<std::endl;
      exit(1);
    };

  err1 = cudaMalloc(&dev_x3, N3*sizeof(double));
  err2 = cudaMalloc(&dev_y3, N3*sizeof(double));
  err3 = cudaMalloc(&dev_z3, N3*sizeof(double));
 

  if(err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess)
    {
      std::cerr<<"Could not allocate memory"<<std::endl;
      exit(1);
    };

  // Copy values
  cudaMemcpy(dev_x1, x1.data(), N1*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_y1, y1.data(), N1*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_z1, z1.data(), N1*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_x2, x2.data(), N2*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_y2, y2.data(), N2*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_z2, z2.data(), N2*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_x3, x3.data(), N3*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_y3, y3.data(), N3*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_z3, z3.data(), N3*sizeof(double), cudaMemcpyHostToDevice);

  // Reading in Comoving Distance
 
  //Declare array for comoving distance
  std::vector<double> com_distance;
  double* dev_com_distance=NULL;

  double z_min, z_max;
  
  // Try reading and stop if not successful
  if(g3lcong::readFunction2Dev(filename_com_distance, num_bins,
			       z_min, z_max, com_distance))
    return 1;


  err1 = cudaMalloc(&dev_com_distance, num_bins*sizeof(double));
  if(err1 != cudaSuccess)
    {
      std::cerr<<"Could not allocate memory"<<std::endl;
      exit(1);
    };
  cudaMemcpy(dev_com_distance, com_distance.data(), num_bins*sizeof(double), cudaMemcpyHostToDevice);


  std::cerr<<"Finished Copying to Device"<<std::endl;
  
  // Calculate Triplecount

  // Declare container for Triplecount
  int* dev_triplecount;

  // Allocate memory for Triplecount on device
  err1 = cudaMalloc(&dev_triplecount, num_bins*num_bins*num_bins*num_bins*num_bins*sizeof(int));
  if(err1 != cudaSuccess)
    {
      std::cerr<<"Could not allocate memory"<<std::endl;
      exit(1);
    };

  // Set Triplecount to 0
  cudaMemset(dev_triplecount, 0, num_bins*num_bins*num_bins*num_bins*num_bins*sizeof(int));

  double r_binwidth=log(r_max/r_min)/num_bins;
  double z_binwidth=(z_max-z_min)/num_bins;
  std::cerr<<"Started Calculating Triplecount"<<std::endl;
  
  g3lcong::addToTriplecount<<<BLOCKS, THREADS>>>(dev_x1, dev_y1, dev_z1, dev_x2, dev_y2, dev_z2, dev_x3, dev_y3, dev_z3, N1, N2, N3, num_bins, r_min, r_binwidth,  dev_com_distance, z_min, z_binwidth, dev_triplecount); 

  std::cerr<<"Finished calculating Triplecount"<<std::endl;
 // Free memory on device
  cudaFree(dev_x1);
  cudaFree(dev_y1);
  cudaFree(dev_z1);
  cudaFree(dev_x2);
  cudaFree(dev_y2);
  cudaFree(dev_z2);
  cudaFree(dev_x3);
  cudaFree(dev_y3);
  cudaFree(dev_z3);
  cudaFree(dev_com_distance);
  
  // Declare arry for Triplecount on host
  int* triplecount;

  // Allocate memory for Triplecount on host
  triplecount = (int*) malloc(num_bins*num_bins*num_bins*num_bins*num_bins*sizeof(int));

  //Copy Triplecount from device to host
  cudaMemcpy(triplecount, dev_triplecount, num_bins*num_bins*num_bins*num_bins*num_bins*sizeof(int),
	     cudaMemcpyDeviceToHost);


  cudaFree(dev_triplecount);

  std::cerr<<"Started output"<<std::endl;

  // Output
  for(int i=0; i<num_bins; i++)
    {
      double r12=exp(log(r_min)+i*r_binwidth);
      for(int j=0; j<num_bins; j++)
	{
	  double r13=exp(log(r_min)+j*r_binwidth);
	  for(int k=0; k<num_bins; k++)
	    {
	      double r23=exp(log(r_min)+k*r_binwidth);
	      for(int l=0; l<num_bins; l++)
		{
		  double pi12=exp(log(r_min)+l*r_binwidth);
		  for(int m=0; m<num_bins; m++)
		    {
		      double pi23=exp(log(r_min)+m*r_binwidth);
		      int index=i*num_bins*num_bins*num_bins*num_bins
			+j*num_bins*num_bins*num_bins
			+k*num_bins*num_bins
			+l*num_bins
			+m;
		      std::cout<<r12<<" "
			       <<r13<<" "
			       <<r23<<" "
			       <<pi12<<" "
			       <<pi23<<" "
			       <<triplecount[index]<<" "
			       <<std::endl;
		    };
		};
	    };
	};
    };
  std::cerr<<"Finished output"<<std::endl;
  return 0;
}
