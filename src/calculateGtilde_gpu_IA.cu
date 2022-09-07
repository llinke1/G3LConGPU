#include "helpers.h"
#include "constants.h"
#include "Function.h"
#include "kernelfunctions.cuh"
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

int main(int argc, char* argv[])
{
  // Checking Command Line

  int n_params=13; // Expected number of params
  std::string usage="calculateGtilde_gpu.x filename_sources filename_lenses1 filename_lenses2 filename_omega theta_min theta_max r_min r_max num_bins sigmaZ filename_sigma_crit Physical?"; //Usage description

  std::string example="calculateGtilde_gpu.x  ../../Data/KIDS_129.0_-0.5.sources.dat ../../Data/KIDS_129.0_-0.5.objects.dat ../../Data/KIDS_129.0_-0.5.objects.dat ../products/omega_allTiles/all.omega.dat  0.15 79.9 0.1 300 100 0 none 0"; //Example usage
  
  // Check Number of CMD Line arguments
  g3lcong::checkCmdLine(argc, n_params, usage, example);


  // Read in command line arguments
  std::string filename_sources = argv[1]; // File with Source Galaxies
  std::string filename_lenses1 = argv[2]; // File with Lens1 Galaxies
  std::string filename_lenses2 = argv[3]; // File with Lens2 Galaxies

  std::string filename_omega = argv[4]; // File with Angular Correlation Func
  
  double theta_min=std::stod(argv[5]); // Min Theta for calc of Gtilde [arcmin]
  double theta_max=std::stod(argv[6]); // Max Theta for calc of Gtilde [arcmin]
  double r_min=std::stod(argv[7]); // Min R for calc of Gtilde [Mpc] (ignored if not physical)
  double r_max=std::stod(argv[8]); // Max R for calc of Gtilde [Mpc] (ignored if not physical)
  int num_bins=std::stoi(argv[9]); // Number of Bins for Gtilde on ONE axis

  // Width of Redshift Weighting Gaussian, if 0:  no weighting
  double sigmaZ=std::stod(argv[10]);
  // File with source-averaged Sigma Crit, if "none": SigCrit is 1
  std::string filename_sigma_crit = argv[11];
  // File with angular diameter distance, if "none": distance is 1
  std::string filename_angular_distance = argv[12];

  bool physical=std::stoi(argv[13]); //Is 1 if physical Gtilde, 0 if angular Gtilde
  
  double phi_min=0.0; // Min Phi for calc of Gtilde [radians]
  double phi_max=2*g3lcong::pi; // Max Phi for calc of Gtilde [radians]

  // Reading in galaxies and copying to device

 
  // x,y, z vectors
  std::vector<double> x1, y1, z1, x2, y2, z2, xS, yS, zS, e1, e2, w, tmp;
  

  if(g3lcong::readSources2Dev(filename_sources, 6, 1, 2, 3, 4, 6, xS,
			      yS, e1, e2, w)) return 1;
  
  if(g3lcong::readLenses2Dev(filename_sources, 6, 1, 2, 5, tmp,
			      tmp, zS)) return 1;

  if(g3lcong::readLenses2Dev(filename_lenses1, 6, 1, 2, 5, x1, y1, z1)) return 1;
  if(g3lcong::readLenses2Dev(filename_lenses2, 6, 1, 2, 5, x2, y2, z2)) return 1;
  

 //Declare arrays for coordinates of galaxies on device
  double *dev_x1, *dev_y1, *dev_x2, *dev_y2, *dev_z1, *dev_z2, *dev_zS;
  double *dev_xS, *dev_yS, *dev_e1, *dev_e2, *dev_w;

    //Numbers of sources and lenses
  int N1, N2, NS;
  N1=x1.size(); //Number of galaxies
  N2=x2.size(); 
  NS=xS.size();
  
  // Allocate memory on device
  cudaError_t err1 = cudaMalloc(&dev_x1, N1*sizeof(double));
  cudaError_t err2 = cudaMalloc(&dev_y1, N1*sizeof(double));
  cudaError_t err3 = cudaMalloc(&dev_z1, N1*sizeof(double));

  if(err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess)
    {
      std::cerr<<"Could not allocate memory 1"<<std::endl;
      exit(1);
    };

  err1 = cudaMalloc(&dev_x2, N2*sizeof(double));
  err2 = cudaMalloc(&dev_y2, N2*sizeof(double));
  err3 = cudaMalloc(&dev_z2, N2*sizeof(double));

  if(err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess)
    {
      std::cerr<<"Could not allocate memory 2"<<std::endl;
      exit(1);
    };

  err1 = cudaMalloc(&dev_xS, NS*sizeof(double));
  err2 = cudaMalloc(&dev_yS, NS*sizeof(double));
  err3 = cudaMalloc(&dev_e1, NS*sizeof(double));
  cudaError_t err4 = cudaMalloc(&dev_e2, NS*sizeof(double));
  cudaError_t err5 = cudaMalloc(&dev_w, NS*sizeof(double));
  cudaError_t err6 = cudaMalloc(&dev_zS, NS*sizeof(double));

  std::cerr<<NS<<std::endl;

  if(err1 != cudaSuccess || err2 != cudaSuccess || err3 != cudaSuccess ||
     err4 != cudaSuccess || err5 != cudaSuccess || err6 != cudaSuccess)
    {
      std::cerr<<"Could not allocate memory 3"<<std::endl;
      exit(1);
    };

  // Copy values
  cudaMemcpy(dev_x1, x1.data(), N1*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_y1, y1.data(), N1*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_z1, z1.data(), N1*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_x2, x2.data(), N2*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_y2, y2.data(), N2*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_z2, z2.data(), N2*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_xS, xS.data(), NS*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_yS, yS.data(), NS*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_e1, e1.data(), NS*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_e2, e2.data(), NS*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_w, w.data(), NS*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_zS, zS.data(), NS*sizeof(double), cudaMemcpyHostToDevice);


  // Reading in Functions

  
  //Declare array for omega, sigma_crit, and D_A
  std::vector<double> omega, sigma_crit, angular_distance;
  double* dev_omega=NULL;
  double* dev_sigma_crit=NULL;
  double* dev_angular_distance=NULL;

  double omega_theta_min, omega_theta_max, sigma_crit_z_min, sigma_crit_z_max,
    angular_distance_z_min, angular_distance_z_max;
  
  // Try reading and stop if not successful
  if(g3lcong::readFunctionLog2Dev(filename_omega, num_bins, omega_theta_min,
		      omega_theta_max, omega)) return 1;
  if(g3lcong::readFunction2Dev(filename_sigma_crit, num_bins, sigma_crit_z_min,
		      sigma_crit_z_max, sigma_crit)) return 1;
  if(g3lcong::readFunction2Dev(filename_angular_distance, num_bins,
		      angular_distance_z_min, angular_distance_z_max,
		      angular_distance))
    return 1;

  err1 = cudaMalloc(&dev_omega, num_bins*sizeof(double));
  err2 = cudaMalloc(&dev_sigma_crit, num_bins*sizeof(double));
  err3 = cudaMalloc(&dev_angular_distance, num_bins*sizeof(double));

  cudaMemcpy(dev_omega, omega.data(), num_bins*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_sigma_crit, sigma_crit.data(), num_bins*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dev_angular_distance, angular_distance.data(), num_bins*sizeof(double), cudaMemcpyHostToDevice);



  // Calculate Gtilde

  //Declare container for Greal, Gimag and weight on device
  double *dev_Greal, *dev_Gimag, *dev_weight;

  //Allocate memory for Greal, Gimag and weight on device
  cudaMalloc(&dev_Greal, num_bins*num_bins*num_bins*sizeof(double));
  cudaMalloc(&dev_Gimag, num_bins*num_bins*num_bins*sizeof(double));
  cudaMalloc(&dev_weight, num_bins*num_bins*num_bins*sizeof(double));

  //Check if memory for Gtilde could be allocated on device
  if(0==dev_Greal || 0==dev_Gimag || 0==dev_weight)
    {
      std::cout<<"Couldnt allocate Space for Gtilde"<<std::endl;
      return 1;
    };

  //Set Gtilde to 0 on device
  cudaMemset(dev_Greal, 0, num_bins*num_bins*num_bins*sizeof(double));
  cudaMemset(dev_Gimag, 0, num_bins*num_bins*num_bins*sizeof(double));
  cudaMemset(dev_weight, 0, num_bins*num_bins*num_bins*sizeof(double));



  double sigma2=sigmaZ*sigmaZ; //Square of stddev of z weighting
  if(physical)
    {      
      return 1;
    }
  else
    {
      g3lcong::addToGtildeSourceRedshift<<<BLOCKS, THREADS>>>(dev_x1, dev_y1, dev_x2, dev_y2,
						dev_xS, dev_yS, dev_z1, dev_z2, dev_zS,
						dev_e1, dev_e2, dev_w,
						dev_omega,
						sigma2, omega_theta_min,
						omega_theta_max, num_bins, N1,
						N2, NS, theta_min, theta_max,
						dev_Greal, dev_Gimag,
						dev_weight);
    };


  
  //Declare arrays for total Greal, Gimag and weight on host
  double *Greal_tot, *Gimag_tot, *weight_tot;
  
  //Allocate memory for total Greal, Gimag and weight on host
  Greal_tot=(double *) malloc(num_bins*num_bins*num_bins*sizeof(double));
  Gimag_tot=(double *) malloc(num_bins*num_bins*num_bins*sizeof(double));
  weight_tot=(double *) malloc(num_bins*num_bins*num_bins*sizeof(double));

  if(Greal_tot==NULL || Gimag_tot==NULL || weight_tot==NULL)
    {
      std::cerr<<"calculateGtilde_gpu: Couldn't allocate memory"<<std::endl;
      exit(1);
    };

    //Copy Gtilde from device to host
  cudaMemcpy(Greal_tot, dev_Greal, num_bins*num_bins*num_bins*sizeof(double),
	     cudaMemcpyDeviceToHost);
  cudaMemcpy(Gimag_tot, dev_Gimag,  num_bins*num_bins*num_bins*sizeof(double),
	     cudaMemcpyDeviceToHost);
  cudaMemcpy(weight_tot, dev_weight, num_bins*num_bins*num_bins*sizeof(double),
	     cudaMemcpyDeviceToHost);



  
  // Output

    // Print out vartheta1, vartheta2, psi, binsize, Gtilde
  double phi_binsize=(phi_max - phi_min)/num_bins;
  double theta_binsize=log(theta_max/theta_min)/num_bins;
  
  for(int i=0; i<num_bins; i++)
    {
      //Theta1 Center of this bin
      double theta1=0.5*(exp(log(theta_min)+theta_binsize*(i+1))
			 + exp(log(theta_min)+theta_binsize*i));
      //Theta1 Binsize of this bin
      double deltaTheta1= exp(log(theta_min)+theta_binsize*(i+1))
	- exp(log(theta_min)+theta_binsize*i);

      for(int j=0; j<num_bins; j++)
	{
	  //Theta2 Center of this bin
	  double theta2=0.5*(exp(log(theta_min)+theta_binsize*(j+1))
			 + exp(log(theta_min)+theta_binsize*j));
	  //Theta2 Binsize of this bin
	  double deltaTheta2= exp(log(theta_min)+theta_binsize*(j+1))
	    - exp(log(theta_min)+theta_binsize*j);
	  
	  for(int k=0; k<num_bins; k++)
	    {
	      // Phi Center of this bin
	      double phi=(k+0.5)*phi_binsize + phi_min;

	      int index=i*num_bins*num_bins + j*num_bins + k;


	      // Weight
	      double weight=weight_tot[index];
	      //Greal
	      double Greal=Greal_tot[index];
	      //Gimag
	      double Gimag=Gimag_tot[index];

	      if(weight!=0)
		{
		  Greal/=weight;
		  Gimag/=weight;
		}
	      else
		{
		  Greal=0;
		  Gimag=0;
		}
	       
	      
	      // Output
	      std::cout
		<< theta1 <<" " //bin center theta 1 [arcmin or Mpc]
		<< theta2 <<" " //bin center theta 2 [arcmin or Mpc]
		<< phi <<" " //phi center [radians]
		<< deltaTheta1 <<" " //bin size theta 1[arcmin or Mpc]
		<< deltaTheta2 <<" " //bin size theta 2[arcmin or Mpc]
		<< phi_binsize <<" " //phi bin size [radians]
		<< Greal <<" "	//Real part of Gtilde [dimensionless or Msun/Mpc²]
		<< Gimag <<" "	//Imaginary part of Gtilde [dimensionless or Msun/Mpc²]
		<< weight // Weight of Gtilde [dimensionless or Msun/Mpc²]
		<<std::endl;
	    };
	};
    };

  


  
  return 0;
}
