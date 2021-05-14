#include "kernelfunctions.cuh"
#include "constants.h"

__global__ void g3lcong::addToGtilde(double* x1, double* y1, double* x2, double* y2,
				     double* xS, double* yS, double* z1, double* z2,
				     double* e1, double* e2,
				     double *w, double *omega, double sigma2,
				     double omega_theta_min,
				     double omega_theta_max,
				     int num_bins, int N1, int N2, int NS,
				     double theta_min, double theta_max,
				     double *Greal, double *Gimag,
				     double *weight)
{
  // global ID of this thread
  int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
  
  // Binwidth of omega
  double om_binwidth=log(omega_theta_max/omega_theta_min)/num_bins;
  double theta_binwidth=log(theta_max/theta_min)/num_bins;
  for(int k=thread_index; k<NS; k+=blockDim.x*gridDim.x) 
    {
      //Get source positions in arcmin
      double x_galS=xS[k];
      double y_galS=yS[k];
      double eps1=e1[k];
      double eps2=e2[k];
      double w_galS=w[k];
      
      //Go through all lens1 galaxies in submatrix
      for(int i=0; i<N1; i++)
	{

	  //Get positions of lens 1 in arcmin
	  double x_gal1=x1[i];
	  double y_gal1=y1[i];
	  
	  double dx1=x_galS-x_gal1;//x-distance betw lens1 & source [arcmin]
	  double dy1=y_galS-y_gal1;//y-distance betw lens1 & source [arcmin]
	  
	  double theta1=sqrt(dx1*dx1+dy1*dy1); //theta 1 [arcmin]

	  // Go through all lens2 galaxies in submatrix
	  for(int j=0; j<N2; j++)
	    {
	      //Get positions of lens 2 in arcmin
	      double x_gal2=x2[j];
	      double y_gal2=y2[j];
		  
	      double dx2=x_galS-x_gal2; //x-distance betw lens2 & source [arcmin]
	      double dy2=y_galS-y_gal2; //y-distance betw lens2 & source [arcmin]

	      double theta2=sqrt(dx2*dx2+dy2*dy2); //theta 2 [arcmin]
	      
	      //Get phi
	      double phi;

	      //dot = \vec{theta1}\dot\vec{theta2} = theta1*theta2*cos(phi)
	      double dot=dx1*dx2+dy1*dy2;
	      //det = \vec{theta1}\cross\vec{theta2} = theta1*theta2*sin(phi)
	      double det=dx1*dy2-dx2*dy1;
     
	      //atan2(x,y) returns angle alpha for which x=sin(alpha) and
	      //y=cos(alpha)
	      phi=atan2(det, dot);

	      //atan2 returns angle in [-pi, pi], we want [0,2pi]
	      if(phi<0) phi+=2*g3lcong::pi;

	       //Get index for Gtilde in logarithmic binning
	      int indexX=0;

	      if(theta1>theta_min) indexX=floor(log(theta1/theta_min)/theta_binwidth);
	     	      
	      int indexY=0;
	      if(theta2>theta_min) indexY=floor(log(theta2/theta_min)/theta_binwidth);     
	    
	      unsigned int index= indexX*num_bins*num_bins + indexY*num_bins
		+ floor(0.5*phi*num_bins/g3lcong::pi);
	
	      
	      if(indexX<num_bins && indexY<num_bins && index<num_bins*num_bins*num_bins)
		{
		  //Get Omega
		  double dx=x_gal1-x_gal2; //x-distance between lenses [arcmin]
		  double dy=y_gal1-y_gal2; //y-distance between lenses [arcmin]
		  
		  double r2=(dx*dx+dy*dy); //Squared distance between lenses [arcmin^2]
		  
		  int index_omega=0; //Index for omega
	
		  if(r2>omega_theta_min*omega_theta_min) index_omega=int((0.5*log(r2/omega_theta_min/omega_theta_min))/om_binwidth);
		  
		  double omega_triplet=0; //omega for this lens pair
		  if(index_omega<num_bins) omega_triplet=omega[index_omega];

		  //Get redshift weighting
		  double weightZ=1;
		  if(sigma2>0)
		    {
		      double dz=z1[i]-z2[j]; //redshift distance between lenses
		      weightZ = exp(-0.5*dz*dz/sigma2);
		    };
	      
		  //Get Phase Angle (Phi_i+Phi_j)
	      
		  double phi1=atan2(dy1, dx1); //phase angle of theta1
		  double phi2=atan2(dy2, dx2); //phase angle of theta2
	      
		  double cos_phase, sin_phase;
		  sincos(phi1+phi2, &sin_phase, &cos_phase);
	      
		  //Get Contribution of Triplet (multiplied by 0.01 to avoid overflow)
		  double Greal_triplet=(1+omega_triplet)
		    *(-eps1*cos_phase-eps2*sin_phase)*w_galS*weightZ*0.01; 
		  double Gimag_triplet=(1+omega_triplet)
		    *(eps1*sin_phase-eps2*cos_phase)*w_galS*weightZ*0.01;
		  

		  // Add Triplet contribution
		  atomicAdd(&Greal[index], Greal_triplet);
		  atomicAdd(&Gimag[index], Gimag_triplet);
		  atomicAdd(&weight[index],w_galS*0.01*weightZ);
		}
	    }
	}
    }
  return;
}


__global__ void g3lcong::addToGtildePhysical(double* x1, double* y1, double* x2, double* y2,
					     double* xS, double* yS, double* z1, double* z2,
					     double* e1, double* e2, double *w,
					     double *omega, double sigma2,
					     double *SigCrit, double *D_A,
					     double omega_theta_min,
					     double omega_theta_max,
					     double sigma_crit_z_min,
					     double sigma_crit_z_max,
					     double angular_distance_z_min,
					     double angular_distance_z_max,
					     int num_bins, int N1, int N2, int NS,
					     double r_min, double r_max,
					     double *Greal, double *Gimag,
					     double *weight)
{
  // global ID of this thread
  int thread_index = blockIdx.x * blockDim.x + threadIdx.x;

  // binsize for D_A and SigCrit
  double sigcrit_binwidth=(sigma_crit_z_max - sigma_crit_z_min)/num_bins;
  double angular_distance_binwidth=(angular_distance_z_max - angular_distance_z_min)/num_bins;
  // Binwidth of omega
  double om_binwidth=log(omega_theta_max/omega_theta_min)/num_bins;
  double r_binwidth=log(r_max/r_min)/num_bins;
  
  for(int k=thread_index; k<NS; k+=blockDim.x*gridDim.x)  // check that we do not go outside of submatrix border
    {
      //Get source positions in arcmin
      double x_galS=xS[k];
      double y_galS=yS[k];
      double eps1=e1[k];
      double eps2=e2[k];
      double w_galS=w[k];
      
      //Go through all lens1 galaxies in submatrix
      for(int i=0; i<N1; i++)
	{
	  //Get positions of lens 1 in arcmin
	  double x_gal1=x1[i];
	  double y_gal1=y1[i];
	  
	  double dx1=x_galS-x_gal1;//x-distance betw lens1 & source [arcmin]
	  double dy1=y_galS-y_gal1;//y-distance betw lens1 & source [arcmin]

	  // Go through all lens2 galaxies in submatrix
	  for(int j=0; j<N2; j++)
	    {
	      //Get positions of lens 2 in arcmin
	      double x_gal2=x2[j];
	      double y_gal2=y2[j];
		  
	      double dx2=x_galS-x_gal2; //x-distance betw lens2 & source [arcmin]
	      double dy2=y_galS-y_gal2; //y-distance betw lens2 & source [arcmin]

	      //Get distance and SigCrit	  
	      double z_ave=0.5*(z1[i]+z2[j]); //average z of lenses
	      
	      int sigma_crit_z_ind = (int) floor((z_ave-sigma_crit_z_min)/sigcrit_binwidth); //Index of z in SigCrit
	      int angular_distance_z_ind = (int) floor((z_ave-angular_distance_z_min)/angular_distance_binwidth); //Index of z in D_A
		  
	      double D_A_triplet=0;
	      double SigCrit_triplet=0;

	      //Do linear interpolation
	      D_A_triplet =D_A[angular_distance_z_ind]+(D_A[angular_distance_z_ind+1]-D_A[angular_distance_z_ind])/angular_distance_binwidth*(z_ave-angular_distance_z_min-angular_distance_z_ind*angular_distance_binwidth); //D_A of this triplet [Mpc]
	      SigCrit_triplet=SigCrit[sigma_crit_z_ind]+(SigCrit[sigma_crit_z_ind+1]-SigCrit[sigma_crit_z_ind])/sigcrit_binwidth*(z_ave-sigma_crit_z_min-sigma_crit_z_ind*sigcrit_binwidth); //SigCrit of this triplet [Msun/Mpc^2]
					  
	      //Get r1
	      double r1=D_A_triplet*sqrt(dx1*dx1+dy1*dy1)/g3lcong::rad_in_arcmin; //r1 [pc]
	      
	      //Get r2    
	      double r2=D_A_triplet*sqrt(dx2*dx2+dy2*dy2)/g3lcong::rad_in_arcmin; //r2[ pc]
	      
	      //Get phi
		 	      
	      //dot = \vec{theta1}\dot\vec{theta2} = theta1*theta2*cos(phi)
	      double dot=dx1*dx2+dy1*dy2;
	      //det = \vec{theta1}\cross\vec{theta2} = theta1*theta2*sin(phi)
	      double det=dx1*dy2-dx2*dy1;
	      
	      //atan2(x,y) returns angle alpha for which x=a*sin(alpha) and y=a*cos(alpha)
	      double phi=atan2(det, dot);
		      
	      //atan2 returns angle in [-pi, pi], we want [0,2pi]
	      if(phi<0) phi+=2*g3lcong::pi;

	      //Get index for Gtilde in logarithmic binning
	      int indexX=0;
	      if(r1>r_min) indexX=floor(log(r1/r_min)/r_binwidth);
		      
	      int indexY=0;
	      if(r2>r_min) indexY=floor(log(r2/r_min)/r_binwidth);
		      
	      
	      unsigned int index= indexX*num_bins*num_bins + indexY*num_bins + floor(0.5*phi*num_bins/g3lcong::pi);
		      

	      if(index < num_bins*num_bins*num_bins && indexX < num_bins && indexY < num_bins && phi < 2*g3lcong::pi
		 && SigCrit_triplet > 0)
		{
		  //Get Omega
		  double dx=x_gal1-x_gal2; //x-distance between lenses [arcmin]
		  double dy=y_gal1-y_gal2; //y-distance between lenses [arcmin]
			  
		  double r2=(dx*dx+dy*dy); //Squared distance between lenses [arcmin^2]
			  
		  int index_omega=0; //Index for omega

		  if(r2>0) index_omega=int((0.5*log(r2/omega_theta_min/omega_theta_min))/om_binwidth);
			  
		  double omega_triplet=0; //omega for this lens pair
		  if(index_omega<num_bins) omega_triplet=omega[index_omega];

		  //Get redshift weight
		  double weightZ=1;
		  if(sigma2!=0)
		    {
		      double dz=z1[i]-z2[j];
		      weightZ=exp(-0.5*dz*dz/sigma2);
		    };
			  
		  //Get Phase Angle (Phi_i+Phi_j)     
		  double phi1=atan2(dy1, dx1); //phase angle of theta1
		  double phi2=atan2(dy2, dx2); //phase angle of theta2
			  
		  double cos_phase, sin_phase;
		  sincos(phi1+phi2, &sin_phase, &cos_phase);

		  //Get Contribution of Triplet
		  double Greal_triplet=(1+omega_triplet)
		    *(-eps1*cos_phase-eps2*sin_phase)*w_galS*weightZ/SigCrit_triplet; 
		  double Gimag_triplet=(1+omega_triplet)
		    *(eps1*sin_phase-eps2*cos_phase)*w_galS*weightZ/SigCrit_triplet;     
		  
		  //Add Triplet contribution	      
		  atomicAdd(&Greal[index], Greal_triplet);
		  atomicAdd(&Gimag[index], Gimag_triplet);
		  atomicAdd(&weight[index],w_galS*weightZ/SigCrit_triplet/SigCrit_triplet);
		};
	    };
	};
    };
  return;
}

__global__ void g3lcong::addToPaircount(double* x1, double* y1,
					double* x2, double* y2,
					double* z1, double* z2,
					int N1, int N2, int num_bins,
					double theta_min, double binwidth,
					double sigma2, double* paircount)
{
 
  //Global ID of this thread
  int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
  
  for(int i=thread_index; i<N1; i+=blockDim.x*gridDim.x) 
    {
      double x_gal1=x1[i];
      double y_gal1=y1[i];
      double z_gal1=z1[i];
      
      //Go through all lens2 galaxies in submatrix
      for(int j=0; j<N2; j++)
	{
	  //X-difference between galaxy 1 and galaxy 2 [arcmin]
	  double dx=x_gal1-x2[j];
	  //Y-difference between galaxy 1 and galaxy 2 [arcmin]
	  double dy=y_gal1-y2[j];
	  //Squared distance between galaxy 1 and galaxy 2 [arcmin^2]
	  double r2=(dx*dx+dy*dy);

	  //Get index of this galaxy pair in paircount (logarithmic binning!)
	  int index=-1; //Set this to zero so the same galaxy is not counted twice!
	  if(r2 > 0) index=int((0.5*log(r2)-log(theta_min))/binwidth);

	  //Update paircount in shared memory
	  if(index>=0 && index<num_bins)
	    {
	      //Get redshift weight
	      double weightZ=1;
	      if(sigma2!=0)
		{
		  double dz=z_gal1-z2[j];
		  weightZ=exp(-0.5*dz*dz/sigma2);
		};
	      
	      atomicAdd(&paircount[index],weightZ*0.01); //Add only 0.01*w per pair to avoid overflow!   
	    };
	}     
    }

  
  return;
}


__global__ void g3lcong::addToTriplecount(double* x1, double* y1, double* z1, double* x2, double* y2, double* z2, double* x3, double* y3, double* z3, int N1, int N2, int N3, int num_bins, double r_min, double r_binwidth, double *Dcom, double z_min, double z_binwidth, int* triplecount)
{
  // global ID of this thread
  int thread_index = blockIdx.x * blockDim.x + threadIdx.x;
  
  for(int i=thread_index; i<N1; i+=blockDim.x*gridDim.x)
    {
      // get physical coordinates of galaxies
      double zgal1=z1[i]; //redshift
      int index_z = int(zgal1-z_min/z_binwidth); //index of redshift
      double d1=Dcom[index_z]; //Comoving distance to galaxy [Mpc]

      double xgal1=d1*x1[i]; //Physical x coordinate [Mpc]
      double ygal1=d1*y1[i]; //Physical y coordinate [Mpc]

      for(int j=0; j<N2; j+=blockDim.x*gridDim.x)
	{
	  double zgal2=z2[i]; //redshift
	  index_z = int(zgal2-z_min/z_binwidth); //index of redshift
	  double d2=Dcom[index_z]; //Comoving distance to galaxy [Mpc]

	  double xgal2=d2*x2[j]; //Physical x coordinate [Mpc]
	  double ygal2=d2*y2[j]; //Physical y coordinate [Mpc]

	  // Projected Separation between gal 1 und gal 2 [Mpc]
	  double rp12=sqrt((xgal2-xgal1)*(xgal2-xgal1)+(ygal2-ygal1)*(ygal2-ygal1));
	  // Line-of-Sight distance between gal1 und gal 2 [Mpc]
	  double pi12=abs(d2-d1);

	  int index_rp12=int(log(rp12/r_min)/r_binwidth);
	  int index_pi12=int(log(pi12/r_min)/r_binwidth);

	  if(index_rp12 < num_bins && index_pi12 < num_bins)
	    {

	      for(int k=0; k<N3; k+=blockDim.x*gridDim.x)
		{
		  double zgal3=z3[i]; //redshift
		  index_z = int(zgal3-z_min/z_binwidth); //index of redshift
		  double d3=Dcom[index_z]; //Comoving distance to galaxy [Mpc]
		  
		  double xgal3=d3*x3[j]; //Physical x coordinate [Mpc]
		  double ygal3=d3*y3[j]; //Physical y coordinate [Mpc]
		  
		  // Projected Separation between gal 1 und gal 3 [Mpc]
		  double rp13=sqrt((xgal3-xgal1)*(xgal3-xgal1)+(ygal3-ygal1)*(ygal3-ygal1));
		  // Projected Separation between gal 2 und gal 3 [Mpc]
		  double rp23=sqrt((xgal3-xgal2)*(xgal3-xgal2)+(ygal3-ygal2)*(ygal3-ygal2));
		  
		  // Line-of-Sight distance between gal2 und gal 3 [Mpc]
		  double pi23=abs(d3-d2);
		  
		  int index_rp13=int(log(rp13/r_min)/r_binwidth);
		  int index_rp23=int(log(rp23/r_min)/r_binwidth);
		  int index_pi23=int(log(pi23/r_min)/r_binwidth);
		  
		  if(index_rp13 < num_bins && index_rp23<num_bins && index_pi23 > num_bins)
		    {
		  
		      int index=index_rp12*num_bins*num_bins*num_bins*num_bins
			+index_rp13*num_bins*num_bins*num_bins
			+index_rp23*num_bins*num_bins
			+index_pi12*num_bins
			+index_pi23;
		  
		      atomicAdd(&triplecount[index], 1);
		    };
		};
	    };
	}
    }
}
