#ifndef KERNELFUNCTIONS_CUH
#define KERNELFUNCTIONS_CUH

namespace g3lcong
{

/**
 * Kernel function, which updates Gtilde in angular units
 * Function is executed for each thread and updates Gtilde and weight 
 * Gtilde and weight are stored in global (not shared) memory
 *
 * @param x1 x-coordinates of lens1 galaxies [arcmin], stored on device
 * @param y1 y-coordinates of lens1 galaxies [arcmin], stored on device
 * @param x2 x-coordinates of lens2 galaxies [arcmin], stored on device
 * @param y2 y-coordinates of lens2 galaxies [arcmin], stored on device
 * @param z1 redshift of lens1 galaxies, stored on device
 * @param z2 redshift of lens2 galaxies, stored on device
 * @param xS x-coordinates of source galaxies [arcmin], stored on device
 * @param yS y-coordinates of source galaxies [arcmin], stored on device
 * @param e1 ellipticity part1 of sources, stored on device
 * @param e2 ellipticity part2 of sources, stored on device
 * @param w weight of sources, stored on device
 * @param omega omega omega values, stored on device
 * @param sigma2 Square of stddev of z weighting
 * @param omega_theta_min Minimum theta for omega array [arcmin]
 * @param omega_theta_max Maximum theta for omega array [arcmin]
 * @param num_bins Number of bins for omega
 * @param N1 Number of lens 1 galaxies
 * @param N2 Number of lens 2 galaxies
 * @param NS Number of source galaxies
 * @param theta_min Minimal Theta [arcmin]
 * @param theta_max Maximal Theta [arcmin]
 * @param Greal, Storage for real part of Gtilde on global memory
 * @param Gimag, Storage for imag part of Gtilde on global memory
 * @param weight, Storage for weight on global memory
 */
__global__ void addToGtilde(double* x1, double* y1, double* x2, double* y2,
			    double* xS, double* yS, double* z1, double* z2,
			    double* e1, double* e2,
			    double *w, double *omega, double sigma2,
			    double omega_theta_min, double omega_theta_max,
			    int num_bins, int N1, int N2, int NS,
			    double theta_min, double theta_max,
			    double *Greal, double *Gimag,
			    double *weight);

__global__ void addToGtildeSourceRedshift(double* x1, double* y1, double* x2, double* y2,
				     double* xS, double* yS, double* z1, double* z2, double* zS,
				     double* e1, double* e2,
				     double *w, double *omega, double sigma2,
				     double omega_theta_min,
				     double omega_theta_max,
				     int num_bins, int N1, int N2, int NS,
				     double theta_min, double theta_max,
				     double *Greal, double *Gimag,
				     double *weight);

/**
 * Kernel function, which updates Gtilde in physical units w/o z weighting
 * Function is executed for each thread and updates Gtilde and weight 
 * Gtilde and weight are stored in global (not shared) memory
 *
 * @param x1 x-coordinates of lens1 galaxies [arcmin], stored on device
 * @param y1 y-coordinates of lens1 galaxies [arcmin], stored on device
 * @param x2 x-coordinates of lens2 galaxies [arcmin], stored on device
 * @param y2 y-coordinates of lens2 galaxies [arcmin], stored on device
 * @param z1 redshift of lens1 galaxies, stored on device
 * @param z2 redshift of lens2 galaxies, stored on device
 * @param xS x-coordinates of source galaxies [arcmin], stored on device
 * @param yS y-coordinates of source galaxies [arcmin], stored on device
 * @param e1 ellipticity part1 of sources, stored on device
 * @param e2 ellipticity part2 of sources, stored on device
 * @param w weight of sources, stored on device
 * @param omega omega omega values, stored on device
 * @param sigma2 Square of stddev of z weighting
 * @param sigma_crit Critical surface mass density (ave over sources) [Msun/Mpc²]
 * @param angular_distance Angular Diameter distance [Mpc]
 * @param omega_theta_min Minimum theta for omega array [arcmin]
 * @param omega_theta_max Maximum theta for omega array [arcmin]
 * @param sigma_crit_z_min Minimum z for Sigma Crit
 * @param sigma_crit_z_max Maximum z for Sigma Crit
 * @param angular_distance_z_min Minimum z for D_A
 * @param angular_distance_z_max Maximum z for D_A
 * @param num_bins Number of bins for omega
 * @param N1 Number of lens 1 galaxies
 * @param N2 Number of lens 2 galaxies
 * @param NS Number of source galaxies
 * @param r_min Minimal r [Mpc]
 * @param r_max Maximal r [Mpc]
 * @param Greal, Storage for real part of Gtilde on global memory [Msun/Mpc²]
 * @param Gimag, Storage for imag part of Gtilde on global memory [Msun/Mpc²]
 * @param weight, Storage for weight on global memory
 */
__global__ void addToGtildePhysical(double* x1, double* y1, double* x2, double* y2,
				    double* xS, double* yS, double* z1, double* z2,
				    double* e1, double* e2, double *w,
				    double *omega, double sigma2,
				    double *sigma_crit, double *angular_distance,
				    double omega_theta_min,
				    double omega_theta_max,
				    double sigma_crit_z_min,
				    double sigma_crit_z_max,
				    double angular_distance_z_min,
				    double angular_distance_z_max,
				    int num_bins, int N1, int N2, int NS,
				    double r_min, double r_max,
				    double *Greal, double *Gimag,
				    double *weight);




/**
 * Kernel function, which updates the paircount with or w/o z weighting
 * Function is executed for each thread and updates paircount 
 * Paircount is stored in global (not shared) memory
 *
 * @param x1 x-coordinates of lens1 galaxies [arcmin], stored on device
 * @param y1 y-coordinates of lens1 galaxies [arcmin], stored on device
 * @param x2 x-coordinates of lens2 galaxies [arcmin], stored on device
 * @param y2 y-coordinates of lens2 galaxies [arcmin], stored on device
 * @param z1 redshift of lens1 galaxies, stored on device
 * @param z2 redshift of lens2 galaxies, stored on device
 * @param N1 Number of lens 1 galaxies
 * @param N2 Number of lens 2 galaxies
 * @param theta_min Minimal Theta [arcmin]
 * @param theta_max Maximal Theta [arcmin]
 * @param num_bins Number of bins for omega
 * @param sigma2 Square of stddev of z weighting
 * @param paircount Storage for paircount on global memory
 */
__global__ void addToPaircount(double* x1, double* y1, double* x2, double* y2,
			       double* z1, double* z2, int N1, int N2,
			       int num_bins,
			       double theta_min, double binwidth, double sigma2,
			       double* paircount);


/**
 * Kernel function, which updates the triplecount 
 * Function is executed for each thread and updates triplecount 
 * Triplecount is stored in global (not shared) memory
 *
 * @param x1 x-coordinates of lens1 galaxies [arcmin], stored on device
 * @param y1 y-coordinates of lens1 galaxies [arcmin], stored on device
 * @param z1 redshift of lens1 galaxies, stored on device
 * @param x2 x-coordinates of lens2 galaxies [arcmin], stored on device
 * @param y2 y-coordinates of lens2 galaxies [arcmin], stored on device
 * @param z2 redshift of lens3 galaxies, stored on device
 * @param x3 x-coordinates of lens3 galaxies [arcmin], stored on device
 * @param y3 y-coordinates of lens3 galaxies [arcmin], stored on device
 * @param N1 Number of lens 1 galaxies
 * @param N2 Number of lens 2 galaxies
 * @param N3 Number of lens 3 galaxies
 * @param num_bins Number of bins for Triplecount
 * @param r_min minimal rp or pi [same unit as Dcom]
 * @param r_binwidth logarithmic binwidth of rp/pi
 * @param Dcom Precomputed comoving distance [Probably Mpc] (on device!)
 * @param z_min Minimal z for Dcom
 * @param z_binwidth binwidth of z for Dcom
 * @param triplecount Storage for triplecount on global memory
 */
  __global__ void addToTriplecount(double* x1, double* y1, double* z1, double* x2, double* y2, double* z2, double* x3, double* y3, double* z3, int N1, int N2, int N3, int num_bins, double r_min, double r_binwidth, double *Dcom, double z_min, double z_binwidth, int* triplecount); 

  
}
#endif //KERNELFUNCTIONS_CUH
