#ifndef GGGL_THREEPTSTAT_H
#define GGGL_THREEPTSTAT_H

#include <vector>
#include <complex>

#include "Kdtree.h"
#include "Function.h"

#define THREADS_KDTREE 8 //Number of threads for KdTree Parallelization

namespace g3lcong
{

  /**
   * Class for the computation of Gtilde 
   * Based on Code by Patrick Simon
   *
   * @author Laila Linke llinke@astro.uni-bonn.de
   * @date August 2018
   */
  class ThreePtStat
  {
  public:

    ///Binning for Distance Lens-Source [arcmin]
    double theta_min_, theta_max_;

    ///Logarithmic Binsize for theta [arcmin]
    double theta_binsize_;

    ///Binning for Physical Distance Lens-Source [Mpc]
    double r_min_, r_max_;

    ///Logarithmic Binsize for r [Mpc]
    double r_binsize_;

    
    ///Binning for Angle between Lens1-Source and Lens2-Source [rad]
    double phi_min_, phi_max_;

    ///Linear Binsize for Phi [rad]
    double  phi_binsize_;

    ///Number of Bins
    int num_bins_;

    ///Container for Real part of Gtilde
    std::vector<double> G_real_total_;

    ///Container for Imaginary part of Gtilde
    std::vector<double> G_imag_total_;

    ///Container for Weight
    std::vector<double> weight_total_;

    
    std::vector<double> theta1_com_total_;
    std::vector<double> theta2_com_total_;
    std::vector<double> phi_com_total_;



    
    
    ///Empty Constructor
    ThreePtStat(){};

    /**
     * Constructor from Binning
     *
     * @param theta_min minimal theta [arcmin]
     * @param theta_max maximal theta [arcmin]
     * @param phi_min minimal phi [radians]
     * @param phi_max maximal phi [radians]
     * @param num_bins number of bins
     */
    ThreePtStat(const double& theta_min, const double& theta_max,
		const double& phi_min, const double& phi_max,
		const int& num_bins);

    /**
     * Constructor from Binning for physical distances
     *
     * @param theta_min minimal theta [arcmin]
     * @param theta_max maximal theta [arcmin]
     * @param r_min minimal theta [Mpc]
     * @param r_max maximal theta [Mpc]
     * @param phi_min minimal phi [radians]
     * @param phi_max maximal phi [radians]
     * @param num_bins number of bins
     */
    ThreePtStat(const double& theta_min, const double& theta_max,
		const double& r_min, const double& r_max,
		const double& phi_min, const double& phi_max,
		const int& num_bins);
  

    
 

    /**
     * Computes Gtilde
     * Uses Parallelization with as many threads as processors
     *
     * @param sources Kdtree containing source galaxies
     * @param lenses1 Kdtree containing first lens galaxies
     * @param lenses2 Kdtree containing second lens galaxies
     * @param omega Function for angular two point correlation function
     * @param sigmaZ Stddev of Gaussian redshift weight (if 0: no weighting)
     */
    void tripleTreeCount(Kdtree sources, Kdtree lenses1, Kdtree lenses2,
			 Function* omega, const double& sigmaZ);


    /**
     * Computes GtildePhys
     * Uses Parallelization with as many threads as processors
     *
     * @param sources Kdtree containing source galaxies
     * @param lenses1 Kdtree containing first lens galaxies
     * @param lenses2 Kdtree containing second lens galaxies
     * @param omega Function for angular two point correlation function
     * @param sigmaZ Stddev of Gaussian redshift weight (if 0: no weighting)
     * @param sigma_crit Function with averaged critical surface mass density
     * @param angular_distance Function for angular diameter distance
     */
    void tripleTreeCount(Kdtree sources, Kdtree lenses1, Kdtree lenses2,
			 Function* omega, double sigmaZ, Function* sigma_crit,
			 Function* angular_distance);
    

    
    /**
     * Recursive subroutine for the computation of Gtilde
     *
     * @param source Kdnode with source galaxies
     * @param lens1 Kdnode with first lens galaxies
     * @param lens2 Kdnode with second lens galaxies
     * @param omega Function for angular two point correlation function
     * @param sigmaZ Stddev of Gaussian redshift weight (if 0: no weighting)
     * @param G_real Container for real part of Gtilde for this thread
     * @param G_imag Container for imaginary part of Gtilde for this thread
     * @param weight Container for weights in this thread
     * @param theta1_com Theta 1 of Center of Mass of Bin [arcmin]
     * @param theta2_com Theta 2 of Center of Mass of Bin [arcmin]
     * @param phi_com Phi of Center of Mass of Bin
     */
    void tripleTreeCountSub(Kdnode* source, Kdnode* lens1, Kdnode* lens2,
			    Function* omega, double sigmaZ,
			    std::vector<double>& G_real,
			    std::vector<double>& G_imag,
			    std::vector<double>& weight,
			    std::vector<double>& theta1_com,
			    std::vector<double>& theta2_com,
			    std::vector<double>& phi_com);


        /**
     * Recursive subroutine for the computation of GtildePhys
     *
     * @param source Kdnode with source galaxies
     * @param lens1 Kdnode with first lens galaxies
     * @param lens2 Kdnode with second lens galaxies
     * @param omega Function for angular two point correlation function
     * @param sigmaZ Stddev of Gaussian redshift weight (if 0: no weighting)
     * @param sigma_crit Function for critical surface density [M_sun/Mpc^2]
     * @param angular_distance Function for angular diameter distances [Mpc]
     * @param G_real_sigCrit Container for real part of GtildePhys in [M_sun/Mpc^2]
     * @param G_imag_sigCrit Container for imag part of GtildePhys in [M_sun/Mpc^2]
     * @param weight_sigCrit Container for weights multiplied with SigCrit
     * @param r1_com R 1 of Center of Mass of Bin [Mpc]
     * @param r2_com R 2 of Center of Mass of Bin [Mpc]
     * @param phi_com Phi of Center of Mass of Bin
     */
    void tripleTreeCountSub(Kdnode* source, Kdnode* lens1, Kdnode* lens2,
			    Function* omega, double sigmaZ,
			    Function* sigma_crit, Function* angular_distance,
			    std::vector<double>& G_real_sigCrit,
			    std::vector<double>& G_imag_sigCrit,
			    std::vector<double>& weight_sigCrit,
			    std::vector<double>& r1_com,
			    std::vector<double>& r2_com,
			    std::vector<double>& phi_com);


    

    /**
     * Calculates triangle distances within the system
     * 
     * @param source Kdnode with source galaxies
     * @param lens1 Kdnode with first lens galaxies
     * @param lens2 Kdnode with second lens galaxies
     * @param a distance lens1-source
     * @param b distance lens2-source
     * @param phi angle between lens1-source and lens2-source
     * @param c1 vector lens1-source
     * @param c2 vector lens2-source
     * @param accept is true if triangle fits to system
     */
    void defineTriangle(Kdnode* source, Kdnode* lens1, Kdnode* lens2,
			double& a, double& b, double& phi,
			std::complex<double>& c1, std::complex<double>& c2,
			bool& accept);


    

    
    ///Splits Kdnode List for parallel processing
    void split(std::vector<Kdnode*>& list);

    
  };
}
#endif //GGGL_THREEPTSTAT_H
