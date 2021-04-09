#ifndef G3LCONG_APERTURESTATISTICS_H
#define G3LCONG_APERTURESTATISTICS_H

#include <complex>
#include <string>
#include <tuple>
#include <vector>

namespace g3lcong
{

  /**
   * Class for the calculation of the aperture statistics from Gtilde
   * Uses the exponential filter from Schneider & Watts (2005)
   * Depending on Gtilde computes either the dimensionless <N2M> or 
   * the <N2M>_phys in [Msun/Mpc²]
   *
   * @author Laila Linke llinke@astro.uni-bonn.de
   */
  class ApertureStatistics
  {
  private:

    /// Container for first bin axis of Gtilde [arcmin or Mpc]
    std::vector<double> vartheta1_;
    /// Container for second bin axis of Gtilde [arcmin or Mpc]
    std::vector<double> vartheta2_;
    /// Container for third bin axis of Gtilde [rad]
    std::vector<double> phi_;
    /// Container for bin volumes [arcmin^2*rad or Mpc^2*rad]
    std::vector<double> V_;    
    /// Container for Gtilde
    std::vector<std::complex<double>> Gtilde_;


    /**
     * Helping Parameters a and Theta8 
     * \f$ a=2\frac{\theta_1^2 \theta_2^2 \theta_3^2}{\theta_1^2\theta_2^2 + \theta_1^2\theta_3^2 + \theta_2^2\theta_3^2}\f$
     * \f$ \Theta^8 = \frac{(\theta_1^2\theta_2^2 + \theta_1^2\theta_3^2 + \theta_2^2\theta_3^2)^2}{9} \f$
     * @param theta1 first coordinate [arcmin or Mpc]
     * @param theta2 second coordinate [arcmin or Mpc]
     * @param theta3 third coordinate [arcmin or Mpc]
     * @return Tuple containing a and Theta8
     */
    std::tuple<double, double> a_Theta8(const double& theta1, const double& theta2, const double& theta3);

    
    /**
     * Integral Kernel for N2Map
     * @param vartheta1 Distance lens1-source [arcmin or Mpc]
     * @param vartheta2 Distance lens2-source [arcmin or Mpc]
     * @param psi Angle between lens1-source and lens2-source in [radians]
     * @param theta1 first aperture radius [arcmin or Mpc]
     * @param theta2 second aperture radius [arcmin or Mpc]
     * @param theta3 third aperture radius [arcmin or Mpc]
     * @param a2 a-parameter, precomputed for the thetas
     * @param bigT Theta8 parameter, precomputed for the thetas
     */
    std::complex<double> A_NNM(const double& vartheta1, const double& vartheta2, const double& phi,
			       const double& theta1, const double& theta2, const double& a2, const double& bigT);

    /**
     * Reads in Gtilde file to Gtilde_, set bins and computes bin volumes
     * @param filename Filename for Gtilde
     * @param tesselated True, if Gtilde has been tesselated, false if Gtilde is on regular grid
     */
    void readGtilde(std::string filename, bool tesselated);
    
  public:

    ///Empty Constructor
    ApertureStatistics(){};

    /**
     * Constructor from Coordinates
     * Automatically computes \f$a\f$ and \f$\Theta^8\f$
     * @param theta1 first coordinate [arcmin or Mpc]
     * @param theta2 second coordinate [arcmin or Mpc]
     * @param theta3 third coordinate [arcmin or Mpc]
     */
    //  ApertureStatistics(const double& theta1, const double& theta2,
    //		       const double& theta3);

    /**
     * Constructor from Gtilde filename
     * Reads in Gtilde, bins and computes bin volumes
     * @param filenameGtilde Filename for Gtilde
     * @param tesselated True, if Gtilde has been tesselated, false if Gtilde is on regular grid
     */
    ApertureStatistics(std::string filenameGtilde, bool tesselated);
    
    /**
     * Calculates N2Map for the given aperture radii
     * \f$ \langle N^2M_{\textrm{ap}} \rangle = \sum_i \sum_j \sum_k V_{\text{bin}} A_{NNM}(\vartheta_i, \vartheta_j, \psi_k) G(\vartheta_i, \vartheta_j, \psi_k) \f$
     * @param theta1 first aperture radius [arcmin or Mpc]
     * @param theta2 second aperture radius [arcmin or Mpc]
     * @param theta3 third aperture radius [arcmin or Mpc]
     * @retval NNMap and NNMperp in complex number [Units of Gtilde]
     */ 
    std::complex<double> NNM(const double& theta1, const double& theta2, const double& theta3);

 



  };
}
#endif //G3LCONG_APERTURESTATISTICS_H
