#ifndef G3LCONG_TWOPTSTAT_H
#define G3LCONG_TWOPTSTAT_H

#include <vector>

#include "Kdtree.h"
#include "Function.h"

#define THREADS 8 //Number of Threads for KdTree parallelization

namespace g3lcong
{
  /**
   * Class for the computation of the angular galaxy-correlation function
   *
   * @author Laila Linke llinke@astro.uni-bonn.de
   */
  class TwoPtStat
  {
  public:

    ///Binning for the distance between galaxies
    double theta_min_, theta_max_, theta_binsize_;

    ///Number of Bins
    int num_bins_;

    ///Container for the Paircount
    std::vector<double> counts_total_;

    ///Empty Constructor
    TwoPtStat(){};

    /**
     * Constructor from Binning
     *
     * @param theta_min minimal theta [arcmin]
     * @param theta_max maximal theta [arcmin]
     * @param num_bins number of bins
     */
    TwoPtStat(const double& theta_min, const double& theta_max,
	      const int& num_bins);

    /**
     * Calculate distance between two nodes
     *
     * @param pointer1 First Node
     * @param pointer2 Second Node
     * @param d_min will contain minimal distance between nodes [arcmin]
     * @param d_max will contain maximal distance between nodes [arcmin]
     * @retval distance between centers of nodes [arcmin]
     */
    double boxDistance(Kdnode* pointer1, Kdnode* pointer2, double& d_min,
		       double& d_max);

    ///Splits Kdnode List for parallel processing
    void split(std::vector<Kdnode*>& list);

    ///Brute Force Implementation for Dual Tree Count
    void dualTreeCount_bruteForce( Kdtree& tree1, Kdtree& tree2,
				   const double& sigmaZ);

    /**
     * Computes Paircount
     * Uses Parallelization with as many threads as processors
     *
     * @param pointer1 Starting Node of first Kdtree
     * @param pointer2 Starting Node of second Kdtree
     * @param sigmaZ Width of redshift weighting Gaussian
     */
    void dualTreeCount(Kdnode* pointer1, Kdnode* pointer2,
		       const double& sigmaZ);

    /**
     * Recursive subroutine for the computation of Paircount
     *
     * @param pointer1 Node with first galaxies
     * @param pointer2 Node with second galaxies
     * @param weight_redshift Function for redshift weighting
     * @param sigmaZ Width of redshift weighting Gaussian
     */
    void dualTreeCountSub(Kdnode* pointer1, Kdnode* pointer2,
			  const double& sigmaZ, 
			  std::vector<double>& counts);

  
  };

 
}



#endif //G3LCONG_TWOPTSTAT_H
