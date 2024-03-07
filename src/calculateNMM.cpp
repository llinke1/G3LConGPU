#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <complex>

#include "ApertureStatistics.h"
#include "helpers.h"

/**
 * Program computing <N M_ap^2> for predefined thetas
 *
 * Output is in units of Gplus / Gminus (unitless or Msun/Mpc^2)
 * Usage: calculateApertureStatistics.x filename_GplusGminus filename_Thetas
 *
 * Example: calculateApertureStatistics.x
 * ../products/gtilde_allTiles/all.gtilde.dat ../products/NNMap/thetas.dat 0
 *
 * @author Laila Linke, llinke@astro.uni-bonn.de
 */

int main(int argc, char *argv[])
{

    int n_params = 3;                                                                         // Expected number of params
    std::string usage = "calculateNMM.x filename_Gplus / Gminus filename_Thetas Tesselated?"; // Usage description

    std::string example = "calculateNNM.x gplusgminus_allTiles/all.gplusgminus.dat thetas.dat 0"; // Example usage

    // Check Number of CMD Line arguments
    g3lcong::checkCmdLine(argc, n_params, usage, example);

    // Reading in from commandline
    std::string filename_GplusGminus = argv[1]; // Name of File containing Gtilde
    std::string filename_Thetas = argv[2];      // Name of File with thetas for which NNMap shall be calculated
    bool tesselation = std::stoi(argv[3]);      // Bool if Gtilde was tesselated or not (0 if not tesselated, 1 else)

    // Initializing Aperture Statistics and read in Gtilde
    g3lcong::ApertureStatistics aperture_statistics(filename_GplusGminus, tesselation, "gplusgminus");

    // Reading in Thetas for which NMap² should be calculated and do calculation
    std::ifstream input(filename_Thetas);
    double theta1, theta2, theta3;
    while (input >> theta1 >> theta2 >> theta3)
    {
        // Calculate Aperture Statistics
        std::complex<double> NMM = aperture_statistics.NMM(theta1, theta2, theta3);
        std::complex<double> NMMstar = aperture_statistics.NMMstar(theta1, theta2, theta3);

        double NMapMap = aperture_statistics.NMapMap(NMM, NMMstar);
        double NMperpMperp = aperture_statistics.NMperpMperp(NMM, NMMstar);
        double NMapMperp = aperture_statistics.NMapMperp(NMM, NMMstar);

        // Print out Aperture Statistics
        std::cout << theta1 << " " << theta2 << " " << theta3 << " "
                  << NMapMap << " " << NMperpMperp << " " << NMapMperp <<" "
                  << real(NMM) << " " << imag(NMM) << " "<< real(NMMstar) << " " <<imag(NMMstar) <<" "
                  << std::endl;
    };

    return 0;
}
