#ifndef G3LCONG_HELPERS_H
#define G3LCONG_HELPERS_H

#include <iostream>
#include <string>
#include <cstddef>
#include <vector>


namespace g3lcong
{
  /**
   * @brief Checks if number of CMD line parameters is correct
   * If number of CMD line parameters is not correct, process exits with error
   * @param argc Number of CMD line parameters
   * @param n_params Number of expected parameters (aside from exec name)
   * @param usage Usage description of exec
   * @param example Example usage description
   */
  void checkCmdLine(int argc, int n_params, std::string usage, std::string example);

  /**
   * @brief Structure for lens galaxies (i.e. no shear)
   */
  typedef struct
  {
    ///Theta1 [arcmin]
    double x_;
    ///Theta2 [arcmin]
    double y_;
    ///Redshift [dimensionless]
    double z_; 
  } lensGalaxy;

  /**
   * @brief Structure for source galaxies (i.e. with shear)
   */
  typedef struct
  {
    ///Theta1 [arcmin]
    double x_;
    ///Theta2 [arcmin]
    double y_;
    ///Redshift [dimensionless]
    double z_; 
    ///Ellipticity Coordinate 1 [dimensionless]
    double epsilon1_;
    ///Ellipticity Coordinate 2 [dimensionless]
    double epsilon2_;
    ///Weight of Galaxy [dimensionless]
    double weight_;
  } sourceGalaxy;
      

  /**
   * Reads in lens galaxies and writes them to device (including redshift)
   *
   * @param filename File with Galaxy Catalog
   * @param total_columns Number of Columns in Galaxy Catalog
   * @param col_x Column with X in arcmin
   * @param col_y Column with Y in arcmin
   * @param col_z Column with Redshift
   * @param x vector that should contain x
   * @param y vector that should contain y
   * @param z vector that should contain z
   * @param N Number of lenses read in
   * @return 0 if reading was successful
   */
  int readLenses2Dev(std::string filename,  const int& total_columns,
		     const int& col_x, const int& col_y, const int& col_z,
		     std::vector<double>& x, std::vector<double>& y,
		     std::vector<double>& z);




    /**
   * Reads in source galaxies and writes them to device
   *
   * @param filename File with Galaxy Catalog
   * @param total_columns Number of Columns in Galaxy Catalog
   * @param col_x Column with X in arcmin
   * @param col_y Column with Y in arcmin
   * @param col_epsilon1 Column with Epsilon1
   * @param col_epsilon2 Column with Epsilon2
   * @param col_weight Column with weight
   * @param x vector that should contain x
   * @param y vector that should contain y
   * @param e1 vector that should contain Epsilon 1
   * @param e2 vector that should contain Epsilon 2
   * @param w vector that should contain weight
   * @return 0 if reading was successful
   */
  int readSources2Dev(std::string filename,  const int& total_columns,
		      const int& col_x, const int& col_y, const int& col_e1,
		      const int& col_e2, const int& col_weight, std::vector<double> &x,
			     std::vector<double> &y, std::vector<double> &e1,
		      std::vector<double> &e2, std::vector<double> &w);

  
  /**
   * Reads in Function from file and writes to device array
   * Assumes linear binning
   *
   * @param filename File with Function
   * @param num_bins Number of Bins
   * @param min Will contain minimum value of x
   * @param max Will contain maximum value of x
   * @param values Vector that will contain function values
   * @return 0 if reading was successful
   */
  int readFunction2Dev(std::string filename,  const int& num_bins, double& min,
		       double& max, std::vector<double> &values);


    /**
   * Reads in Function from file and writes to device array
   * Assumes Log binning
   *
   * @param filename File with Function
   * @param num_bins Number of Bins
   * @param min Will contain minimum value of x
   * @param max Will contain maximum value of x
   * @param values Vector that will contain function values
   * @return 0 if reading was successful
   */
  int readFunctionLog2Dev(std::string filename,  const int& num_bins, double& min,
		       double& max, std::vector<double> &values);
};


#endif // G3LCONG_HELPERS_H
