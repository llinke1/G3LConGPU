#ifndef G3LCONG_FUNCTION_H
#define G3LCONG_FUNCTION_H

#include <vector>
#include <string>
#include <iostream>

namespace g3lcong
{
  /**
   * Class for tabularized functions
   *
   * @author Laila Linke llinke@astro.uni-bonn.de
   * @date August 2018
   */
  
  class Function
  {
  private:

    /**
     * @brief Returns the y value between two indices by linear interpolation
     * @param index1 First index
     * @param index2 Second index
     * @param x x-value corresponding to y
     */
    float interpolate(const int& index1, const int& index2, const float& x);

    ///Default Value of the Function
    float default_value_;
    
  public:

    ///X-Values
    std::vector<float> x_values_;

    ///Y-Values
    std::vector<float> y_values_;

    ///Empty Constructor
    Function(){};

    /**
     * @brief Constructor from file
     * @detail Reads in an ascii file and stores first column as x and second 
     * column as y. If filename="none", the default value is used for all x
     *
     * @warning Requires a two-column ascii file
     * @param filename File with tabularized function or "none"
     * @param default_value Default Value of the function
     */
    Function(std::string filename, const float& default_value);

    /**
     * @brief Returns y-value corresponding to x
     * @param x x-value
     */
    float at(const float& x);

    /**
     * @brief Returns the lower border of the x-bin corresponding to y
     * @param y y-Value
     */
    float x_lower(const float& y);
    
  };

  
}






#endif //G3LCONG_FUNCTION_H
