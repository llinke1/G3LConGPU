#include "helpers.h"
#include "Function.h"
#include <vector>
#include <fstream>
#include <sstream>

// Number of blocks/Threads (adjust for GPU)
#define BLOCKS
#define THREADS



int g3lcong::readLenses2Dev(std::string filename,  const int& total_columns,
			    const int& col_x, const int& col_y,
			    const int& col_z, std::vector<double>& x,
			    std::vector<double>& y, std::vector<double>& z)
{


  //Reading in
  std::ifstream input(filename.c_str());

  if(input.fail())
    {
      std::cerr<<"readLenses2Dev: Couldn't open "<<filename<<std::endl;
      return 1;
    };
  
  
  std::string line;

  while(getline(input, line))
    {
      std::vector<double> values(total_columns);
      std::stringstream ss;
      ss.str(line);
      for(int i=0; i<total_columns; i++)
	{
	  ss>>values.at(i);
	  
	};
      x.push_back(values.at(col_x-1));
      y.push_back(values.at(col_y-1));
      z.push_back(values.at(col_z-1));

    };


  // Success
  return 0;
}



int g3lcong::readSources2Dev(std::string filename,  const int& total_columns,
			     const int& col_x, const int& col_y,
			     const int& col_e1, const int& col_e2,
			     const int& col_w, std::vector<double> &x,
			     std::vector<double> &y, std::vector<double> &e1,
			     std::vector<double> &e2, std::vector<double> &w,
           bool flipE1, bool flipE2)
{


  
  //Reading in
  std::ifstream input(filename.c_str());

  if(input.fail())
    {
      std::cerr<<"readSources2Dev: Couldn't open "<<filename<<std::endl;
      return 1;
    };
  
  std::string line;

  while(getline(input, line))
    {
      std::vector<double> values(total_columns);
      std::stringstream ss;
      ss.str(line);
      for(int i=0; i<total_columns; i++)
	{
	  ss>>values.at(i);
	  
	};
      x.push_back(values.at(col_x-1));
      y.push_back(values.at(col_y-1));
      e1.push_back(values.at(col_e1-1));
      e2.push_back(values.at(col_e2-1));
      w.push_back(values.at(col_w-1));

    };

    if (flipE1)
    {
      for (int i=0; i<e1.size(); i++ )
      {
        e1[i]*=(-1);
      }
    };

    if (flipE2)
    {
      for (int i=0; i<e2.size(); i++ )
      {
        e2[i]*=(-1);
      }
    };

  // Success
  return 0;
}


int g3lcong::readFunction2Dev(std::string filename,  const int& num_bins,
			      double& min,  double& max, std::vector<double> &values)
{
  // Read in Function
  Function func(filename, 0.);
  min=func.x_values_.at(0);
  max=func.x_values_.back();
  double bin=(max-min)/num_bins; //Width of a bin
  for(int i=0; i<num_bins; i++)
    {
      double x=min+i*bin;
      double val=func.at(x);
      values.push_back(val);
    };

  // Success
  return 0;
}


int g3lcong::readFunctionLog2Dev(std::string filename,  const int& num_bins,
			      double& min,  double& max, std::vector<double> &values)
{
  // Read in Function
  Function func(filename, 0.);
  min=func.x_values_.at(0);
  max=func.x_values_.back();
  double bin=log(max/min)/num_bins; //Logwidth of a bin
  for(int i=0; i<num_bins; i++)
    {
      double x=exp(log(min)+i*bin);
      double val=func.at(x);
      values.push_back(val);
    };

  // Success
  return 0;
}
