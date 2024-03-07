#include "helpers.h"


void g3lcong::checkCmdLine(int argc, int n_params, std::string usage, std::string example)
  {
    if(argc != n_params+1)
      {
	std::cerr<<"Wrong number of CMD Line arguments"<<std::endl;
  std::cerr<<"Expected:"<<n_params<<std::endl;
  std::cerr<<"Got:"<<argc-1<<std::endl;
	std::cerr<<"Usage:"<<usage<<std::endl;
	std::cerr<<"Example:"<<example<<std::endl;
	exit(1);
      };
  }


