#include "TwoPtStat.h"

#include <cmath>

#include <iostream>

#include <omp.h>

g3lcong::TwoPtStat:: TwoPtStat(const double& theta_min, const double& theta_max,
			    const int& num_bins)
{

  theta_min_=theta_min;
  theta_max_=theta_max;
  num_bins_=num_bins;
  theta_binsize_=(log(theta_max)-log(theta_min))/num_bins;

  for(int i=0; i<num_bins; i++)
    {
      counts_total_.push_back(0);
    };
}

double g3lcong::TwoPtStat::boxDistance(Kdnode* pointer1, Kdnode* pointer2,
				    double& d_min, double& d_max)
{
  double dx=pointer1->center_x_ - pointer2->center_x_;
  double dy=pointer1->center_y_ - pointer2->center_y_;

  double r = sqrt(dx*dx + dy*dy);

  d_max = r + 0.5*(pointer1->diameter_ + pointer2->diameter_);
  d_min = r - 0.5*(pointer1->diameter_ + pointer2->diameter_);

  d_min = (d_min < 0 ? 0:d_min);

  return r;
    
}


void g3lcong::TwoPtStat::split(std::vector<Kdnode*>& list)
{
  for(unsigned int i=0; i<list.size(); i+=2)
    {
      Kdnode* tmp = list.at(i);
      list.at(i) = tmp->child_left_;
      if (i < list.size()-1)
	list.insert(list.begin()+i+1, 1, tmp->child_right_);
      else
	list.push_back(tmp->child_right_);
    };
}

void g3lcong::TwoPtStat::dualTreeCount_bruteForce(Kdtree& tree1, Kdtree& tree2, const double& sigmaZ)
{

 
  for(unsigned int i=0; i<tree1.getSize(); i++)
    {
      for(unsigned int j=0; j<tree2.getSize(); j++)
	{
	  double dx=tree1.galaxies_.at(i).x_-tree2.galaxies_.at(j).x_;
	  double dy=tree1.galaxies_.at(i).y_-tree2.galaxies_.at(j).y_;
	  double dz=tree1.galaxies_.at(i).z_-tree2.galaxies_.at(j).z_;
	  double theta=sqrt(dx*dx+dy*dy);

	  int index = floor((log(theta)-log(theta_min_))/(log(theta_max_)-log(theta_min_))*num_bins_);
	  
	  if (index>=0 && index < num_bins_)
	    {
	      double weight_redshift=exp(-dz*dz/2/sigmaZ/sigmaZ);
	      counts_total_.at(index)+=weight_redshift;
	    };
	};
    };
}




void g3lcong::TwoPtStat::dualTreeCount(Kdnode* pointer1, Kdnode* pointer2,
				       const double& sigmaZ)
{
  // Set number of threads to number of processors
  const int threads = THREADS; ///< number of threads
  omp_set_num_threads(threads);

  
  std::vector< std::vector<double>> count_list; ///<List of counts for threads

  // Initialization of counts for threads
  for(int i=0; i<threads; i++)
    {
      count_list.push_back(counts_total_);
    };

  // Starting Nodes
  std::vector<Kdnode*> node_list1;
  std::vector<Kdnode*> node_list2;

  node_list1.push_back(pointer1);
  node_list2.push_back(pointer2);


  // Jobmodels for the Threads
  switch(threads)
    {
    case 1: break;
    case 2:
      split(node_list1);
      break;
    case 4:
      split(node_list1);
      split(node_list2);
      break;
    case 8:
      split(node_list1);
      split(node_list1);
      split(node_list2);
      break;
    case 16:
      split(node_list1);
      split(node_list1);
      split(node_list2);
      split(node_list2);
      break;
    case 32:
      split(node_list1);
      split(node_list1);
      split(node_list1);
      split(node_list2);
      split(node_list2);
    default:
      std::cerr<<"TwoPtStat::dualTreeCount Wrong number of threads"<<std::endl;
      exit(1);
    }

  int jobindices[2][16]; ///< Starting Points for the jobs

   for(unsigned int i=0,count=0;i<node_list1.size();i++)
     {
       for(unsigned int j=0;j<node_list2.size();j++)
	 {
	   jobindices[0][count]   = i;
	   jobindices[1][count++] = j;
	 }
     };

   int tid; ///< Thread ID

#pragma omp parallel private(tid)
   {
     // Assign Thread ID
     tid = omp_get_thread_num();
     
     // Start Thread
     dualTreeCountSub(node_list1.at(jobindices[0][tid]),
		      node_list2.at(jobindices[1][tid]), sigmaZ,
		      count_list.at(tid));
   }

   // Wait for all jobs to finish
#pragma omp barrier
  
   // Combine job results
   for(unsigned int i=0;i<counts_total_.size();i++)
     {
       for(int j=0;j<threads;j++)
	 {
	   counts_total_.at(i)+=count_list.at(j).at(i);
	 }
     }
  
   return;   
}


void g3lcong::TwoPtStat::dualTreeCountSub(Kdnode* pointer1, Kdnode* pointer2,
					  const double& sigmaZ,
					  std::vector<double>& counts)
{
  // Stop if one of the pointers is a leaf
  if (pointer1==NULL || pointer2==NULL) return;

  // Calculate Maximum and Minimum Distances Between Node Boxes
  double dmax, dmin;
  double theta = boxDistance(pointer1, pointer2, dmin, dmax);

  // Stop if the pointers are too close or too far apart
  if(dmin)
    {
      if((log(dmin)>=log(theta_max_) || log(dmax)<log(theta_min_)))
	{
	  return;
	}
    };

  // Short Cut for Pointers which are the same
  if(pointer1 == pointer2)
    {
      dualTreeCountSub(pointer1->child_left_, pointer2->child_left_,
		       sigmaZ, counts);
      dualTreeCountSub(pointer1->child_right_, pointer2->child_right_,
		       sigmaZ, counts);
      dualTreeCountSub(pointer1->child_left_, pointer2->child_right_,
		       sigmaZ, counts);
      return;
    }
  

  if(sigmaZ!= 0) //Check if Redshift Weight has to be considered
    {
      // Split pointer 1 until only one galaxy is left
      if(pointer1->N_ > 1)
	{
	  dualTreeCountSub(pointer1->child_left_, pointer2, sigmaZ,
			   counts);
	  dualTreeCountSub(pointer1->child_right_, pointer2, sigmaZ,
			   counts);
	  return;
	};
      if(pointer2->N_ > 1) //split pointer 2
	{
	  dualTreeCountSub(pointer1, pointer2->child_left_, sigmaZ,
			   counts);
	  dualTreeCountSub(pointer1, pointer2->child_right_, sigmaZ,
		       counts);
	  return;
	}
    };

  
  // Add to Statistics, if the nodes are small enough
  if((dmin||theta==0)&&(pointer1->diameter_ + pointer2->diameter_)<=g3lcong::RECTOLERANCE*dmin)
    {
      // Calculate Index in Logbins
      int index = floor((log(theta)-log(theta_min_))/
			(log(theta_max_)-log(theta_min_))*num_bins_);

      // Calculate Redshift distance between galaxies
      double dz=pointer1->central_gal_.z_ - pointer2->central_gal_.z_;

      
      if(index>=0 && index<num_bins_)
	{
	  double weight_redshift=exp(-dz*dz/2/sigmaZ/sigmaZ);
	  counts.at(index)+=weight_redshift*pointer1->central_gal_.weight_*pointer2->central_gal_.weight_;
	  return;
	};
      return;
    };

  //Split if no redshift weight
  if(sigmaZ == 0) //Check if Redshift Weight has to be considered
    {
      // Split pointer 1 if it is larger than pointer 2
      if(pointer1->diameter_ > pointer2->diameter_)
	{
	  dualTreeCountSub(pointer1->child_left_, pointer2, sigmaZ,
			   counts);
	  dualTreeCountSub(pointer1->child_right_, pointer2, sigmaZ,
			   counts);
	  return;
	}
      else
	{
	  dualTreeCountSub(pointer1, pointer2->child_left_, sigmaZ,
			   counts);
	  dualTreeCountSub(pointer1, pointer2->child_right_, sigmaZ,
			   counts);
	  return;
	};
    };
  
  return;
  
}



