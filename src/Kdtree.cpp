#include "Kdtree.h"

#include <cmath>

#include <fstream>
#include <vector>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "helpers.h"

g3lcong::Kdtree:: Kdtree(const std::string& filename, const int& total_columns, const int& col_x, const int& col_y, const int& col_z, const int& col_epsilon1, const int& col_epsilon2, const int& col_weight )
{
  
  //Reading in
  std::ifstream input(filename.c_str());

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
      sourceGalaxy galaxy;
      galaxy.x_=values.at(col_x-1);
      galaxy.y_=values.at(col_y-1);
      galaxy.z_=values.at(col_z-1);

      //Storing Ellipticity and Weight if column number not 0
      if(col_epsilon1!=0) galaxy.epsilon1_ = values.at(col_epsilon1-1);
      if(col_epsilon2!=0) galaxy.epsilon2_ = values.at(col_epsilon2-1);
      if(col_weight!=0) galaxy.weight_ = values.at(col_weight-1);

      galaxies_.push_back(galaxy); //push back galaxy into container
    };

 
  galaxies_.pop_back();//Delete last Galaxy
  
  // Initializing starting node
  root_=NULL;

 
  // Building Tree
  build();
      
  
}

g3lcong::Kdtree::~Kdtree()
{
  if(root_!=NULL) deleteTree(root_);
}

long g3lcong::Kdtree::getSize()
{
  return galaxies_.size();
}

g3lcong::Kdnode* g3lcong::Kdtree::getRoot()
{
  return root_;
}


void g3lcong::Kdtree::build()
{
  //delete Tree if it is already built
  if(root_!=NULL) deleteTree(root_);

  root_ = new Kdnode;

  build(root_, 0, galaxies_.size()-1);
}


void g3lcong::Kdtree::build(Kdnode* pointer, const long& n1, const long& n2)
{
  // Initialization
  pointer->child_left_ = NULL;
  pointer->child_right_ = NULL;
  pointer->start_ = n1;
  pointer->end_ = n2;

  // Average Catalog between n1 and n2
  averageCatalog(pointer, n1, n2);

  // Stop, if no galaxies between n1 and n2 or nodes cannot be further divided
  if(n1==n2 || (pointer->dx_<=0. && pointer->dy_<=0.) )
    {
      return;
    }

  // Subdivide if only two galaxies are left
  if(n2-n1==1)
    {
      pointer->child_left_ = new Kdnode;
      pointer->child_right_ = new Kdnode;
      build(pointer->child_left_, n1, n1);
      build(pointer->child_right_, n2, n2);
      return;
    }

  // Subdivision into two ranges with sizes in1 and in2
  long in1, in2;
  subdivideCatalog(pointer, in1, in2);

  // Create Subnodes for in1 range
  if(in1)
    {
      pointer->child_left_ = new Kdnode;
      build(pointer->child_left_, n1, n1+in1-1);
    };

  // Create Subnodes for in2 range
  if(in2)
    {
      pointer->child_right_ = new Kdnode;
      build(pointer->child_right_, n1+in1, n2);
    }

  
}

void g3lcong::Kdtree::averageCatalog(Kdnode* pointer, const long& n1,
				  const long& n2)
{
  // Get Dimensions of Boundary box and average z, epsilons and weights
  
  double  maxx=0.,minx=0.;
  double  maxy=0.,miny=0.;

  double x=0.,y=0.,z=0.,epsilon1=0, epsilon2=0, weight=0; 

  
  for(long n=n1; n<=n2; n++)
    {
      x += galaxies_.at(n).x_;
      y += galaxies_.at(n).y_;
      z += galaxies_.at(n).z_;
      epsilon1 += galaxies_.at(n).epsilon1_*galaxies_.at(n).weight_;
      epsilon2 += galaxies_.at(n).epsilon2_*galaxies_.at(n).weight_;
      weight += galaxies_.at(n).weight_;


      minx = (n==n1||minx>galaxies_.at(n).x_ ? galaxies_.at(n).x_:minx);
      maxx = (n==n1||maxx<galaxies_.at(n).x_ ? galaxies_.at(n).x_:maxx);
      miny = (n==n1||miny>galaxies_.at(n).y_ ? galaxies_.at(n).y_:miny);
      maxy = (n==n1||maxy<galaxies_.at(n).y_ ? galaxies_.at(n).y_:maxy);

      
    }

  pointer->N_ = n2-n1+1;

  pointer->central_gal_.x_ = x/pointer->N_;
  pointer->central_gal_.y_ =y/pointer->N_;
  pointer->central_gal_.z_ = z/pointer->N_;
  pointer->central_gal_.epsilon1_ = epsilon1;
  pointer->central_gal_.epsilon2_ = epsilon2;
  pointer->central_gal_.weight_ = weight;

  pointer->x1_ = minx;
  pointer->x2_ = maxx;

  pointer->y1_ = miny;
  pointer->y2_ = maxy;

  pointer->dx_ = maxx-minx;
  pointer->dy_ = maxy-miny;

  pointer->center_x_     = pointer->x1_+.5*pointer->dx_;
  pointer->center_y_     = pointer->y1_+.5*pointer->dy_;

  pointer->diameter_ = sqrt(pointer->dx_*pointer->dx_
			    + pointer->dy_*pointer->dy_);
  
  
}


void g3lcong::Kdtree::subdivideCatalog(Kdnode* pointer, long& in1, long& in2)
{
  // Get Start and End of Node
  long head = pointer->start_;
  long tail = pointer->end_;
 
  // Calculate Direction of best separation
  double xx=0., yy=0., xy=0., count=0.;

  for(long i=head; i<= tail; i++)
    {
      xx = xx*count/(1.+count) +
	(galaxies_.at(i).x_-pointer->central_gal_.x_) *
	(galaxies_.at(i).x_-pointer->central_gal_.x_)/(1.+count);
      yy = yy*count/(1.+count) +
	(galaxies_.at(i).y_-pointer->central_gal_.y_) *
	(galaxies_.at(i).y_-pointer->central_gal_.y_)/(1.+count);
      xy = xy*count/(1.+count) +
	(galaxies_.at(i).x_-pointer->central_gal_.x_) *
	(galaxies_.at(i).y_-pointer->central_gal_.y_)/(1.+count);
      count++;
    }

  
  double eigenvalue = 0.5*(xx+yy+sqrt((xx+yy)*(xx+yy)-4.*(xx*yy-xy*xy)));

  double eigenvector_x = eigenvalue - yy;
  double eigenvector_y = xy;

  if(fabs(xy)<1E-10 && fabs(xx)>1E-10)
    {
      eigenvector_x=1.;
      eigenvector_y=0.;
    }
  else if (fabs(xy)<1E-10 && fabs(yy)>1E-10)
    {
      eigenvector_x=0.;
      eigenvector_y=1.;
    }

  // If only one galaxy: seperate along connecting line
  if(tail-head == 1)
    {
      eigenvector_x = galaxies_.at(head).x_ - galaxies_.at(tail).x_;
      eigenvector_y = galaxies_.at(head).y_ - galaxies_.at(tail).y_;
    }


  // Do separation

  sourceGalaxy dummy; 
 
  while(head != tail)
    {
      if(!isLeftChild(head,pointer,eigenvector_x,eigenvector_y))
	{
	  dummy=galaxies_.at(head);
	  galaxies_.at(head) = galaxies_.at(tail);
	  galaxies_.at(tail--)=dummy;
	}
      else
	{
	  head++;
	}

  
    };

  in1 = head - pointer->start_ + isLeftChild(head, pointer, eigenvector_x,eigenvector_y);
  in2 = pointer->end_ - pointer->start_ + 1 - in1;

  
}


bool g3lcong::Kdtree::isLeftChild(const long& n, Kdnode* pointer, const double& eigenvector_x, const double& eigenvector_y)
{
  double dx = galaxies_.at(n).x_ - pointer->central_gal_.x_;
  double dy = galaxies_.at(n).y_ - pointer->central_gal_.y_;

  return ((dx*eigenvector_x + dy*eigenvector_y)>=0);
}


void g3lcong::Kdtree::deleteTree(Kdnode* pointer)
{
  // Check if both children are Leaves
  if (pointer->child_left_==NULL && pointer->child_right_==NULL) return;

  if (pointer->child_left_==NULL)
    {
      deleteTree(pointer->child_left_);
      delete pointer->child_left_;
    }
  if (pointer->child_right_==NULL)
    {
      deleteTree(pointer->child_right_);
      delete pointer->child_right_;
    }
  return;  
}

void g3lcong::Kdtree::randomizeRedshifts(const double& z_max)
{
  std::mt19937 random_number_generator;

  random_number_generator.seed(std::random_device()());

  std::uniform_int_distribution<std::mt19937::result_type> distribution (0,1000);


  for(unsigned int i=0; i<galaxies_.size(); i++)
    {
      galaxies_.at(i).z_ = (double)distribution(random_number_generator)*z_max/1000.;
    };
  
}
