#ifndef G3LCONG_KDTREE_H
#define G3LCONG_KDTREE_H

#include <vector>
#include <string>

#include "helpers.h"



namespace g3lcong
{

  
  /**
   * Struct for Nodes of the Kdtree
   */
  struct Kdnode
  {
    ///beginning of node in catalog
    long start_;
    ///end of node in catalog
    long end_;
    ///number of galaxies in node
    long N_;

    ///size of node
    double x1_, x2_, y1_, y2_, dx_, dy_, diameter_;
    double center_x_, center_y_;

    ///Average galaxy in node (has shear)
    sourceGalaxy central_gal_;

    ///Subnodes
    Kdnode *child_left_, *child_right_;
  };
    
  const double RECTOLERANCE = 0.1;

   /**
   * Class for Kdtrees (two-dimensional)
   * Based on "kdTree.h" from Patrick Simon
   *
   * @author Laila Linke llinke@astro.uni-bonn.de
   */
  class Kdtree
  {
  public:

    ///Galaxy catalog
    std::vector<sourceGalaxy> galaxies_;

    /**
     * Constructor from File
     *
     * @param filename File with Galaxy Catalog
     * @param total_columns Number of Columns in Galaxy Catalog
     * @param col_x Column with X in arcmin
     * @param col_y Column with Y in arcmin
     * @param col_z Column with Redshift
     * @param col_epsilon1 Column with Epsilon 1
     * @param col_epsilon2 Column with Epsilon 2
     * @param col_weight Column with Weight
     */
    Kdtree(const std::string& filename, const int& total_columns,
	   const int& col_x, const int& col_y, const int& col_z,
	   const int& col_epsilon1, const int& col_epsilon2,
	   const int& col_weight, bool flipE1=false, bool flipE2=false );

    ///Destructor
    ~Kdtree();

    ///Returns number of galaxies in tree
    long getSize();

    ///Returns Root of Tree
    Kdnode* getRoot();

    void randomizeRedshifts(const double& z_max);
    
  private:

    ///First Node of Tree
    Kdnode* root_;

    ///Build Trees
    void build();

    /**
     * Recursively builds Tree starting from pointer
     *
     * @param pointer Starting Node
     * @param n1 Index of first galaxy in tree
     * @param n2 Index of last galaxy in tree
     */
    void build(Kdnode* pointer, const long& n1, const long& n2);

    /**
     * Average Galaxies in Node
     *
     * @param pointer Node that is to be averaged
     * @param n1 Index of first galaxy in node
     * @param n2 index of last galaxz in node
     */
    void averageCatalog(Kdnode* pointer, const long& n1, const long& n2);

    /**
     * Subdivides Catalog
     * @param pointer Node from which we subdivide
     * @param in1 Number of Galaxies in first subnode
     * @param in2 Number of Galaxies in second subnode
     */
    void subdivideCatalog(Kdnode* pointer, long& in1, long& in2);


    /**
     * Returns whether galaxy at index is left
     * @param n index of galaxy
     * @param pointer node to be divided
     * @param ex x-coordinate of axis
     * @param ey y-coordinate of axis
     */
    bool isLeftChild(const long& n, Kdnode* pointer, const double& ex, const double& ey);

    /**
     * Recursively deletes Tree
     * @pointer Starting Node
     */
    void deleteTree(Kdnode* pointer);


   
    
  };

  




  
}









#endif //G3LCONG_KDTREE_H

