#ifndef POINT_FUNCTIONS_H
#define POINT_FUNCTIONS_H

#include <algorithm>
#include <iostream>
#include <istream>
#include <sstream>
#include <fstream>
#include <vector>
#include <list>
#include <random>
#include <string>


////////////////////////////////////////////////////////    POINT   ////////////////////////////////////////////////////////

class ClassPoint
{
public:
    std::string itemID; //index id of point
    std::vector<double> vpoint; //point-vector of coordinates
    
};

class Vector_of_points
{
public:
    std::vector<ClassPoint> points; //vector of points
};



////////////////////////////////////////////////////////    CURVE   ////////////////////////////////////////////////////////

class ClassCurve
{
public:
    std::string curveID; //index id of curve
    std::vector<ClassPoint> cpoints; //time series contains points
};

class Vector_of_curves
{
public:
    std::vector<ClassCurve> curves; //vector of curves
};



////////////////////////////////////////////////////////    CLUSTERING   ////////////////////////////////////////////////////////

class Cluster_of_points //class to implement clusters easier
{
    public:
    std::vector<ClassPoint> centroids; //vector of centroid points
    std::vector<Vector_of_points> points; //vector of points in each cluster
};

class Cluster_of_curves //class to implement clusters easier
{
    public:
    std::vector<ClassCurve> centroids; //vector of centroid curves
    std::vector<Vector_of_curves> curves; //vector of curves in each cluster
};


class kplusplus_helper //helping class
{
    public:
    std::vector<double> Additive_Square_Sums;
    std::vector<std::vector<double>> Dist_From_Centroids;
    std::vector<double> Minimum_Distances;
    std::vector<ClassPoint> Centroids;
    std::vector<bool> IsCentroid;
};

class kplusplus_helper_curves //helping class
{
    public:
    std::vector<double> Additive_Square_Sums;
    std::vector<std::vector<double>> Dist_From_Centroids;
    std::vector<double> Minimum_Distances;
    std::vector<ClassCurve> Centroids;
    std::vector<bool> IsCentroid;
};



////////////////////////////////////////////////////////    FUNCTIONS   ////////////////////////////////////////////////////////

// Parse dataset and return a vector of dataset's vectors.
Vector_of_curves curve_parsing(std::string fileName, int dim);

// Compute discrete frechet distance between curves.
double discrete_frechet_distance (ClassCurve c1, ClassCurve c2);

// Returns the optimal traversal between two curves.
std::list<std::pair<int,int>> FindOptimalTraversal( ClassCurve curve1, ClassCurve curve2 );

// Return the mean of two curves.
ClassCurve Mean2Curves( ClassCurve c1, ClassCurve c2 );

// Parse dataset and return a vector of dataset's vectors.
Vector_of_points parsing(std::string fileName);

// Compute distance between vectors. L=1 for manhattan, L=2 for euclidian.
double distance(std::vector<double> v1, std::vector<double> v2, int L);

// Compute the inner product between two vectors.
double inner_prod(std::vector<int> v1, std::vector<double> v2);
int inner_prod(std::vector<int> v1, std::vector<int> v2);
double inner_prod(std::vector<double> v1, std::vector<double> v2);
double inner_prod(std::vector<double> v1, std::vector<int> v2);

// Returns a random integer in the specified range.
int random_number(int begin, int end);

// Returns a random DOUBLE in the specified range.
double random_double(double n1, double n2);

// Returns modulo of two numbers.
long int modulo(long int a, long long int b);


#endif