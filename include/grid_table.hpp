#ifndef GRID_TABLE_H
#define GRID_TABLE_H

#include <iostream>
#include <list>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include "point_functions.hpp"
#include "hash_functions.hpp"


class GridNode
{
public:
    long int ID;
    ClassCurve *curve;
    GridNode(ClassCurve *p, long int ID);
};

class GridTable
{
    int bucketsNumber;
    std::vector<std::list<GridNode>> lists; //a list for every bucket
    
public:
    int pointDim; //dimension of grid / points of curves
    int curveDim; //dimension of curves (number of points)
    double tShiftGrid; //shift grid's y/x by this number
    double epsilon; //for filtering (1D grid)
    double delta;
    std::vector<double> t;
    std::vector<std::vector<double> > v; //k vectors to use to compute every h
    std::vector<int> r;

    GridTable(int bucketsNumber, double delta, int curveDim, int pointDim);
    
    // Insert curve in grid
    void GridInsert(ClassCurve *p, LSH_hash_info *hInfo);

    // Display grid (debug)
    void GridDisplay();

    int get_bucketsNumber();

    std::list<GridNode> get_bucketList(int g);
};

#endif