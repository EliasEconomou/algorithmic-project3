#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <utility>
#include <climits>
#include <float.h>
#include <chrono>
#include "point_functions.hpp"
#include "hash_table.hpp"
#include "cube_table.hpp"
#include "hash_functions.hpp"
#include "grid_table.hpp"

#include "frechet.hpp"


// Compare function to use in set.
struct CompDist
{
    constexpr bool operator()(std::pair<ClassPoint, double> const& a, std::pair<ClassPoint, double> const& b)
    const noexcept
    {
        return a.second < b.second;
    }
};


//TRUE

std::pair<ClassPoint,double> true_NN(ClassPoint q, Vector_of_points inputData, double &time);

std::set<std::pair<ClassPoint,double>, CompDist> true_nNN(ClassPoint q, int N, Vector_of_points inputData, double &time);


//LSH

std::pair<ClassPoint,double> lsh_approximate_NN(ClassPoint q, std::vector<HashTable> hashTables, LSH_hash_info *hInfo, double &time);

std::set<std::pair<ClassPoint,double>, CompDist> lsh_approximate_nNN(ClassPoint q, int N, std::vector<HashTable> hashTables, LSH_hash_info *hInfo, double &time);

std::unordered_map<std::string,double> lsh_approximate_range_search(ClassPoint q, double R, std::vector<HashTable> hashTables, LSH_hash_info *hInfo);


//HYPERCUBE

std::pair<ClassPoint,double> cube_approximate_NN(ClassPoint q, CubeTable cubeTable, CUBE_hash_info *hInfo);

std::set<std::pair<ClassPoint,double>, CompDist> cube_approximate_nNN(ClassPoint q, int N, CubeTable cubeTable, CUBE_hash_info *hInfo, double &time);

std::unordered_map<std::string,double> cube_approximate_range_search(ClassPoint q, double R, CubeTable cubeTable, CUBE_hash_info *hInfo);


//DISCRETE FRECHET

std::pair<ClassCurve,double> lsh_approximate_NN(ClassCurve q, std::vector<GridTable> gridTables, LSH_hash_info *hInfo, double &time);

std::pair<ClassCurve,double> true_NN(ClassCurve q, Vector_of_curves inputData, double &time);

std::unordered_map<std::string,double> lsh_approximate_range_search(ClassCurve q, double R, std::vector<GridTable> gridTables, LSH_hash_info *hInfo);

#endif