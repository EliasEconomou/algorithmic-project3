#ifndef SEARCH_METHODS_H
#define SEARCH_METHODS_H

#include <iostream>
#include <string>
#include "point_functions.hpp"
#include "hash_table.hpp"
#include "cube_table.hpp"
#include "hash_functions.hpp"
#include "algorithms.hpp"
#include "grid_table.hpp"


// Function for i) assigment using LSH.
void time_series_LSH(std::string inputFile, std::string queryFile, std::string outputFile, int k, int L);

// Function for i) assigment using Hypercube.
void time_series_Hypercube(std::string inputFile, std::string queryFile, std::string outputFile, int k, int M, int probes);

// Function for ii) assigment using Discrete Ferchet.
void time_series_DiscreteFrechet(std::string inputFile, std::string queryFile, std::string outputFile, int k, int L, double delta, bool FrechetBrute);

// Function for iii) assigment using Continuous Frechet.
void time_series_ContinuousFrechet(std::string inputFile, std::string queryFile, std::string outputFile, int k, int L, double delta, bool FrechetBrute);

#endif