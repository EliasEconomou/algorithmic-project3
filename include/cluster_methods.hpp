#ifndef CLUSTER_METHODS_H
#define CLUSTER_METHODS_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <string.h>
#include <vector>
#include "point_functions.hpp"
#include "hash_table.hpp"
#include "hash_functions.hpp"
#include "algorithms.hpp"


using namespace std;

double silhuette(Cluster_of_points cluster, int i);

double silhuette(Cluster_of_curves cluster, int i);

void calculate_centroids(Cluster_of_points &cluster);

Cluster_of_points initialize_kplusplus(Vector_of_points &Data, Cluster_of_points &cluster, int number_of_clusters);

Cluster_of_points lloyds(Vector_of_points &Data, Cluster_of_points &cluster, int iter_num_input);

Cluster_of_curves lloyds(Vector_of_curves &Data, Cluster_of_curves &cluster, int iter_num_input);

Cluster_of_points cluster_Classic(Vector_of_points &Data, Cluster_of_points &cluster, int number_of_clusters);

Cluster_of_curves cluster_Classic(Vector_of_curves &Data, Cluster_of_curves &cluster, int number_of_clusters);

Cluster_of_points cluster_LSH(Vector_of_points &Data, Cluster_of_points &cluster, int number_of_clusters, int L_of_LSH, int k_of_LSH);

Cluster_of_curves cluster_LSH_Frechet(Vector_of_curves &Data, Cluster_of_curves &cluster, int number_of_clusters, int L_of_LSH, int k_of_LSH);

Cluster_of_points cluster_Hypercube(Vector_of_points &Data, Cluster_of_points &cluster, int number_of_clusters, int M_of_Hypercube, int k_of_Hypercube, int probes_of_Hypercube);

#endif