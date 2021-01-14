#pragma once

#include "rwgraph.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <cstring>
#include <random>
#include "mappedHeap.hpp"
#include "HeapData.hpp"

#if defined(_OPENMP)

#include <omp.h>

#else
typedef int omp_int_t;
inline omp_int_t omp_set_num_threads(int t) { return 1;}
inline omp_int_t omp_get_thread_num() { return 0;}
#endif

// using namespace std;

/*
* building the hypergraph procedure which generates hyperedges following LT model
*/
long long addHyperedge(Graph &g, HyperGraph &hg, int t, long long num, bool lt);

/*
* find seed nodes procedure using greedy algorithm
*/
void buildSeedSet(HyperGraph &hg, vector<int> &seeds, unsigned int n, int k, vector<double> &degree);
