/*
 * main.cpp
 *
 *  Created on: Jun 3, 2010
 *      Author: thang
 */

#define IGRAPH
extern "C" {
#include <igraph/igraph.h>
}

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits>
#include <string>
#include <sstream>
#include <ctime>
#include "generalGraph.hpp"
#include "graphAlgorithm.hpp"
#include "igraphConversion.hpp"
#include "DynGraph.hpp"



// Quick utilities
template <class T>
inline std::string toString(T t) {
  std::ostringstream oss;
  oss << t;
  return oss.str();
} 

#ifndef	 DCHECK
#ifdef	_DEBUG

#define	 DCHECK(a, b) cout<<"Checking "<<__FILE__<<":"<<__LINE__<<": ";cout.flush();\
		if (!(a) ) { cout << b << endl; exit(1); }\
		cout<<"...passed!"<<endl;
#else
#define		DCHECK( a, b)  {}

#endif
#endif

using namespace std;

//#ifdef	_DEBUG
//#include "unweightedGraphTest.h"
//#endif

void test_igraph() {
	igraph_real_t avg_path;
	  igraph_t graph;
	  igraph_vector_t dimvector;
	  igraph_vector_t edges;
	  int i;

	  igraph_vector_init(&dimvector, 2);
	  VECTOR(dimvector)[0]=30;
	  VECTOR(dimvector)[1]=30;
	  igraph_lattice(&graph, &dimvector, 0, IGRAPH_UNDIRECTED, 0, 1);

	  srand(100);
	  igraph_vector_init(&edges, 20);
	  for (i=0; i<igraph_vector_size(&edges); i++) {
	    VECTOR(edges)[i] = rand() % (int)igraph_vcount(&graph);
	  }

	  igraph_average_path_length(&graph, &avg_path, IGRAPH_UNDIRECTED, 1);
	  printf("Average path length (lattice):            %f\n", (double) avg_path);

	  igraph_add_edges(&graph, &edges, 0);
	  igraph_average_path_length(&graph, &avg_path, IGRAPH_UNDIRECTED, 1);
	  printf("Average path length (randomized lattice): %f\n", (double) avg_path);

	  igraph_vector_destroy(&dimvector);
	  igraph_vector_destroy(&edges);
	  igraph_destroy(&graph);

	  return ;
}

int head[] = { 1, 1, 3, 1, 2, 4, 5, 0 };
int tail[] = { 2, 3, 1, 2, 1, 5, 6, 1 };

void	test_DynGraph() {
	int len = 8;
	vector<int> h(head, head+len),
				t(tail, tail+len);
	vector<int> w(len, 1);
	DynGraph<int>	dg(7, h, t, w);
	cout << dg;
	dg.simplify();
	cout << dg;
	dg.delete_edge(5, 6, 1);
	dg.delete_edge(1, 2, 1);
	dg.simplify();
	cout << dg;
	MY_FOR(i, 5)
		dg.insert_edge(1, 4, 1);
	MY_FOR(i, 5)
		dg.insert_edge(1, 5, 1);
	MY_FOR(i, 5)
		dg.insert_edge(5, 6, 1);
	cout << dg;
	MY_FOR(i, dg.size() )
			cout <<i<<"\t";
	cout << endl;
	MY_FOR(i, dg.size() )
		cout << dg.degree(i)<<"\t";
	cout << endl;
	MY_FOR(i, dg.size() )
		cout << dg.out_degree(i)<<"\t";

}
int	main() {
	GeneralGraph<> gg(7, 8, head, tail);
	IGraph<> ig(gg);
	gg = ig.getBack();
	vector<vector<int>  > cc = connected_components(gg);
	cout <<"Hello world: " << cc.size() <<" "<< gg.size()<<endl;
	MY_FOR(i, cc.size() ) {
		cout << cc[i].size()<<" ";
		MY_FOR(j, cc[i].size())
			cout <<cc[i][j] <<" ";
		cout << endl;
	}
	test_igraph();
	test_DynGraph();
    return 0;
}

