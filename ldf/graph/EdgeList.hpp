/*
 * GeneralGraph.hpp
 *
 *  Created on: Oct 26, 2010
 *      Author: Thang N. Dinh
 *      Version: 0.1
 */

#ifndef EDGE_LIST_H_
#define EDGE_LIST_H_

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
#include <fstream>

using namespace std;

#ifndef	MY_FOR
	#define	ALL(v)	(v).begin(), (v).end()
	#define	MY_FOR(i,n) for(int i =0;i<(int)(n);++i)
	#define CHECK( a, b)  if (!(a) ) { cerr << b << endl; exit(1); }
	#ifdef	_DEBUG
		#define	 DCHECK(a, b) 	if (!(a) ) { cout << "Failed "<<__FILE__<<":"\
							<<__LINE__<<": "<< b << endl; cout.flush();exit(1); }
		#define	IFDEBUG(a)	  a;
	#else
	#define		DCHECK( a, b)  {}
		#define		IFDEBUG(a) {}
	#endif
#endif


// Delete an array and set the pointer to NULL
#define		free_array(x)		if (x) { delete[] (x); (x) = NULL;}




/**
 * @brief	Class that represent a general graph using
 * adjacency list. The underlying graph is  DIRECTED i.e. if you want
 * to represent an undirected graph you have to maintain the symmetry of
 * the "weight matrix" by yourself.
 *
 * The toUndirected() method add the reversed arcs of existing arcs to
 * the graph and make sure the "weight matrix" is symmetric
 *
 * Basic properties:
 * 	+ Vertices:	Number of vertices in the graph. Vertices are numbered 0..Vertices-1
 *  + Arcs:		Number of directed arcs
 *  + Edges = Arcs/2: The number of edge if the graph is undirected
 *
 * Assumption:
 *  + Vertices in the graph are numbered from 0 to n-1 (the # of vertices)
 * Optimized for COMPATIBILITY with unweightedGraph & PERFORMANCE & MEMORY
 *
 */

template<class T = int>
struct EdgeList {
	int n, m; // Number of vertices and edges
	vector<int>	h, t;
	vector<T>	w;
};


template<class T>
ostream& operator<<(ostream &os, const EdgeList<T> &g) {
	os << g.n << " " << g.m;
	os << endl;
	MY_FOR(i, g.m) {
			os << g.h[i]<<" "<<g.t[i]
			   << " " << g.w[i] << endl;
	}
	return os;
}


template<class T>
istream& operator>>(istream &is, EdgeList<T> &g) {
	is >> g.n >> g.m;
	g.h.resize(g.m); g.t.resize(g.m);
	g.w.resize(g.m);
	MY_FOR(i, g.m) {
			is >> g.h[i] >> g.t[i];
			is >> g.w[i];
	}
	return is;
}


/*
template <>
istream& operator>> (istream &is, EdgeList<bool> &g) {
	is >> g.n >> g.m;
	g.h.resize(g.m); g.t.resize(g.m);
	g.w.resize(0);
	MY_FOR(i, g.m) {
			is >> g.h[i] >> g.t[i];
	}
	return is;
}
*/



#endif /* EDGE_LIST_H_ */

