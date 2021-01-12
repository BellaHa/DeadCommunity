/*
 * EdgeMap.h
 *
 *  Created on: Nov 1, 2010
 *      Author: Thang N. Dinh
 *		Version: 0.1 11/01/2010
 */

#ifndef EDGEMAP_H_
#define EDGEMAP_H_

#include "generalGraph.hpp"
#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <cmath>

#ifndef	MY_FOR
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

/**
 * @brief	Provide an ordering of directed edges in a general graph.
 * Each edge is associated with an id in range 0..Arcs - 1 where Arcs is the number
 * of directed edge in the graph
 *
 * The class implement a hash table (with linear probing + double hashing) to maintain
 * the association between edge (u, v) and its eid. All operations have expected time
 * complexity O(1)
 */

template <class T=int>
class EdgeMap {
protected:
	int	**idMap;		// idMap[u]: Id hash table for edges incident at u
	int *mapSize;		// mapSize[u]: Size of hash table for edge id of u
	int	*edgeHead;
	GeneralGraph<T> *g;	// Pointer to the original graph
	double hashLoad;	// Maximum hash table load, default is 0.50
	T*  weights;

	/*
	 * @brief	Return the smallest prime larger than x
	 * 			Expected complexity log(x) \sqrt(x)
	 */
	int	 nextPrime(int x);
public:
	const static	int	NA = -1;
	/* @brief	Create the edge map of an graph ignoring the weight
	 */
	EdgeMap(GeneralGraph<T> &g, double load = 0.5);

	/* @brief	Set the load factor for the edge hash table
	 */
	double	setHashLoad(double load);

	/* @brief	Return the id/order associated with the directed edge (u, v)
	 * 			The id is in range 0..(Arcs -1). Return -1 if the edge does not exist.
	 */
	int		edgeId(int u, int v) const;

    bool	isEdge(int u, int v) const { return (edgeId(u, v) != NA); }

	/**
	 * @brief	Pointer to the weight/attributes of an edge
	 * @param u	The head of the directed edge
	 * @param v The tail of the directed edge
	 * @return	The pointer to the edge's weight/attribute or NA (-1)
	 * if the edge does not exist.
	 */
	T* operator() (int u, int v) const;

	/* @brief   Locate the edge by the id given in eid and point to the value
	 * of the head.
	 */
	const int*	getHead(int eid) const;

	/**
	 * @brief	Return the number of edges/arcs in the original graph
	 */
	int		size() const {return g->Arcs();}

	/* @brief   Locate the edge by the id given in eid and point to the value
	 * of the head. Please be aware that if the underlying graph doesn't have
	 * edge attributes/weights of type T, the behavior is undefined.
	 */
	T*  getWeight(int eid) const;

	/* @brief   Locate the edge by the id given in eid and point to the value
	 * of the head.
	 */
	const int*	getTail(int eid) const;

	virtual ~EdgeMap() {
		delete[] edgeHead;
		delete[] mapSize;
		MY_FOR(i, g->Vertices())
			delete[] idMap[i];
		delete[] idMap;
	}
};




/*
 * @brief	Return the smallest prime larger than x
 * 			Expected complexity log(x) * sqrt(x)
 */
template <class T>
int	EdgeMap<T>::nextPrime(int x) {
	int sq = (int)(sqrt((float)x));
	bool prime;
	do {
		x++;
		if (sq * sq < x)
				sq++;
		prime = true;
		for(int i = 2; i <= sq;++i)
			if (x % i == 0) {
				prime = false;
				break;
			}
	} while (!prime);
	return x;
}

/* @brief	Create the edge map of an graph ignoring the weight
	 */
template <class T>
EdgeMap<T>::EdgeMap(GeneralGraph<T> &gr, double load): g(&gr), hashLoad(load) {
	// Check if the graph is weighted
	if (gr.isWeighted() ) {
		weights = gr.Weight(0); // pointer to the attr array
	} else
		weights = NULL;
	int n = g->Vertices();
	edgeHead = new int[g->Arcs()];
	idMap = new int*[n];
	mapSize = new int[n];
	// Initialize heads of directed edges
	int *head = edgeHead;
	MY_FOR(u, n)
		MY_FOR(i, g->degree(u) ) {
			*head = u;
			head++;
		}
	// Create hash tables
#ifdef	_DEBUG
	int hashLookup = 0, maxLookup=0, rLookup;
#endif
	MY_FOR(u, n) {
		int m = (g->degree(u) <= 0)? g->degree(u)
		                            : nextPrime((int) (g->degree(u) / hashLoad)+1);
		mapSize[u] = m;
		idMap[u] = new int[m];
//		vector<int>  mc(m, 0);
		MY_FOR(i, m)
			idMap[u][i] = -1;
//		for (int id = g->neighbors[u] - g->edgeList;
//				id < g->neighbors[u + 1]	- g->edgeList; ++id) {
//					int v = g->edgeList[id];
//					int x = MurmurHash2(&v, sizeof(int), 0) % m;
//					int y = MurmurHash2(&v, sizeof(int), u) % m;
//					mc[x]++;
//					mc[y]++;
//		}
		// TODO: Improve hashing, handle multiple edge, measure variance of hash lookup
		int* edges = (*g)[0];
		for (int id = (*g)[u] - (*g)[0]; id < (*g)[u + 1]- (*g)[0]; ++id) {
			int v = edges[id];
			int x = MurmurHash2(&v, sizeof(int), 0) % m;
			IFDEBUG(hashLookup++;rLookup=1)
			if (idMap[u][x] < 0) {
				idMap[u][x] = id;
			} else {
				x = ( (x << 13) + 11-x) % m; /* changed from 43 to 5003 on Apr 5, 2011 and to 8191 on May 17, 2011 */
				while (idMap[u][x] >= 0) {
					//x = (x + 1) % m;
                    x =  (x + 1) % m;
					IFDEBUG(hashLookup++;rLookup++)
				}
				idMap[u][x] = id;
			}
//			IFDEBUG(hashLookup+=2;rLookup=2)
//			if (idMap[u][y] < 0 && (idMap[u][x]>=0 || mc[y] < mc[x]) ) {
//				idMap[u][y] = id;
//			} else if (idMap[u][x] < 0) {
//				idMap[u][x] = id;
//				IFDEBUG(hashLookup--;rLookup--)
//			} else {
//				int step = 1, r=1, x = y + 1;
//				IFDEBUG(hashLookup++;rLookup++)
//				while (idMap[u][x] >= 0) {
//					x = (x + r) % m;
//					step++; r = (r + step * 2 - 1) % m;
//					IFDEBUG(hashLookup++;rLookup++)
//				}
//				idMap[u][x] = id;
//			}
			IFDEBUG(maxLookup = max(maxLookup, rLookup))
		}
	}	
}

/**
 * @brief	Pointer to the weight/attributes of an edge
 * @param u	The head of the directed edge
 * @param v The tail of the directed edge
 * @return	The pointer to the edge's weight/attribute
 */
template<class T>
inline T* EdgeMap<T>::operator() (int u, int v) const {
	return getWeight( edgeId(u, v));
}

/* @brief	Set the load factor for the edge hash table
 *
 */
template <class T>
double	EdgeMap<T>::setHashLoad(double load) {
	if (load <= 0.8 && load > 0)
		hashLoad = load;
	return hashLoad;
}

/* @brief	Return the id/order associated with the directed edge (u, v)
 * 			The id is in range 0..(Arcs -1). Return -1 if the edge does not exist.
 */
template <class T>
inline int	EdgeMap<T>::edgeId(int u, int v) const {
	DCHECK( (u>=0 && v>=0 && u < g->Vertices() && v < g->Vertices() ),"Error:edgeMap::edgeId: Invalid values specified "<<u <<" "<<v);
	unsigned int m = mapSize[u];
	if (m <= 0)
		return NA;
	int *edges = (*g)[0];
	unsigned int x=MurmurHash2(&v,sizeof(int),0) % m;
    if (idMap[u][x] < 0)
        return NA;
    if (edges[ idMap[u][x] ] == v)
        return idMap[u][x];
    x = ( (x<<13) + 11 -x) % m; /* changed from 43 to 5003 on Apr 5, 2011 */
	while (idMap[u][x] >= 0) {
        if ( edges[ idMap[u][x] ] == v)
				return idMap[u][x];	
		x =  ( x+1) % m;					
    }
	  
	return NA;
}

/* @brief   Locate the edge by the id given in eid and point to the value
 * of the head.
 */
template <class T>
inline const int*	EdgeMap<T>::getHead(int eid) const {
	DCHECK(eid >=0 && eid < g->Arcs(),"edgeMap::edgeHead: Invalid value specified "<<eid)
	return edgeHead + eid;
}

/* @brief   Locate the edge by the id given in eid and point to the value
 * of the tail.
 */
template <class T>
inline const int*	EdgeMap<T>::getTail(int eid) const {
	DCHECK(eid >=0 && eid < g->Arcs(),"edgeMap::edgeTail: Invalid value specified "<<eid)
	return (*g)[0] + eid;
}

/* @brief   Locate the edge by the id given in eid and point to the value
 * of the head. Please be aware that if the underlying graph doesn't have
 * edge attributes/weights of type T, the behavior is undefined.
 */
template <class T>
inline T* EdgeMap<T>::getWeight(int eid) const {
	DCHECK(eid >=0 && eid < g->Arcs() && weights != NULL,"edgeMap::edgeWeight: Invalid value specified or graph is unweighted. Eid: "<<eid)
	return &(g->edw(eid));
}

typedef	EdgeMap<> GraphMap;
#endif /* EDGEMAP_H_ */
