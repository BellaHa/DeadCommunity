/*
 * igraphConversion.hpp
 *
 *  Created on: May 22, 2011
 *      Author: Thang N. Dinh (tdinh@cise.ufl.edu)
 */

#ifndef IGRAPHCONVERSION_HPP_
#define IGRAPHCONVERSION_HPP_

#ifdef	IGRAPH
extern "C" {
#include <igraph/igraph.h>
}


#define	iGraph		igraph_t
#define	iVector		igraph_vector_t
#define	iInteger	igraph_integer_t
#define iBool		igraph_bool_t
#define	iReal		igraph_real_t
#define	iMatrix		igraph_matrix_t
#define iSelector 	igraph_vs_t

#include <iostream>
#include "generalGraph.hpp"
#include "edgeMap.hpp"
#include <algorithm>
#include <vector>


using namespace std;

/**
 * @brief Convert a GeneralGraph to an igraph
 *
 */
template <class T=int>
class	IGraph {
public:
	int n;	
    iGraph		g;
    //vector<int> rid; // Reversed id
	//vector<int> id;  // id of a node in the original graph
	vector<iReal> wv; // weight vector
	iVector		  w;
	bool		weighted;
	bool		directed;
	/**
	 * @brief Convert a GeneralGraph to a ListDigraph (IGRAPH graph)
	 * @param defEdgeVal	Value for edges in case the graph is unweighted
	 * @return	A ListDigraph WITH the weight map
	 */
    IGraph(const GeneralGraph<T> &gg) {
    	// TODO: check
    	n = gg.size();
    	weighted = gg.isWeighted();
    	directed = gg.directed;

    	vector<iReal>  edges;
		ALL_EDGES(gg,h,t) {
			if (directed || (h <= t) ){
				edges.push_back(h);
				edges.push_back(t);
				if (weighted)
					wv.push_back( gg.edw(id) );
			}
		}
		iVector	v;
		igraph_vector_view(&v,(iReal*) (&edges[0]), edges.size() );
		igraph_create(&g, &v, 0, directed?IGRAPH_DIRECTED:IGRAPH_UNDIRECTED);

		if (weighted)
			igraph_vector_view(&w, (iReal*)(&wv[0]), wv.size() );
	}

    /**
     * Default constructor: Init an empty graph.
     */
    IGraph() {
    	n = 0;
    	weighted =false;
    	directed =false;
    }


	/**
	 * @brief Construct a GeneralGraph from the igraph
	 * @return
	 */
    GeneralGraph<T>	getBack() {
		// TODO: check
		n = (int)(igraph_vcount(&g));
		vector<int> head, tail;
		vector<T>	ww;
		MY_FOR(i, igraph_ecount(&g) ) {
				iInteger u, v;
				igraph_edge(&g, i,&u, &v);
				//n = max(n, ((int)max(u,v) );
				head.push_back(u);
				tail.push_back(v);
				if (weighted)
					ww.push_back( VECTOR(w)[i]);
		}
		return  GeneralGraph<T> (n, head.size(), &head[0], &tail[0],
				weighted?(&ww[0]):NULL);
	}

	~IGraph() {
		igraph_destroy(&g);
	}
};


template <class T>
inline vector<T> 	iVector_to_vector(iVector *v) {
	int s = igraph_vector_size(v);
	vector<T>	result( s );
	MY_FOR(i, s)
		result[i] = VECTOR(*v)[i];
	return result;
}

template <class T>
inline void vector_to_iVector(const vector<T> &v, iVector *iv) {
	int s = v.size();
	igraph_vector_resize(iv, s);
	MY_FOR(i, s)
		VECTOR(*iv)[i] = v[i];
}

#endif /* IGRAPH */
#endif /* IGRAPHCONVERSION_HPP_ */

