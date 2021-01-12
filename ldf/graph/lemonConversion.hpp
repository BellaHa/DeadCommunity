/*
 * lemonConversion.hpp
 *
 *  Created on: May 16, 2011
 *      Author: thang
 */

#ifndef LEMONCONVERSION_HPP_
#define LEMONCONVERSION_HPP_

#define	_LEMON
#include <iostream>
#include <lemon/list_graph.h>
#include "generalGraph.hpp"
#include "edgeMap.hpp"
#include <algorithm>
#include <vector>


using namespace lemon;

/**
 * @brief Convert a GeneralGraph to a ListDigraph (LEMON graph)
 *
 */
template <class T>
class	LemonGraph {
public:
	int n;	
    ListDigraph		g;
    ListDigraph::NodeMap<int> rid; // Reversed id
	ListDigraph::ArcMap<T>	w;	// Edge weight map
	vector<ListDigraph::Node> id; // id of a node in the original graph

	

	/**
	 * @brief Convert a GeneralGraph to a ListDigraph (LEMON graph)
	 * @param defEdgeVal	Value for edges in case the graph is unweighted
	 * @return	A ListDigraph WITH the weight map
	 */
    LemonGraph(const GeneralGraph<T> &gg, T defEdgeVal=1): rid(g), w(g){
		n = gg.size();
		id.resize(n);
		MY_FOR(i, n)
			id[ i ] =  g.addNode();			    
        MY_FOR(i, n) 
        	rid.set( id[i], i);        
		ALL_EDGES(gg,h,t) {
			ListDigraph::Arc e = g.addArc( this->id[h], this->id[t]);
			if (gg.isWeighted() )
				w[e] = gg.edw(id);
			else
				w[e] = defEdgeVal;
		}
	}

	/**
	 * @brief Construct a GeneralGraph from the ListDigraph
	 * @return
	 */
	GeneralGraph<T>	getBack() {
		// Check for new nodes and assign new ids if necessary
		for (ListDigraph::NodeIt u(g); u != INVALID;++u)
			if ( rid[ u ] <= 0 && u != id[0]) {
				rid[ u ] = n;
				++n;
			}
		int m = countArcs(g), i=0;
		vector<int> head(m), tail(m);
		vector<T>	wg(m);
		for(ListDigraph::ArcIt e(g); e !=  INVALID;++e,++i) {
			head[i] = rid[ g.source(e)];
			tail[i] = rid[ g.target(e)];
			wg		= w[ e ];
		}
		return GeneralGraph<T>(n, m, &head[0], &tail[0], &wg[0]);
	}
};


#endif /* LEMONCONVERSION_HPP_ */
