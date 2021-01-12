/*
 * graphAlgorithm.hpp
 *
 *  Created on: Nov 7, 2010
 *      Author: Thang N. Dinh
 */

#ifndef GRAPHALGORITHM_HPP_
#define GRAPHALGORITHM_HPP_

#include  "generalGraph.hpp"
#include  "mappedHeap.hpp"
#include  "edgeMap.hpp"
#include  <vector>
#include  <algorithm>
#include  <iostream>
#include  <sstream>
#include  <limits>

using namespace std;


template <class T>
struct largerDistance { // Compare operator for Dijkstra's algorithm
    static T *value;
	bool operator()(int &x, int &y) const {
		return value[x] > value[y];
	}
};

template <class T>
T* largerDistance<T>::value= NULL;

/**
 * Finding the shortest path between two vertices s, t using
 * Dijkstra's algorithm.
 * @param	s	The source node
 * @param	t	The destination, if t < 0, then find shortest
 * paths from s to all vertices in the graph
 * @param	g	The graph with n vertices and non-negative edge weights
 * @param	d	Pointer to the distance array of size n
 * If NULL then a temporary array will be allocated.
 * @param	prev	Pointer to the trace array of size n
 * If NULL then a temporary array will be allocated.
 * @return	Distance from s to t or, if t <0, the distance to the
 * farthest vertex from s
 */
template <class T>
T dijkstra(const GeneralGraph<T> &g, int s, int t,  T* d=NULL, int *prev=NULL) {
	return dijkstra<T, T>(g, &(g.edw(0)), s, t, d, prev);
}

/**
 * Finding the shortest path between two vertices s, t using
 * Dijkstra's algorithm.
 * @param	s	The source node
 * @param	t	The destination, if t < 0, then find shortest
 * paths from s to all vertices in the graph
 * @param	g	The graph with n vertices and non-negative edge weights
 * @param	wg	The weight vector.
 * @param	d	Pointer to the distance array of size n
 * If NULL then a temporary array will be allocated.
 * @param	prev	Pointer to the trace array of size n
 * If NULL then a temporary array will be allocated.
 * @return	Distance from s to t or, if t <0, the distance to the
 * farthest vertex from s
 */
template <class TG, class T>
T dijkstra(const GeneralGraph<TG> &g, const T *w, int s, int t,  T* d=NULL, int *prev=NULL) {
	bool  td = (d), tprev = (prev); // Memory allocation
	int  n = g.Vertices();
	if (!td)
		d = new T[ n ];
	if (!tprev)
		prev = new int[n];
	vector<int> vt(n, 0);
	MY_FOR(u, n) {
		d[u] = numeric_limits<T>::max();
		prev[u] = -1;
		vt[u ] = u;
	}
	d[s] = 0;
	largerDistance<T>::value = &d[0];
	MappedHeap<largerDistance<T> >	 queue(vt);
	int u;
	while ( !queue.empty()) {
        u = queue.pop();
        if (u == t || d[u] >= numeric_limits<T>::max())
            break;
        NEIGHBORS(g, u, v) {
            if (d[v] > d[u] + w[id]) {
                DCHECK(wg > 0,"Error: Negative weight detected on edge "<<u <<"-"<<v<<": "<<wg[id]<<endl);
                d[v] = d[u] + w[id];
                prev[v] = u;
                queue.heapify(v);
           }
        }
	}
    if (!td)
        delete[] d;
    if (!tprev)
        delete[] prev;
    if (t >= 0 && t < n)
    	return  d[t];
    else return d[u]; // Last considered vertex
}


/**
 * Remove edges from a SIMPLE graph. This is NOT a memory efficient
 * implementation. It requires as much space as 4 copies of the original graph
 * in memory.
 * @param g		Original graph
 * @param edges Edges/arcs to delete
 * @return	Graph after removing edges
 */
template <class T>
GeneralGraph<T> removeEdges(GeneralGraph<T> &g, vector< pair<int, int> > edges) {
	EdgeMap<T> gm(g);
	int m = g.Arcs();
	int	*tHead = new int[m];
	memcpy(tHead, gm.getHead(0), sizeof(int) * m);
	int *tTail = new int[m];
	memcpy(tTail, gm.getTail(0), sizeof(int) * m);
	T	*tW	   = new  T[m];
	memcpy(tW, gm.getWeight(0), sizeof(int) * m);
	// Mark deleted edges
	MY_FOR(i, edges.size() ) {
		int eid = gm.edgeId( edges[i].first, edges[i].second);
		DCHECK(eid >=0,"Warning: remove non-existing edge ("
				        << edges[i].first<<","<<edges[i].second<<")"<<endl);
		if (eid >= 0) {
			tHead[eid] = tTail[eid] = -1;
		}
	}
	// Move deleted edges to the end
	for(int i = m-1;i>=0;--i)
	if (tHead[i] < 0){
		--m;
		tHead[i] = tHead[m];
		tTail[i] = tTail[m];
		tW[i]	 = tW[m];
	}
	// Create a new graph from remaining edges
	GeneralGraph<T> result(g.Vertices(), m, tHead, tTail, tW);
	delete[] tHead;
	delete[] tTail;
	delete[] tW;
	return result;
}



/**
 * @brief Find the connected component of a node 
 * @param start  Node to find connected component 
 * @param marked Array that marks visited nodes
 * @param[out] A list of nodes in the connected component of 'start'
 */
template <class T>
vector<int>	 bfs(const GeneralGraph<T> &g, int start, vector<bool> &marked) {
	vector<int> queue(1, start);	
	marked[start] = true;
	int top = 0, last = 0;
	while (top <= last) {
		int u = queue[top];top++;
		NEIGHBORS(g, u, v)
			if (!marked[v]) {
				marked[v] = true;
				queue.push_back(v);
				last++;
			}
	}
	return queue;
}

/**
 * @brief Find the connected component of a node 
 * @param start  Node to find connected component 
 * @param black_list  A list of node that should not be considered
 * @param[out] A list of nodes in the connected component of 'start'
 */
template <class T>
vector<int>	 bfs(const GeneralGraph<T> &g, int start, vector<int> black_list=vector<int>(0)) {
	int n = g.size();
	vector<int> queue(1, start);
	vector<bool> marked(n, false);
	MY_FOR(i, black_list.size())
		marked[ black_list[i] ] = true;
	marked[start] = 1;
	int top = 0, last = 0;
	while (top <= last) {
		int u = queue[top];top++;
		NEIGHBORS(g, u, v)
			if (!marked[v]) {
				marked[v] = true;
				queue.push_back(v);
				last++;
			}
	}
	return queue;
}


/**
 * @brief Return the connected components of the given UNDIRECTED graph 'g'.
 * @param[in] visisted Visisted vertices will not be considered
 * @param[out] A two dimensional vector. Vertices inside a
 * connected components is listed as a vector of integers.
 */
template <class T> 
vector< vector<int> > connected_components(const GeneralGraph<T> &g, 
					vector<int> black_list=vector<int>(0)) {
	int n = g.size();
	vector< vector<int> > result;
	vector<int> marked(n, 0);
	MY_FOR(i, black_list.size() )
		marked[ black_list[i] ]  = 1;
	vector<int> q(n, 0);
	MY_FOR(u, n)
		if (!marked[u]) {
			int first = 0, last = 0, w;
			q[ first ] = u; marked[u] = 1;
			while (first <= last) {
					w = q[first]; first++;
					NEIGHBORS(g, w, v) 
						if (!marked[v]) {
							marked[v] = 1;
							last++;
							q[last] = v;
						}
				} // while
			result.push_back( vector<int>(q.begin(), q.begin()+last+1) );
		} // if
	return result;
}


/**
 * @brief Return the connected subgraphs of the given UNDIRECTED graph 'g'.
 * @param[in] visisted Mask for vertices that will  be ignored
 * @param[out] A vector of subgraphs that are connected components of g.
 */
template <class T>
vector< GeneralGraph<T> > connected_subgraphs(const GeneralGraph<T> &g,
					vector<int> black_list=vector<int>(0)) {
	// TODO: Debug this!
	int n = g.size();
	vector< GeneralGraph<T> > result;
	vector<int> marked(n, 0);
	MY_FOR(i, black_list.size() )
		marked[ black_list[i] ]  = 1;
	vector<int> q(n, 0);
	MY_FOR(u, n)
	if (!marked[u]) {
		vector<int> head, tail;
		vector<T> weight;
		int first = 0, last = 0, w, vid=0;
		q[ first ] = u; marked[u] = 1;
		idLookupTable vertexId;
		vector<int> names;
		while (first <= last) {
			w = q[first]; first++;
			if (vertexId.find(w) == vertexId.end()) {
				vertexId[w] = vid;
				names.push_back(w);
				vid++;
			}
			NEIGHBORS(g, w, v) {
				head.push_back(w);
				tail.push_back(v);
				// Add new vertex id if necessary
				if (vertexId.find(v) == vertexId.end()) {
					vertexId[v] = vid;
					names.push_back(v);
					vid++;
				}
				head.push_back(vertexId[w]);
				tail.push_back(vertexId[v]);
				if (g.isWeighted() )
					weight.push_back(g.edw(id));
				if (!marked[v]) {
					marked[v] = 1;
					last++;
					q[last] = v;
				}
			}
		} // while
		T* wp = (g.isWeighted())? (&weight[0]):NULL;
		GeneralGraph<T> gs(vid, head.size(), &head[0], &tail[0], wp);
		gs.vid = names;
		result.push_back( gs);
	} // if
	return result;
}



#endif /* GRAPHALGORITHM_HPP_ */
