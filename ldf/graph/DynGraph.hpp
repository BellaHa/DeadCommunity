/*
 * DynGraph.hpp
 *
 *  Created on: Oct 26, 2010
 *      Author: Thang N. Dinh
 *      Version: 0.1
 */

#ifndef DYNAMIC_GRAPH_H_
#define DYNAMIC_GRAPH_H_

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

template<class T>
struct SumHolder {
	typedef T TSum;
};

template<>
struct SumHolder<char> {
	typedef unsigned int TSum;
};

template<>
struct SumHolder<int> {
	typedef long long TSum;
};

template<>
struct SumHolder<float> {
	typedef double TSum;
};

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
template<class T = char>
class DynGraph {
private:
	vector<int> mask; 	// Temporary mask array
	int			rMarker;	// Rotating marker for compression
	vector<int> tList;// Temporary weight list
	bool		updating_flag; // If true, then updating is in progress. No auto-compression.
protected:
	int n, nArcs; // Number of vertices, number of arcs
	vector<vector<int> > nb; // Outgoing neighbors
	vector<vector<T> > wg; // Weight vectors
	vector< typename SumHolder<T>::TSum > in_deg, out_deg; // Weighted degrees
	vector<int> next_compress; // The size of adjacency list that triggers auto-compression

	/**
	 * Change the number of vertices
	 */
	void resize(int ns);
public:

#define	auto_compress  1.5 // Threshold for  adjacency list auto-compression
	/*
	 * Default constructor: initialize an empty graph;
	 */
	DynGraph<T>():n(0), nArcs(0), rMarker(0), updating_flag(false) {	resize(0); };

	/*
	 * Construct the graph from list of edges (u, v, wg)
	 * The vectors head, tail, weight are assumed to be of same length.
	 */
	DynGraph<T>& build_graph(int nVertices, const vector<int> &head,
			const vector<int> &tail, const vector<T> &weight);

	/*
	 * Construct the graph from list of edges (u, v, wg)
	 * The vectors head, tail, weight are assumed to be of same length.
	 */
	DynGraph<T>(int nVertices, const vector<int> &head, const vector<int> &tail,
			const vector<T> &weight):updating_flag(false) {				
		build_graph(nVertices, head, tail, weight);
	}

	/*
	 * Read the adjacent list of an unweighted graph from a file
	 * File format:
	 * 		Firstline:			#Vertices	#Edges
	 * 		Next #Edges lines:
	 * 			- if weighted = false, each line contains two integer number representing an edge
	 * 			- if weighted = true, each line contains  3 numbers (u, v, wg) that mean an edge
	 * 			(u, v) with REAL weight wg.
	 * Return false if error
	 */
	DynGraph<T>& read_file(const char* fileName, bool weighted = false, bool directed=true);

	/*
	 * return number of vertices in graph
	 */
	inline int size() const {	return n;	}

	void	begin_update() { updating_flag = true;}
	void	end_update() { updating_flag = false;}


	/*
	 * return number of arcs in graph
	 */
	inline int Arcs() const {	return nArcs;	}

	/*
	 * Return the number of adjacent vertices to a vertex u
	 */
	inline int degree(int u) const { return nb[u].size();	}


	/*
	 * Return the list of adjacent vertices of a vertex u
	 */
	inline vector<int>& operator[](int u) { return nb[u];	}

	/*
	 * Return the list of adjacent vertices of a vertex u
	 */
	inline  const vector<int>& operator[](int u) const { return nb[u];	}



	/*
	 * Return the list of adjacent vertices of a vertex u
	 */
	inline vector<T>& weight(int u) {	return wg[u];	}

	/*
	 * Return the list of adjacent vertices of a vertex u
	 */
	const vector<T>& weight(int u) const {	return wg[u];	}

	/**
	 * Aggregate multiple edges into one by summing up the weights
	 * then remove all edges with zero weights
	 * Time complexity: O( deg(u) )
	 */
	void compress_adjacent_list(int u, bool aggregate=true, bool loop=true);

	/*
	 * Add a new edge (u, v) with weight w to the graph
	 */
	void insert_edge(int u, int v, T w);

	/*
	 * Remove edge (u, v) with weight w from the graph
	 */
	void delete_edge(int u, int v, T w);

	/*
	 * Add all the edges in list from the graph
	 * Note that: if (u,v) has weight w1 in the graph and
	 * we add edge (u,v) with weight w2 the new weight will be w1 + w2.
	 *
	 */
	void insert_edges(const vector<int> &head, const vector<int> &tail,
			const vector<T> &weight);

	/*
	 * Remove duplicate edges and loops in graph
	 */
	void simplify(bool aggregate = true, bool loop = true);

	/**
	 * Return the edge list
	 * @param head
	 * @param tail
	 * @param wg
	 */
    void get_edge_list(vector<int> &head, vector<int> &tail, vector<T> &wg) const;
    T out_degree(int u) { return out_deg[u]; }

	/*
	 * Remove a vertex and all incident edges
	 * from the graph
	 */
	void delete_vertex(int	n) { cerr<<"TODO:"; }

	/**
	 * Input the graph into an outstream
	 */
	friend ostream& operator<<(ostream& os, const DynGraph<T>& g) {
		os << g.n << "\t" << g.nArcs << endl;
		MY_FOR(u, g.n) {
			os << u;
			MY_FOR(i, g.nb[u].size() )
				//os << "\t" << g.nb[u][i] << "\t" << g.wg[u][i];
				 os << "\t(" << g.nb[u][i] << ", " << g.wg[u][i]<<") ";
			os << endl;
		}
		return os;

	}
};



template<class T> void DynGraph<T>::get_edge_list(
		vector<int> &head, vector<int> &tail, vector<T> &w) const
{
    head.resize(nArcs); tail.resize(nArcs); w.resize(nArcs);
	int it_bg=0, it_ed;
    MY_FOR(u, n) {
		it_ed = it_bg + nb[u].size();
		std::fill(head.begin()+it_bg, head.begin()+it_ed, u);
		std::copy(nb[u].begin(), nb[u].end(), tail.begin()+it_bg);
		std::copy(wg[u].begin(), wg[u].end(), w.begin()+it_bg);    	
		it_bg = it_ed;
	}
/*
	head.clear();tail.clear(); w.clear();
	head.reserve(nArcs); tail.reserve(nArcs); w.reserve(nArcs);	
    MY_FOR(u, n) {
		MY_FOR(i, nb[u].size() ) {
    		head.push_back(u);
    		tail.push_back(nb[u][i]);
    		w.push_back(wg[u][i]);
    	}		
	}
*/
}



template <class T>
inline void DynGraph<T>::insert_edge(int u, int v, T w) {
#undef max
	if (u >= n || v >= n)
		resize(std::max(u, v) + 1);
	nb[u].push_back(v);
	wg[u].push_back(w);
	out_deg[u] += w;
	in_deg[v] += w;	
	if ( (!updating_flag) && (nb[u].size() > next_compress[u]) )
		compress_adjacent_list(u);
	++nArcs;
}


template <class T>
inline void DynGraph<T>::delete_edge(int u, int v, T w) {
#undef max
	if (u >= n || v >= n)
		resize(std::max(u, v) + 1);
	nb[u].push_back(v);
	wg[u].push_back(-w);
	out_deg[u] -= w;
	in_deg[v] -= w;	
	if ( (!updating_flag) && (nb[u].size() > next_compress[u]) )
		compress_adjacent_list(u);
	++nArcs;
}


/*
 * Construct the graph from list of edges (u, v, wg)
 * The vectors head, tail, weight are assumed to be of same length.
 * Time complexity: O( |V| + |E|)
 */
template <class T>
DynGraph<T>& DynGraph<T>::build_graph(int nVertices, const vector<int> &head,
		const vector<int> &tail, const vector<T> &weight) {
	// Clear any existing data
	*this = DynGraph<T>();
	// Build new graph
	n = nVertices;
	nArcs = 0;
	this->resize(n);
	rMarker = n;
	// Pre-allocate memory
	MY_FOR(i, head.size() )
		++next_compress[ head[i] ];
	MY_FOR(u, n) {
		nb[u].reserve( next_compress[u]);
		wg[u].reserve( next_compress[u]);
	}
	insert_edges( head, tail, weight);
	return *this;
}


/* Read the adjacent list of an unweighted graph from a file
 * File format:
 * 		Firstline:			#Vertices	#Edges
 * 		Next #Edges lines:
 * 			- if weighted = false, each line contains two integer number representing an edge
 * 			- if weighted = true, each line contains  3 numbers (u, v, wg) that mean an edge
 * 			(u, v) with REAL weight wg.
 * Return false if error
 */
template <class T>
DynGraph<T>& DynGraph<T>::read_file(const char* fileName, bool weighted, bool directed) {
	ifstream f(fileName);
	DCHECK(f, "Error: Cannot open file " << fileName )
	int nEdges;
	f >> n >> nEdges;
	IFDEBUG(
			cout << "Reading ... " << fileName << "\t" << nVertices
			<< " vertices, " << nEdges << " edges." << endl;
	)
	nArcs = nEdges;
	if (!directed)
		nArcs = 2*nEdges;
	vector<int> head(nArcs), tail(nArcs);
	vector<T>	w(nArcs);
	int u, v;
	MY_FOR(i, nEdges) {
		f >> head[i] >> tail[i];
		if (weighted)
			f >> w[i];
		else {
			if (numeric_limits<T>::is_bounded)
				w[i] = 1;
			f.ignore(10000, '\n');
		}
		DCHECK(u >= 0 && v >= 0 && u < n && v < n,
				"Error: Invalid vertex id (" << u << "," << v << ") ")
		if (!directed) {
			head[i + nEdges] = tail[i];
			tail[i + nEdges] = head[i];
			w[i + nEdges ] =   w[i];
		}
	}
	f.close();
	return build_graph(head, tail, w);
}

template <class T>
void DynGraph<T>::insert_edges(const vector<int> &head, const vector<int> &tail,
			const vector<T> &weight) {
#undef max
		int ns = std::max(*max_element(head.begin(), head.end()),
				*max_element(tail.begin(), tail.end())	);
		if (ns >= n)
			this->resize(ns+1 );
		MY_FOR(i, head.size() )
			insert_edge(head[i], tail[i], weight[i]);
}


#define		AUTO_COMPRESS(x)	((x)+((x)>>1)+3) // x*1.5+3
/**
 * Aggregate multiple edges into one by summing up the weights
 * then remove all edges with zero weights
 * Time complexity: O( deg(u) )
 */
template <class T>
void DynGraph<T>::compress_adjacent_list(int u, bool aggregate, bool loop) {
	if (rMarker >= n) {
		/*
		MY_FOR(i, n)
				mask[i] = -1;
		*/
		std::fill(ALL(mask), -1);
		rMarker = 0;
	} else ++rMarker;
	if (!loop)
		mask[u] = rMarker;
	int du= nb[u].size();
	MY_FOR(j, du ) {
		int v = nb[u][j];
		if (mask[v] != rMarker) {
			mask[v] = rMarker;
			tList[v] = j;
		} else {
			if (aggregate && (loop || (u != v)))
				wg[u][ tList[v] ] += wg[u][j];
			wg[u][j] = 0; // Set edges' weights to zero
		}
	}
	// Remove all edges with zero weights
	static T EPS = numeric_limits<T>::epsilon();
	for (int j = du - 1; j >= 0; --j)
		if ((wg[u][j] <= EPS) && (wg[u][j] >= -EPS)) {
			wg[u][j] = wg[u].back();
			wg[u].pop_back();
			nb[u][j] = nb[u].back();
			nb[u].pop_back();
		}
	nArcs -= (du - nb[u].size() );
	next_compress[u] = AUTO_COMPRESS(nb[u].size() );
}

/*
 * Remove duplicate edges and loops in graph
 * @param	aggregate If true, multiple edges will be merged into
 * one edge by summing the weights. Otherwise, only the first edge
 * encountered will be kept.
 * @param	loop	If false, self-loops will be discarded
 * Time complexity: O(|E|+|V|)
 */
template <class T>
void DynGraph<T>::simplify(bool aggregate, bool loop) {
	MY_FOR(u, n)		
		if (next_compress[u] < AUTO_COMPRESS(nb[u].size() ))
			compress_adjacent_list(u, aggregate, loop);
}

/**
 * Change the number of vertices in the graph
 * Time complexity: O( |V| )
 */
template <class T>
void DynGraph<T>::resize(int ns) {
	DCHECK(ns >=0," invalid size "<<ns)
	CHECK(ns >= n," TODO: "<<ns <<" "<< n )
	nb.resize(ns);
	wg.resize(ns);
	tList.resize(ns);
	mask.resize(ns, -1);
	in_deg.resize(ns, 0);
	out_deg.resize(ns, 0);
	next_compress.resize(ns, 0);
	n = ns;
}

#endif /* DYNAMIC_GRAPH_H_ */

