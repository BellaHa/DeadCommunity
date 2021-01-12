/*
 * GeneralGraph.hpp
 *
 *  Created on: Oct 26, 2010
 *      Author: Thang N. Dinh
 *      Version: 0.1
 */

#ifndef GENERAL_GRAPH_H_
#define GENERAL_GRAPH_H_

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


/**
 * @brief	Loops through neighbors of a vertex u in the graph g. Pointer v points
 * to the adjacent vertex.
 */
#define	NB(g, u, v)	for(int* v=(g)[u],ev=(g)[u+1]; v != ev;++v)

// Delete an array and set the pointer to NULL
#define		free_array(x)		if (x) { delete[] (x); (x) = NULL;}




#define	ALL_EDGES(graph, h, t) for(int h=0,id=0, t=(graph).tail(id), nvt=(graph).Vertices();h<nvt;++h)\
										for(int _tmp_j=(graph).degree(h);(_tmp_j>0)&&( (t=(graph).tail(id))>=0);--_tmp_j,++id)
#define	NEIGHBORS(graph, h, t) for(int id=(graph)[h]-(graph)[0], t=*((graph)[h]), mid= ((graph)[(h)+1]-(graph)[0]);(id < mid)&&( (t=(graph).tail(id))>=0);++id)
//#define	NEIGHBORS(graph, h, t) for(int id=(graph)[h]-(graph)[0], t=*((graph)[h]);(id < ((graph)[(h)+1]-(graph)[0]))&&( (t=(graph).tail(id))>=0);++id)
//#define	NEIGHBORS(graph, h, t) for(int id=(graph)[h]-(graph)[0];(id < ((graph)[(h)+1]-(graph)[0]))&&( (t=(graph).tail(id))>=0);++id)




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
class GeneralGraph {
protected:
	int vertices, arcs; // Number of vertices, edges, arcs in the graph
	// edges = arcs/2
	int *vDegrees, *inDegree; // Contain (out) degrees and in degrees of vertices in the graph.
	T *wDegree, *wIndegree;// Contain weighted (out) degrees and in degrees of vertices in the graph.
	int **neighbors; // neighbors[u]: pointer to the beginning of
	// neighbors of u in the edgeList. neighbors[vertices]
	// points to the end of edgeList (edgeList[ edges+1]).
	int *edgeList;
	T *attr; // Arcs' attributes
	T **nbw; // nbw[u]: pointer to the beginning of the weight list
	// of edges incident at u, matching the neighbors[u]
	void graphFromEdgeList(int vertices, int arcs, const int* head,
			const int* tail, const T* attribute = NULL);
	void duplicate(const GeneralGraph<T> &g);
public:
	bool directed;
	vector<int>	vid; // Vertices' ids
	/**
	 * @brief Read the graph from a file
	 * Based on the fileName extension the corresponding format will be used.
	 * Available formats are:
	 * 			Extension			Description
	 * 			.adj				+ Text file, adjacency list. The 1st line contains # vertices,
	 * 								# of edges. Next lines contain heads and tails of edges/arcs followed by
	 * 								the edge's attributes.
	 * 								+ Lines start with a letter will be ignored. Hence,
	 * 								special characters such as #,@,$, .etc can be used
	 * 								for comments purpose
     *                              + Vertices are ZERO-based.   
	 * 			.bin				+ Binary file contains numbers of type int.
	 * 								The first number is #vertices (n).
	 * 								Next n numbers (d(0),d(1),...,d(n-1))
	 * 								are degree of vertices.
	 * 								+ Next d(0) numbers are
	 * 								neighbors of 0; next d(1) numbers are neighbors
	 * 								of 1, and so on.
	 * 								+ Next d(0) records are attributes of edges incident at 0; next
	 * 								d(1) records are attributes of edges incident at 1 and so on.	 *
	 * 			.edge				+Text file. Each line contains two ends of an edge followed by their attributes.
	 * 								Names of vertices can be any strings that does not contain
	 * 								blank characters (space, tab, ,.etc).
	 * 								+ No comments are allowed
	 * 			.matrix				+Text file, adjacency matrix. The first line contains the
	 * 								number of vertices n. Next n lines contain on each lines
	 * 								n records that are edges' attributes.
	 * @param	weighted	The graph is weighted or not, default is weighted
	 */
	GeneralGraph(string fileName, bool weighted = false);

	/**
	 * @brief default constructor
	 */
	GeneralGraph() {
		vertices = arcs = 0;
		vDegrees = new int[1];
		neighbors = new int*[1];
		nbw = new T*[1];
		edgeList = new int[1];
		attr = NULL;
		inDegree = NULL;
		wDegree = NULL;
		wIndegree = NULL;
	}
	/**
	 * @brief Create a graph from a list of edges and weights. Call graphFromEdgeList
	 * @param vertices	Number of vertices in the graph
	 * @param arcs	Number of arcs/directed edges
	 * @param head	Array contains heads of edges
	 * @param tail	Array contains tails of edges
	 * @param attribute Weights of edges, if NULL then graph is unweighted
	 */
	GeneralGraph(int vertices, int arcs, const int* head, const int* tail,
			const T* attribute = NULL);

	/**
	 * @brief	Return number of vertices in the graph
	 */
	int Vertices() const {
		return vertices;
	}

	/**
	 * @brief	Return number of vertices in the graph
	 */
	int size() const {
		return vertices;
	}

	/**
	 * @brief	Return number of edges in the graph
	 */
	int Edges() const {
		return arcs / 2;
	}

	/**
	 * @brief	Return number of edges in the graph
	 */
	int Arcs() const {
		return arcs;
	}

	/**
	 * @brief	Copy constructor
	 */
	GeneralGraph(const GeneralGraph<T> &g);

	/**
	 * @brief	Assignment operator
	 */
	const GeneralGraph<T> &operator=(const GeneralGraph<T> &g);

	/**
	 * @brief	Return a pointer to the array that contains neighbors of a vertex 0 <= u <= n-1
	 * 			Return end of the edge list if u = n
	 * 			Return NULL otherwise
	 */
	int* operator[](int u);

	/**
	 * @brief	Return a pointer to the array that contains attributes of edges incident with u
	 * 			Return end of the attribute array if u = n
	 * 			Return NULL otherwise
	 */
	T* Weight(int u);

	/**
	 * @brief	Return a constant pointer to the array that contains neighbors of a vertex 0 <= u <= n-1
	 * 			Return end of the edge list if u = n
	 * 			Return NULL otherwise
	 */
	const int* operator[](int u) const;

	/**
	 * @brief	Return constant a pointer to the array that contains attributes of edges incident with u
	 * 			Return end of the attribute array if u = n
	 * 			Return NULL otherwise
	 */
	T* Weight(int u) const;

	/**
	 * @brief	Return reference to the weight of an edge with given id
	 * 			Return NULL otherwise
	 */
	T& edw(int id);

	/**
	 * @brief	Return reference to the weight of an edge with given id
	 * 			Return NULL otherwise
	 */
	const T& edw(int id) const;

	/**
	 * @brief	Return reference to the tail of the edge with given id
	 */
	int& tail(int id);

	/**
	 * @brief	Return reference to the tail of the edge with given id
	 */
	const int& tail(int id) const;

	/**
	 * @brief 	Return true if two graphs have same set of vertices and arcs and attributes
	 */
	bool operator==(const GeneralGraph<T> &g) const;


	/**
	 * Copy the topology of the graph, and set the weights to NULL
	 * @param g
	 */
	template <class T2>
	void topology_copy(const GeneralGraph<T2> &g) {
		// Release allocated memory
		free_array(vDegrees);
		free_array(neighbors);
		free_array(edgeList);
		free_array(attr);
		free_array(nbw);
		free_array(wDegree);
		free_array(wIndegree);
		free_array(inDegree);
		// Copy g's topology
		vid = g.vid;
		vertices = g.Vertices();
		arcs = g.Arcs();
		vDegrees = new int[vertices];
		neighbors = new int*[vertices + 1];
		edgeList = new int[arcs];
		directed = g.directed;
		MY_FOR(i, vertices)
			vDegrees[i] = g.degree(i);
		neighbors[0] = edgeList;
		MY_FOR(i, vertices)
			neighbors[i+1] = neighbors[i] + vDegrees[i];
		memcpy(edgeList, g[0], arcs * sizeof(int));
	}

	/**
	 * @brief  Return degree of a vertex u
	 */
	int degree(int u) const;

	/**
	 * @brief Write the graph to a file WITHOUT weights/attributes
	 * 		  Based on the fileName extension the corresponding format will be used.
	 * @see   GeneralGraph(const char *fileName)
	 */
	bool writeUnweighted(const char* fileName);

	/**
	 * @brief Write the graph to a file WITHOUT weights/attributes
	 * 		  Based on the fileName extension the corresponding format will be used.
	 * @see   GeneralGraph(const char *fileName)
	 */
	void readUnweighted(const char* fileName);

	/**
	 * @brief Write the graph to a file WITH weights/attributes
	 * Based on the fileName extension the corresponding format will be used.
	 * @see   GeneralGraph(const char *fileName)
	 */
	bool writeToFile(const char* fileName);

	/**
	 * @brief Remove duplicate edges/arcs and loops from the graph
	 * @param aggregate	If true, then the weights of duplicate edges
	 * will be combined (summed up).
	 * @param self_loop If true, then the self_loop will NOT be removed
	 * from the graph
	 */
	GeneralGraph<T> &simplify(bool aggregate=false, bool self_loop=false);

	/**
	 * Convert the graph into a undirected graph by doubling directed edges
	 * If loop_discarded is true, any existing loops will be discarded.
	 */
	void toUndirected(bool loop_discarded= true);

	/**
	 * @brief	Return an directed graph on the same set of vertices with all of
	 * the edges reversed compared to the orientation of the corresponding
	 * edges in the original graph. If the underlying graph is undirected,
	 * the function will not modify the graph.
	 */
	void get_inverse_graph(GeneralGraph<T> &gi) const {
			int *head = new int[arcs];
			ALL_EDGES(*this, h, t)
				head[id] = h;
			gi.graphFromEdgeList(vertices, arcs, edgeList, head, attr);
			delete[] head;
	}
	//
	//	/**
	//	 * @brief	Return a graph on the same vertices such that two vertices of the new
	//	 * graph are adjacent if and only if they are not adjacent in the
	//	 * original graph. Work both with directed & undirected graph.
	//	 *
	//	 */
	//	const GeneralGraph &complement();


	/**
	 *
	 * @return Whether the graph is weighted
	 */
	bool isWeighted() const {
		return (attr != NULL);
	}

	int in_degree(int u) {
		if (!inDegree) {
			inDegree = new int[vertices];
			for (int u = 0; u < vertices; ++u)
				inDegree[u] = 0;
			for (int i = 0; i < arcs; ++i)
				inDegree[tail[i]]++;
		}
		return inDegree[u];
	}

	inline T wdeg_out(int u) {
		if (!wDegree) {
			wDegree = new T[vertices];
		    MY_FOR(u,vertices)
			    wDegree[u] = 0;            
            ALL_EDGES( *this, h, t)                
			    wDegree[h] += attr[id];                                    
	    }
	    return wDegree[u];
    }

	T wdeg_in(int u) {
		if (!wIndegree) {
			wIndegree = new T[vertices];
			MY_FOR(u, vertices)
				wIndegree[u] = 0;
			ALL_EDGES(*this, h, t)
					wIndegree[t] += attr[id];            
            
		}
		return wIndegree[u];
	}
	/**
	 * @brief Assign all edges a same weight wg.
	 * The weight vector will be allocated if the graph is unweighted.
	 */
	void	setWeight(T w) {        
		if (!attr) {
			nbw		  = new T*[vertices + 1 ];
			attr	  = new T [arcs ];
			nbw[0]	  = attr;
			MY_FOR(u, vertices)
				nbw[u + 1] = nbw[u] + vDegrees[u];
		}
		MY_FOR(e, arcs)
			attr[e] = w;        
        free_array(wDegree);
        free_array(wIndegree);
	}

	/**
	 * @brief Assign a new weight vector for edges.
	 * The weight vector will be allocated if the graph is unweighted.
	 */
	void setWeight(T* vw) {
		if (!attr) {
			nbw = new T*[vertices + 1];
			attr = new T[arcs];
			nbw[0] = attr;
			MY_FOR(u, vertices)
				nbw[u + 1] = nbw[u] + vDegrees[u];
		}
		MY_FOR(e, arcs)
			attr[e] = vw[e];
        inDegree = NULL;
        wDegree = wIndegree= NULL;
	}

	/**
	 * @brief	Deconstructor
	 */
	virtual ~GeneralGraph() {
		free_array(vDegrees);
		free_array(neighbors);
		free_array(edgeList);
		free_array(nbw);
		free_array(attr);
		free_array(inDegree);
		free_array(wDegree);
		free_array(wIndegree);

	}
};


//-----------------------------------------------------------------------------
// MurmurHash2, by Austin Appleby

// Note - This code makes a few assumptions about how your machine behaves -

// 1. We can read a 4-byte value from any address without crashing
// 2. sizeof(int) == 4

// And it has a few limitations -

// 1. It will not work incrementally.
// 2. It will not produce the same results on little-endian and big-endian
//    machines.
inline unsigned int MurmurHash2( const void * key, int len, unsigned int seed )
{
	// 'm' and 'r' are mixing constants generated offline.
	// They're not really 'magic', they just happen to work well.

	const unsigned int m = 0x5bd1e995;
	const int r = 24;

	// Initialize the hash to a 'random' value

	unsigned int h = seed ^ len;

	// Mix 4 bytes at a time into the hash

	const unsigned char * data = (const unsigned char *)key;

	while(len >= 4)
	{
		unsigned int k = *(unsigned int *)data;

		k *= m;
		k ^= k >> r;
		k *= m;

		h *= m;
		h ^= k;

		data += 4;
		len -= 4;
	}

	// Handle the last few bytes of the input array

	switch(len)
	{
	case 3: h ^= data[2] << 16;
	case 2: h ^= data[1] << 8;
	case 1: h ^= data[0];
	        h *= m;
	};

	// Do a few final mixes of the hash to ensure the last few
	// bytes are well-incorporated.

	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;

	return h;
}


struct hash_string {
	std::size_t operator() (string s) const {
		return MurmurHash2((const void *)s.c_str(), s.length(), 0);
	}
};

#ifdef	_MSC_VER
	#include <unordered_map>
    typedef std::tr1::unordered_map<string,int, hash_string > nameIdTable;
    typedef std::tr1::unordered_map<int,int> 	idLookupTable;
#else
	#include <tr1/unordered_map>
	typedef std::tr1::unordered_map<string,int, hash_string > nameIdTable;
	typedef std::tr1::unordered_map<int,int> 	idLookupTable;
#endif



using namespace std::tr1;




inline string	getFileExtension(string fileName) {
	return fileName.substr( fileName.find_last_of(".")+1);
}

inline string	getFileName(string fileName) {
	return fileName.substr(0, fileName.find_last_of("."));
}

template <class T>
void GeneralGraph<T>::graphFromEdgeList(int vertices, int arcs, const int* head,  const int* tail, const T* attribute) {
	this->vertices = vertices;
	this->arcs = arcs;
	vDegrees  = new int[ vertices];
	neighbors = new int*[vertices + 1 ];
	edgeList  = new int[ arcs ];
	inDegree = NULL;
	wDegree = wIndegree = NULL;
    directed = true;
	if (attribute != NULL) {
		nbw		  = new T*[vertices + 1 ];
		attr	  = new T[arcs ];
	} else {
		attr = NULL;
		nbw = NULL;
	}
	// Count degree
	fill( vDegrees, vDegrees + vertices, 0);
	MY_FOR(i, arcs) {
		CHECK(head[i] >=0 && head[i] < vertices && tail[i] >= 0
				&& tail[i] < vertices,
				"Error:generalGraph<T>::graphFromEdgeList: Invalid edges ("
				<<head[i]<<";"<<tail[i]<<") ")
		++vDegrees[ head[i] ];			// Assume the graph is directed
		//++vDegrees[ tail[i] ];
	}
	// Find position to put neighbors
	neighbors[0] = edgeList;
	MY_FOR(i, vertices)
		neighbors[i+ 1]= neighbors[i] + vDegrees[i];
	if (attr) {
		nbw[0] = attr;
		MY_FOR(i, vertices)
			nbw[i+1] = nbw[i] + vDegrees[i];
	}
	MY_FOR(i, arcs) {
		*(neighbors[ head[i] ]) = tail[i];
		++neighbors[ head[i] ];
		if (attr) {
			*(nbw[ head[i] ]) = attribute[i];
			++nbw[ head[i] ];
		}
	}
	MY_FOR(i, vertices) {
		neighbors[i] = neighbors[i] - vDegrees[i];
		if (attr)
			nbw[i] = nbw[i] - vDegrees[i];
	}
}

/**
 * @brief	Print out the graph in human readable form
 * 			Start with number of vertices, edges on the 1st line
 * 			On next lines, each line starts with a vertex id following
 * 			with its neighbors
 */
template<class T>
ostream& operator<<(ostream &os, const GeneralGraph<T> &g) {
	os << g.Vertices() << "\t" << g.Arcs() << endl;
	MY_FOR(v, g.Vertices())
		if (g.degree(v) > 0) {
			os << v;
			MY_FOR(i, g.degree(v)) {
				os << "\t" << g[v][i];
				if (g.isWeighted())
					os << " " << g.Weight(v)[i];
			}
			os << endl;
		}
	return os;
}

template <class T>
GeneralGraph<T> &GeneralGraph<T>::simplify(bool aggregate, bool self_loop) {
	bool weighted = (attr != NULL);
	// Remove duplicate edges
	int *tmpList = edgeList;
	int **tmpNb =  neighbors;
	T *tmpAttr	= attr;
	T **tmpNbw	= nbw;
	int v;
	int *mask = new int[vertices];
	vector<int> wl(vertices, -1);
	fill(mask, mask + vertices, -1);
	MY_FOR(u, vertices) {
		int du = 0;
		if (!self_loop)
			mask[u] = u; // Forbid self-loop
		MY_FOR(i, vDegrees[u]) {
			v = tmpNb[u][i];
			if (mask[v] != u) {
				mask[v] = u;
				wl[v] = du;
				tmpNb[u][du] = v;
				if (weighted)
					tmpNbw[u][du] = tmpNbw[u][i];
				++du;
			} else if (aggregate &&( (v != u)||(self_loop) ) )
				tmpNbw[u][ wl[v] ] += tmpNbw[u][i];
		}
		vDegrees[u] = du;
	}
	delete[] mask;
	// Copy data back
	//arcs = accumulate(vDegrees, vDegrees + vertices, 0);
	arcs = 0;
	MY_FOR(u, vertices)
		arcs += vDegrees[u];
	edgeList  = new int[arcs];
	neighbors = new int*[ vertices + 1];
	neighbors[0] = edgeList;
	MY_FOR(u, vertices) {
		neighbors[u + 1] = neighbors[u] + vDegrees[u];
		memcpy((void*) neighbors[u], (void*) tmpNb[u], vDegrees[u]
				* sizeof(int));
	}
	delete[] tmpNb;
	delete[] tmpList;
	if (weighted) {
		attr = new T[arcs];
		nbw  = new T*[ vertices + 1];
		nbw[0] = attr;
		MY_FOR(u, vertices) {
			nbw[u + 1] = nbw[u] + vDegrees[u];
			MY_FOR(i, vDegrees[u])
				nbw[u][i] = tmpNbw[u][i];
			//memcpy((void*) nbw[u], (void*) tmpNbw[u], vDegrees[u] * sizeof(T));
		}
		delete[] tmpNbw;
		delete[] tmpAttr;
	}
	return *this;
}

template <class T>
void GeneralGraph<T>::readUnweighted(const char* fileName) {
	attr = NULL;// Unweighted
	string ext = getFileExtension(fileName);
	FILE *ft = fopen(fileName, "r");
	CHECK(ft, "Error: generalGraph("<<fileName<<"): Cannot access file "<<fileName<<endl);
	fclose(ft);
    directed = true;
	char buffer[1001];
	if (ext == "adj") {
		FILE* f = fopen(fileName, "r");
		int vertices = -1;
		// Get the number of edges and vertices
		while (!feof(f)) {
			fgets(buffer, 1000, f);
			if (buffer[0] == '#')
				continue;
			istringstream(buffer) >> vertices >> arcs;
			break;
        }
		CHECK(vertices > 0, "Error: generalGraph("<<fileName<<"): invalid file format!"<<endl )
		// Read edges and put into a vector
		vector<int> head, tail;
		head.reserve(arcs);
		tail.reserve(arcs);
		int u, v;
        unsigned int one_flag = 0;
		while (!feof(f)) {
			// Allow comments to interlace with data in the DEBUG mode
            u = -1; v = -1; buffer[0]=0;
#ifdef	_DEBUG
			fgets(buffer, 1000, f);
			if (buffer[0] == '#')
				continue;
			istringstream(buffer) >> u >> v;
#else
			fscanf(f,"%d %d",&u, &v);
#endif
            if (u == -1 && v == -1)
                break;
			head.push_back(u);
			tail.push_back(v);
            if (u == 0 || v == 0) {
                one_flag = one_flag | 2;
            }
            if (one_flag == 0 && (u == vertices || v == vertices) ) {
                cerr <<"Warning: Vertices are NOT ZERO-based. Automatically subtract all vertices indices by 1!"<<endl;
                one_flag = one_flag |  1;
            }
            if (u < 0 || v < 0 || u > vertices || v > vertices) {
                cerr <<"Error: Invalid vertices found: "<< u <<" - "<< v <<"!"<<endl;
                exit(3);
            }
		}
        if (arcs != (int)head.size() ) {
            cout <<"Warning: Number of reported edges/arcs mismatched. Reported "<< arcs
                 <<", actually found "<< head.size() << endl;
        } 
        if (one_flag == 3) {
            cerr <<"Error: Vertex index out of range!"<<endl;
            exit(1);
        }
        if (one_flag == 1) {    // Auto subtract 1
            MY_FOR(i, head.size() ) {
                --head[i];
                --tail[i];
            }
        }
		graphFromEdgeList(vertices, head.size(), &head[0], &tail[0]);
		fclose(f);
	} else if (ext == "bin") {
		FILE* f = fopen(fileName, "rb");
		CHECK(f, "Error: Cannot read from file"<< fileName )
		fread((void *) &vertices, sizeof(int), 1, f);
		// Get degree of each vertex
		vDegrees = new int[vertices];
		fread((void *) vDegrees, sizeof(int), vertices, f);
		arcs = accumulate(vDegrees, vDegrees + vertices, 0);
		// Get all the edges at once
		edgeList = new int[arcs];
		fread((void *) edgeList, sizeof(int), arcs, f);
		fclose(f);
		// Find the starting points of neigbors for each vertex
		neighbors = new int*[vertices + 1];
		neighbors[0] = edgeList;
		MY_FOR(i, vertices)
			neighbors[i + 1] = neighbors[i] + vDegrees[i];
	} else if (ext == "edge") {
		FILE* f = fopen(fileName, "r");
		nameIdTable vertexId;
		vector<int> head, tail;
		int vid = 0;
		vector<string> names;
		while (!feof(f)) {
			if (fgets(buffer, 1000, f)) {
				if (buffer[0] == '#') // Ignore comments
					continue;
				char* b1 = buffer;
				// Get the first vertex
				while (*b1 == ' ' || *b1 == '\t' || *b1 == '\r' || *b1 == '\n')
					b1++;
				char* e1 = b1;
				while (*e1 != 0 && *e1 != ' ' && *e1 != '\t' && *e1 != '\r'
						&& *e1 != '\n')
					e1++;
				if (e1 == b1)
					continue;
				// Get the second vertex
				char* b2 = e1;
				while (*b2 == ' ' || *b2 == '\t' || *b2 == '\r' || *b2 == '\n')
					b2++;
				char* e2 = b2;
				while (*e2 != 0 && *e2 != ' ' && *e2 != '\t' && *e2 != '\r'
						&& *e2 != '\n')
					e2++;
				if (e2 == b2)
					continue;
				*e1 = *e2 = 0;
				// Add new vertex id if necessary
				if (vertexId.find(b1) == vertexId.end()) {
					vertexId[b1] = vid;
					names.push_back(b1);
					vid++;
				}
				if (vertexId.find(b2) == vertexId.end()) {
					vertexId[b2] = vid;
					names.push_back(b2);
					vid++;
				}
				// Add edges to the edge list
				head.push_back(vertexId[b1]);
				tail.push_back(vertexId[b2]);
			} // if (fgets
		} // while
		fclose(f);
		graphFromEdgeList(vid, head.size(), &head[0], &tail[0]);
		// Write out the map file
		FILE *fm = fopen((string(fileName) + ".map").c_str(), "wt");
		if (fm) {
			MY_FOR(i, vid)
				fprintf(fm, "%d\t%s\n", i, names[i].c_str());
			fclose(fm);
		}
	} else if (ext == "matrix") {
		FILE* f = fopen(fileName, "r");
		// Get the number of vertices
		int vertices = -1;
		// Get the number of edges and vertices
		while (!feof(f)) {
			fgets(buffer, 1000, f);
			if (buffer[0] == '#')
				continue;
			istringstream(buffer) >> vertices;
			break;
		}
		CHECK(vertices > 0, "Error: generalGraph("<<fileName<<"): invalid file format!"<<endl )
		vector<int> head, tail;
		head.reserve(vertices);
		tail.reserve(vertices);
		int x;
		MY_FOR(i, vertices)
			MY_FOR(j, vertices ) {
				fscanf(f, "%d", &x);
				if (x) {
					head.push_back(i);
					tail.push_back(j);
				}
			}
		fclose(f);
		graphFromEdgeList(vertices, head.size(), &head[0], &tail[0]);
	} else {
		CHECK(0, "Error: generalGraph(<<"<<fileName<<"): Cannot recognize file name extension!"<<endl)
	}
}


template<class T>
GeneralGraph<T>::GeneralGraph(string fileName, bool weighted):inDegree(NULL),wDegree(NULL),wIndegree(NULL),attr(NULL) {
	if (!weighted) {
		attr = NULL;
		readUnweighted(fileName.c_str());
		return;
	}
    directed = true;
	string ext = getFileExtension(fileName);
	FILE *ft = fopen(fileName.c_str(), "r");
	CHECK(ft, "Error: generalGraph("<<fileName<<"): Cannot access file "<<fileName<<endl);
	fclose(ft);
	if (ext == "adj") {
		ifstream f(fileName.c_str());
		int vertices = -1;
		// Get the number of edges and vertices
		while (!f.eof() && (f.good())) {
			if (f.peek() == '#')
#undef max
				f.ignore(numeric_limits<int>::max(), '\n');
			else {
				f >> vertices >> arcs;
				f.ignore(numeric_limits<int>::max(), '\n');
				break;
			}
		}
		CHECK(vertices > 0, "Error: generalGraph("<<fileName<<"): invalid file format!"<<endl )
		vector<int> head, tail;
		vector<T> weight;
		head.reserve(arcs);
		tail.reserve(arcs);
		weight.reserve(arcs);
		int u, v;
		bool zero_flag = false, v_flag = false;
		T w;
		while (!f.eof()) {
			// Allow comments to interlace with data in the DEBUG mode
			u = -1; v = -1;
			if (f.peek() == '#') {
				f.ignore(numeric_limits<int>::max(), '\n');
				continue;
			}
			f >> u >> v >> w;
			f.ignore(numeric_limits<int>::max(), '\n');
			if (u == -1 || v == -1)
				break;
			if (u < 0 || v < 0 || u > vertices || v > vertices) {
			    cerr <<"Error: Invalid vertices found: "<< u <<" - "<< v <<"!"<<endl;
			    exit(3);
			}
			head.push_back(u);
			tail.push_back(v);
			if (u == 0 || v == 0)
				zero_flag = true;
			if (u == vertices || v == vertices)
				v_flag = true;
			weight.push_back(w);
		}
		if ((!zero_flag) && (v_flag)) {
			cerr <<"Warning: Vertices are NOT ZERO-based. "
				 <<"Automatically subtract all vertices indices by 1!"<< endl;
			MY_FOR(i, head.size() ) {
				--head[i];
				--tail[i];
			}
		}
		if (arcs != (int) head.size()) {
			cout << "Warning: Number of reported "
				 << "edges/arcs mismatched. Reported "
				 << arcs << ", actually found " << head.size() << endl;
		}
		graphFromEdgeList(vertices, head.size(), &head[0], &tail[0], &weight[0]);
		f.close();
	} else if (ext == "bin") {
		FILE* f = fopen(fileName.c_str(), "rb");
		fread((void *) &vertices, sizeof(int), 1, f);
		// Get degree of each vertex
		vDegrees = new int[vertices];
		fread((void *) vDegrees, sizeof(int), vertices, f);
		arcs = accumulate(vDegrees, vDegrees + vertices, 0);
		// Get all the edges at once
		edgeList = new int[arcs];
		fread((void *) edgeList, sizeof(int), arcs, f);
		// Get all the weights/attribute
		attr = new T[arcs];
		fread((void *) attr, sizeof(T), arcs, f);
		// Find the starting points of neighbors for each vertex
		neighbors = new int*[vertices + 1];
		nbw = new T*[vertices + 1];
		neighbors[0] = edgeList;
		nbw[0] = attr;
		MY_FOR(i, vertices) {
			neighbors[i + 1] = neighbors[i] + vDegrees[i];
			nbw[i + 1] = nbw[i] + vDegrees[i];
		}
		fclose(f);
	} else if (ext == "edge") {
		ifstream f(fileName.c_str());
		nameIdTable vertexId;
		vector<int> head, tail;
		vector<T> weight;
		int vid = 0;
		vector<string> names;
		string s1, s2;
		T w;
		while (!f.eof()) {
			s1 = s2 = "";
			f >> s1;
			if (f.eof() ) break;
#undef max
			if (s1.length() < 1 || s1[0] == '#' ) {
				f.ignore(numeric_limits<int>::max(), '\n');
				continue;
			}
			f >> s2 >> w;
			if (f.eof() ) break;
			if (s2.length() < 1) {
				f.ignore(numeric_limits<int>::max(), '\n');
				continue;
			}
			//cout <<"@:"<<s1<<"-"<<s2<<endl;
			// Trim spaces
			//			int st = s1.find_first_not_of(" \t\n");
			//			int len = s1.find_last_not_of(" \t\n") - st + 1;
			//			s1 = s1.substr(st, len);
			//			st = s2.find_first_not_of(" \t\n");
			//			len = s2.find_last_not_of(" \t\n") - st + 1;
			//			s2 = s2.substr(st, len);
			// Add new vertex id if necessary
			if (vertexId.find(s1) == vertexId.end()) {
				vertexId[s1] = vid;
				names.push_back(s1);
				vid++;
			}
			if (vertexId.find(s2) == vertexId.end()) {
				vertexId[s2] = vid;
				names.push_back(s2);
				vid++;
			}
			head.push_back(vertexId[s1]);
			tail.push_back(vertexId[s2]);
			weight.push_back(w);
		} // while
		f.close();
		graphFromEdgeList(vid, head.size(), &head[0], &tail[0], &weight[0]);
		// Write out the map file
		FILE *fm = fopen((string(fileName) + ".map").c_str(), "wt");
		if (fm) {
			MY_FOR(i, vid)
				fprintf(fm, "%d\t%s\n", i, names[i].c_str());
			fclose(fm);
		}
	} else if (ext == "matrix") {
		ifstream f(fileName.c_str());
		static T wz; // Zero element
		// Get the number of vertices
		int vertices = -1;
		// Get the number of edges and vertices
		while (!f.eof() && f.peek()=='#')
			f.ignore(numeric_limits<int>::max(), '\n');
		CHECK(vertices > 0, "Error: generalGraph("<<fileName<<"): invalid file format!"<<endl )
		// vector< pair<int, int> > edgeList;
		vector<int> head, tail;
        vector<T> weight;

		head.reserve(vertices);
		tail.reserve(vertices);
		T x;
		MY_FOR(i, vertices)
			MY_FOR(j, vertices ) {
				f >> x;
				if (!(x == wz)) {
					head.push_back(i);
					tail.push_back(j);
					weight.push_back(x);
				}
			}
		f.close();
		graphFromEdgeList(vertices, head.size(), &head[0], &tail[0], &weight[0]);
	} else {
		CHECK(0, "Error: generalGraph(<<"<<fileName<<"): Cannot recognize file name extension!"<<endl)
	}
}

template <class T>
GeneralGraph<T>::GeneralGraph(int vertices, int arcs, const int* head, const int* tail,
		const T* attribute) {
	graphFromEdgeList(vertices, arcs, head, tail, attribute);
}







/**
 * @brief	Copy constructor
 */
template <class T>
GeneralGraph<T>::GeneralGraph(const GeneralGraph<T> &g):
inDegree(NULL),wDegree(NULL),wIndegree(NULL),attr(NULL), neighbors(NULL), vDegrees(NULL),edgeList(NULL),nbw(NULL)
{
	duplicate(g);
}


/**
 * @brief	Assignment operator
 */
template <class T>
const GeneralGraph<T> &GeneralGraph<T>::operator=(const GeneralGraph<T> &g) {
	if (this != &g)
		duplicate(g);
	return *this;
}

/**
 * @brief	Return a pointer to the array that contains neighbors of a vertex 0 <= u <= n-1
 * 			Return end of the edge list if u = n
 * 			Return NULL otherwise
 */
template <class T>
inline int* GeneralGraph<T>::operator[](int u) {
	DCHECK( (0<=u) && (u<=vertices), "generalGraph<T>::operator[]: index out of range "<<u<<endl);
	return neighbors[u];
}

/**
 * @brief	Return  a pointer to the array that contains attributes of edges incident with u
 * 			Return end of the attribute array if u = n
 * 			Return NULL otherwise
 */
template<class T>
inline T* GeneralGraph<T>::Weight(int u) {
	DCHECK( (0<=u) && (u<=vertices) && (attr!=NULL), "generalGraph<T>::Weight: index out of range or unweighted graph!"<<u<<endl);
	return (nbw[u]);
}

/**
 * @brief	Return constant a pointer to the array that contains attributes of edges incident with u
 * 			Return end of the attribute array if u = n
 * 			Return NULL otherwise
 */
template<class T>
inline T* GeneralGraph<T>::Weight(int u) const {
	DCHECK( (0<=u) && (u<=vertices) && (attr != NULL), "generalGraph<T>::Weight (const): index out of range or unweighted graph!"<<u<<endl);
	return (nbw[u]);
}


/**
 * @brief	Return a constant pointer to the array that contains neighbors of a vertex 0 <= u <= n-1
 * 			Return end of the edge list if u = n
 * 			Return NULL otherwise
 */
template <class T>
inline const int* GeneralGraph<T>::operator[](int u) const {
	DCHECK( (0<=u) && (u<=vertices), "const generalGraph<T>::operator[]: index out of range "<<u<<endl);
	return neighbors[u];
}


/**
 * @brief  Return degree of a vertex u
 */
template <class T>
inline int GeneralGraph<T>::degree(int u) const {
	DCHECK( (0<=u) && (u<=vertices), "generalGraph<T>::degree: index out of range "<<u<<endl);
	return vDegrees[u];
}

template <class T>
bool GeneralGraph<T>::writeUnweighted(const char* fileName) {
	string ext = getFileExtension(fileName);
	if (ext == "adj") {
		FILE* f = fopen(fileName, "wt");
		CHECK( f,"Error:generalGraph<T>::writeToFile: Cannot access file " << fileName << endl)
		fprintf(f, "%d %d\n", vertices, arcs);
		MY_FOR(u, vertices)
			for (int *v = neighbors[u]; v != neighbors[u + 1]; ++v)
				fprintf(f, "%d %d\n", u, *v);
		fclose(f);
	} else if (ext == "bin") {
		FILE* f = fopen(fileName, "wb");
		CHECK( f,"Error:generalGraph<T>::writeUnweighted: Cannot access file " << fileName << endl)
		fwrite((void *) &vertices, sizeof(int), 1, f);
		fwrite((void *) vDegrees, sizeof(int), vertices, f);
		fwrite((void *) edgeList, sizeof(int), arcs, f);
		fclose(f);
	} else if (ext == "edge") {
		FILE* f = fopen(fileName, "wt");
		CHECK( f,"Error:generalGraph<T>::writeToFile: Cannot access file " << fileName << endl)
		fprintf(f, "# %d vertices, %d arcs \n", vertices, arcs);
		MY_FOR(u, vertices)
			for (int *v = neighbors[u]; v != neighbors[u + 1]; ++v)
				fprintf(f, "%d %d\n", u, *v);
		fclose(f);
	} else if (ext == "matrix") {
		FILE* f = fopen(fileName, "wt");
		CHECK( f,"Error:generalGraph<T>::writeToFile: Cannot access file " << fileName << endl)
		fprintf(f, "%d\n", vertices);
		MY_FOR(u, vertices) {
			vector<int> mask(vertices, 0);
			for (int *v = neighbors[u]; v != neighbors[u + 1]; ++v)
				mask[*v] = 1;
			MY_FOR(i, vertices - 1)
				fprintf(f, "%d ", mask[i]);
			fprintf(f, "%d\n", mask[vertices - 1]);
		}
		fclose(f);
	} else {
		CHECK(0,"Error:generalGraph<T>::writeToFile("<<fileName<<"): Cannot recognize file extension."<< endl)
	}
	return true;
}

/**
 * @brief Read the graph from a file
 * Based on the fileName extension the corresponding format will be used.
 * @see   GeneralGraph(const char *fileName)
 */
template<class T>
bool GeneralGraph<T>::writeToFile(const char* fileName) {
	if (attr == NULL)
		return writeUnweighted(fileName);
	string ext = getFileExtension(fileName);
	int*v;
	T *w;
	if (ext == "adj") {
		ofstream f(fileName);
		CHECK( f.good(),"Error:generalGraph<T>::writeToFile: Cannot access file " << fileName << endl)
		f << vertices << " " << arcs << endl;
		MY_FOR(u, vertices)
			for (v = neighbors[u], w = nbw[u]; v != neighbors[u + 1]; ++v, ++w)
				f << u << "\t" << *v << "\t" << *w << endl;
		f.close();
	} else if (ext == "bin") {
		FILE* f = fopen(fileName, "wb");
		CHECK( f,"Error:generalGraph<T>::writeToFile: Cannot access file " << fileName << endl)
		fwrite((void *) &vertices, sizeof(int), 1, f);
		fwrite((void *) vDegrees, sizeof(int), vertices, f);
		fwrite((void *) edgeList, sizeof(int), arcs, f);
		fwrite((void *) attr, sizeof(T), arcs, f);
		fclose(f);
	} else if (ext == "edge") {
		ofstream f(fileName);
		CHECK( f.good(),"Error:generalGraph<T>::writeToFile: Cannot access file " << fileName << endl)
		f << "# " << vertices << " vertices, " << arcs << " arcs " << endl;
		MY_FOR(u, vertices)
			for (v = neighbors[u], w = nbw[u]; v != neighbors[u + 1]; ++v, ++w)
				f << u << "\t" << *v << "\t" << *w << endl;
		f.close();
	} else if (ext == "matrix") {
		ofstream f(fileName);
		CHECK( f.good(),"Error:generalGraph<T>::writeToFile: Cannot access file " << fileName << endl)
		f << vertices << endl;
		MY_FOR(u, vertices) {
			vector<T> mask(vertices, T(0)); // default weights
			for (v = neighbors[u], w = nbw[u]; v != neighbors[u + 1]; ++v, ++w)
				mask[*v] = *w;
			MY_FOR(i, vertices - 1)
				f << mask[i] << endl;
			f << mask[vertices - 1] << endl;
		}
		f.close();
	} else {
		CHECK(0,"Error:generalGraph<T>::writeToFile("<<fileName<<"): Cannot recognize file extension."<< endl)
	}
	return true;
}

/**
 * Convert the graph into a general graph by doubling directed edges
 * Any existing loops will be discarded.
 */
template<class T>
void GeneralGraph<T>::toUndirected(bool loop_discarded) {
	bool weighted = (attr != NULL);
    directed = false;
	// Doubling every edges (u,v) into (u, v) and (v, u)
	MY_FOR(i, arcs)
		++vDegrees[edgeList[i]];
	int *tmpList = new int[arcs * 2];
	int **tmpNb = new int*[vertices + 1];
	T *tmpAttr;
	T **tmpNbw;
	tmpNb[0] = tmpList;
	MY_FOR(u, vertices)
		tmpNb[u + 1] = tmpNb[u] + vDegrees[u];
	MY_FOR(u, vertices) {
		memcpy((void*) tmpNb[u], (void*) neighbors[u], sizeof(int)
				* (neighbors[u + 1] - neighbors[u]));
		tmpNb[u] += neighbors[u + 1] - neighbors[u];
	}
	if (weighted) {
		tmpAttr = new T[arcs * 2];
		tmpNbw = new T*[vertices + 1];
		tmpNbw[0] = tmpAttr;
		MY_FOR(u, vertices)
			tmpNbw[u + 1] = tmpNbw[u] + vDegrees[u];
		MY_FOR(u, vertices) {
			for(int i=nbw[u + 1]- nbw[u]-1;i>=0;--i)
					tmpNbw[u][i] = nbw[u][i];
//			memcpy((void*) tmpNbw[u], (void*) nbw[u], sizeof(T) * (nbw[u + 1]
//					- nbw[u]));
			tmpNbw[u] += nbw[u + 1] - nbw[u];
		}
	}
	MY_FOR(u, vertices) {
		for (int *v = neighbors[u]; v != neighbors[u + 1]; ++v) {
			*tmpNb[*v] = u;
			tmpNb[*v]++;
			if (weighted) {
				*tmpNbw[*v] = *(attr + (v - edgeList));
				tmpNbw[*v]++;
			}
		}
	}
	delete[] edgeList;
	if (weighted)
		delete[] attr;
	MY_FOR(u, vertices) {
		tmpNb[u] -= vDegrees[u];
		if (weighted)
			tmpNbw[u] -= vDegrees[u];
	}
	// Remove duplicate edges
	int v;
	int *mask = new int[vertices];
	fill(mask, mask + vertices, -1);
	MY_FOR(u, vertices) {
		int du = 0;
        if (loop_discarded)
		    mask[u] = u; // Forbid self-loop
		MY_FOR(i, vDegrees[u]) {
			v = tmpNb[u][i];
			if (mask[v] != u) {
				mask[v] = u;
				tmpNb[u][du] = v;
				if (weighted)
					tmpNbw[u][du] = tmpNbw[u][i];
				++du;
			}
		}
		vDegrees[u] = du;
	}
	delete[] mask;
	// Copy data back
	//arcs = accumulate(vDegrees, vDegrees + vertices, 0);
	arcs = 0;
	MY_FOR(u, vertices)
		arcs += vDegrees[u];
	edgeList = new int[arcs];
	neighbors[0] = edgeList;
	MY_FOR(u, vertices) {
		neighbors[u + 1] = neighbors[u] + vDegrees[u];
		memcpy((void*) neighbors[u], (void*) tmpNb[u], vDegrees[u]
				* sizeof(int));
	}
	delete[] tmpNb;
	delete[] tmpList;
	if (weighted) {
		attr = new T[arcs];
		nbw[0] = attr;
		MY_FOR(u, vertices) {
			nbw[u + 1] = nbw[u] + vDegrees[u];
			MY_FOR(i, vDegrees[u])
				nbw[u][i] = tmpNbw[u][i];
			// memcpy((void*) nbw[u], (void*) tmpNbw[u], vDegrees[u] * sizeof(T));
		}
		delete[] tmpNbw;
		delete[] tmpAttr;
	}
    // 
    free_array(inDegree);
    free_array(wDegree);
    free_array(wIndegree);
}


/**
 * @brief 	Return true if two graphs have same set of vertices and arcs
 */
template <class T>
bool GeneralGraph<T>::operator==(const GeneralGraph<T> &g) const {
	if (this == &g)
			return true;
	if ( (vertices != g.vertices) || (arcs != g.arcs))
		return false;
	vector<int>	marker(vertices, 0);
	MY_FOR(u, vertices) {
		if (vDegrees[u] != g.vDegrees[u])
			return false;

		MY_FOR(i, vDegrees[u]) {
			 ++marker[ neighbors[u][i] ];
			 --marker[ g.neighbors[u][i] ];
		}
		MY_FOR(i, vDegrees[u])
			if (marker[ neighbors[u][i]])	// != 0
				return false;
	}
	return true;
}

/**
 * @brief	Return reference to the weight of an edge with given id
 * 			Return NULL otherwise
 */
template <class T>
inline T& 	GeneralGraph<T>::edw(int id) {
	DCHECK( (0<=id) && (id<arcs) && (attr != NULL), "generalGraph<T>::edw: index out of range or unweighted graph!"<<id<<endl);
	return attr[id];
}

/**
 * @brief	Return const reference to the weight of an edge with given id
 * 			Return NULL otherwise
 */
template <class T>
inline const T& 	GeneralGraph<T>::edw(int id) const {
	DCHECK( (0<=id) && (id<arcs) && (attr != NULL), "generalGraph<T>::edw: index out of range or unweighted graph!"<<id<<endl);
	return attr[id];
}


/**
 * @brief	Return reference to the tail of the edge with given id
 */
template <class T>
inline int& GeneralGraph<T>::tail(int id) {
	DCHECK( (0<=id) && (id<arcs),"generalGraph<T>::edw: index out of range or unweighted graph!"<<id<<endl);
	return edgeList[id];
}

/**
 * @brief	Return const reference to the tail of the edge with given id
 */
template <class T>
inline const int& GeneralGraph<T>::tail(int id) const {
	DCHECK( (0<=id) && (id<arcs),"generalGraph<T>::edw: index out of range or unweighted graph!"<<id<<endl);
	return edgeList[id];
}



template <class T>
void GeneralGraph<T>::duplicate(const GeneralGraph<T> &g) {
	topology_copy<T>( g);
	// Copy weights if necessary
	if (g.attr != NULL) {		// Weighted graph
		nbw = new T*[vertices + 1];
		attr = new T[arcs];
		nbw[0] = attr;
		MY_FOR(i, vertices)
			nbw[i + 1] = nbw[i] + vDegrees[i];
		// memcpy(attr, g.attr, (arcs) * sizeof(T));
		MY_FOR(i, arcs)
			attr[i] = g.attr[i];
	} else attr = NULL;
}



typedef	 GeneralGraph<>	 Graph;



#endif /* GENERAL_GRAPH_H_ */

