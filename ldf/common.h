/*
 * common.h
 *
 *  Created on: Dec 24, 2008
 *      Author: Thang
 */

#ifndef COMMON_H_
#define COMMON_H_

#include <string>
#include <sstream>

//#define		IGRAPH	// Enable igraph
//#define		CPPUNIT	// Enable CppUnit Test
//#define		_DEBUG

// Quick macros
#define	byte	unsigned char
#define	FOR(i,n) for(int i =0;i !=(int)n;++i)
#define	FORAB(i, a, b) for(int i = (int)(a); i != (int)(b);i++)
//#define MIN(a,b) ((a)<(b)?(a):(b))
//#define MAX(a,b) ((a)>(b)?(a):(b))
#define	ALL(v)	 v.begin(), v.end()
#define	FOREACH(itr,v)	for( typeof(v.begin()) itr=v.begin(); itr!=v.end();itr++ )
#define	CHECK( a, b)  if (!(a) ) { cerr << b << endl; exit(1); }

#ifdef	_DEBUG
	#define		DCHECK( a, b)  if (!(a) ) { cerr << b << endl; exit(1); }
#else
	#define		DCHECK( a, b)  {}
#endif

// Quick utilities
template <class T>
inline std::string toString(T t) {
  std::ostringstream oss;
  oss << t;
  return oss.str();
}

// Free vector's memory
template <class T>
void	freeVector(T &v) {
	T tmp;
	v.swap(tmp);
}


// CPPUNIT header files
#ifdef	CPPUNIT
// Test Unit
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>

#define ASSERT 	CPPUNIT_ASSERT

#endif

// Igraph short aliases
#ifdef	IGRAPH
extern "C" {
#include <igraph.h>
}
// Short Alias
#define	Graph	igraph_t
#define	Vector	igraph_vector_t
#define	Integer	igraph_integer_t
#define Bool	igraph_bool_t
#define	Real	igraph_real_t
#define	Matrix	igraph_matrix_t
#define Selector igraph_vs_t
#define	VSIZE(a) igraph_vector_size(& (a))	// Size of a vector
#define	GSIZE(g) igraph_vcount(& g) 		// Number of vertices in a graph

#endif

#endif /* COMMON_H_ */
