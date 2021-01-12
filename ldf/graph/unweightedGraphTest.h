#include "edgeMap.hpp"
#include "generalGraph.hpp"
#include "graphAlgorithm.hpp"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

struct	TestUnweightedGraph {
	TestUnweightedGraph() {
//		unweightedGraph g("data/a.edge");
//		g.toUndirected();
//		DCHECK(g.Edges() == 4, "Error: Wrong number of edges! Expected 4, actual: " << g.Edges() )
//		DCHECK(g.Vertices() == 4, "Error: Wrong number of vertices! Expected 4, actual: "<<g.Vertices() )
//		g.writeToFile("data/a.bin");
//		unweightedGraph g2("data/a.bin");
//		DCHECK(g2==g,"Error: Invalid graph writing/reading!")
//		g2.writeToFile("data/a2.adj");
//		cout.flush();
////		unweightedGraph fb("D:/Research/thang/Data/facebook/facebooks.edge");
////		fb.toUndirected();
////		EdgeMap fbm(fb);
////		fb.writeToFile("D:/Research/thang/Data/facebook/facebooks.bin");
//		// Check simplify function
//		int	head[]={1, 1, 3, 1, 2};
//		int	tail[]={2, 3, 1, 2, 1};
//		GeneralGraph<> gg(4, 5, head, tail);
//		DCHECK( gg.Arcs() == 5 && gg.Vertices()==4 && gg.degree(1)==3
//				,"Invalid graph initialization");
//		gg.simplify();
//		DCHECK( gg.Arcs() == 4 && gg.Vertices()==4 && gg.degree(1)==2
//				,"Invalid graph simplification");
//
////		g = fb;
////		EdgeMap  gm(g);
////		const int *u = gm.getHead(0), *v = gm.getTail(0);
////		FOR(i, gm.size()) {
////			DCHECK(gm.edgeId(*u, *v) == i,"Error: Incorrect edgeId returned! ("
////					<<*u<<";"<<*v<<")->"<<i)
////			u++,v++;
////		}
////		g.toUndirected();
////		DCHECK( g.Arcs() == 4 && g.Vertices()==4 && g.degree(1)==2
////				,"Invalid graph undirected conversion");
////		GeneralGraph<int> gw("/data/m1_6.edge");
//
	}
};

struct	TestGeneralGraph {
	TestGeneralGraph() {
		GeneralGraph<> g("data/a.edge");
		g.toUndirected();
		DCHECK(g.Vertices() == 4, "Error: Wrong number of vertices! Expected 4, actual: "<<g.Vertices() )
		DCHECK(g.Edges() == 4, "Error: Wrong number of edges! Expected 4, actual: " << g.Edges() )
		g.writeToFile("data/a.bin");
		GeneralGraph<> g2("data/a.bin");
		DCHECK(g2==g,"Error: Invalid graph writing/reading!")
		g2.writeToFile("data/a2.adj");
		GeneralGraph<> g3("data/a2.adj");
		GeneralGraph<> fb("data/facebooks.edge");
		fb.toUndirected();
		fb.writeToFile("data/facebooks.bin");
		GeneralGraph<> f2("data/facebooks.bin");
		DCHECK(fb==f2,"Error: Invalid graph writing/reading!")
////		// Check simplify function
		int head[] = { 1, 1, 3, 1, 2 };
		int tail[] = { 2, 3, 1, 2, 1 };
		GeneralGraph<> gg(4, 5, head, tail);
		DCHECK( gg.Arcs() == 5 && gg.Vertices()==4 && gg.degree(1)==3
				,"Invalid graph initialization");
		gg.simplify();
		DCHECK( gg.Arcs() == 4 && gg.Vertices()==4 && gg.degree(1)==2
				,"Invalid graph simplification");
		gg.toUndirected();
				DCHECK( g.Arcs() == 8 && g.Vertices()==4 && g.degree(1)==2
						,"Invalid graph undirected conversion");

		EdgeMap<> fbm(fb);
		const int *u = fbm.getHead(0), *v = fbm.getTail(0);
		FOR(i, fbm.size()) {
			CHECK(fbm.edgeId(*u, *v) == i,"Error: Incorrect edgeId returned! ("
					<<*u<<";"<<*v<<")->"<<i)
			u++, v++;
		}
		GeneralGraph<> gw("data/m1_6.edge",true);
		gw.simplify();
		gw.toUndirected();
		gw.writeToFile("data/m1_6.bin");
		gw.writeToFile("data/m1_6.adj");
		GeneralGraph<> g5("data/m1_6.bin", true);
		GeneralGraph<> g6("data/m1_6.adj",true);
		DCHECK(gw==g5,"Error: weight bin");
		EdgeMap<> gwm(gw);
	}
};

TestUnweightedGraph	testUnweightedGraph;
TestGeneralGraph testGeneralGraph;
