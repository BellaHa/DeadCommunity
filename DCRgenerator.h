#pragma once
#include "SocialGraph.h"
#include "DCRgraph.h"

class DCRgenerator
{
public:
	DCRgenerator();
	DCRgenerator(SocialGraph* g);
	~DCRgenerator();

	void setSocialGraph(SocialGraph * g);

	DCRgraph* generateDCRgraph();
    DCRgraph* generateDCRgraphMig();

private:
	SocialGraph *g;

	DCRgraph * generateDCRgraphIC();
    DCRgraph * generateDCRgraphICMig();
	DCRgraph * generateDCRgraphLT();

	void dfs(int u, vector<int> * reachable, map<int, vector<int>> * mapNeighbors);
	Common * commonInstance;

};

