#include "Graph.h"
const vector<int> &Graph::operator[](int u) const {
    return adjList[u];
}

const vector<int> &Graph::operator[](int u) {
    return adjList[u];
}

const vector<UI> &Graph::getWeight(int u) const {
    return weights[u];
}

const vector<UI> &Graph::getWeight(int u) {
    return weights[u];
}

/*
* get degree of node u
*/
int Graph::getDegree(int u) const {
    return adjList[u].size();
}

/*
* get the number of nodes
*/
int Graph::getSize() const {
    return numNodes;
}

/*
* get the number of edges
*/
int Graph::getEdge() const {
    return numEdges;
}

/*
* read binary graph input for LT model
* difference between LT and IC is for LT we accumulate the weights for fast choosing a random node
*/
void Graph::readGraphLT(const char *filename) {
    FILE *pFile;
    pFile = fopen(filename, "rb");
    fread(&numNodes, sizeof(int), 1, pFile);
    fread(&numEdges, sizeof(long long), 1, pFile);
    node_deg = vector<int>(numNodes + 1);
    dart = vector<bool>(numNodes + 1, false);
    fread(&node_deg[1], sizeof(int), numNodes, pFile);

    vector<int> a;
    vector<UI> b;
    adjList.push_back(a);
    weights.push_back(b);

    for (unsigned int i = 1; i <= numNodes; ++i) {
        vector<int> tmp(node_deg[i]);
        fread(&tmp[0], sizeof(int), node_deg[i], pFile);

        adjList.push_back(tmp);
    }

    for (unsigned int i = 1; i <= numNodes; ++i) {
        vector<float> tmp(node_deg[i] + 1, 0);
        vector<UI> tmp1(node_deg[i] + 1, 0);
        fread(&tmp[1], sizeof(float), node_deg[i], pFile);

        for (int j = 1; j < node_deg[i] + 1; ++j) {
            tmp[j] += tmp[j - 1];
            if (tmp[j] >= 1) {
                tmp1[j] = UI_MAX;
            } else {
                tmp1[j] = tmp[j] * UI_MAX;
            }
        }

        weights.push_back(tmp1);
        node_deg[i]++;
    }
}

/*
* read input graph for IC model
*/
void Graph::readGraphIC(const char *filename) {
    FILE *pFile;
    pFile = fopen(filename, "rb");
    fread(&numNodes, sizeof(int), 1, pFile);
    fread(&numEdges, sizeof(long long), 1, pFile);
    node_deg = vector<int>(numNodes + 1);
    fread(&node_deg[1], sizeof(int), numNodes, pFile);

    testV = vector<int>(10);
    testV.resize(1);
    vector<int> a;
    vector<UI> b;
    adjList.push_back(a);
    weights.push_back(b);

    for (unsigned int i = 1; i <= numNodes; ++i) {
        vector<int> tmp(node_deg[i]);
        fread(&tmp[0], sizeof(int), node_deg[i], pFile);
        adjList.push_back(tmp);
    }

    for (unsigned int i = 1; i <= numNodes; ++i) {
        vector<float> tmp(node_deg[i] + 1, 0);
        vector<UI> tmp1(node_deg[i] + 1, 0);
        fread(&tmp[1], sizeof(float), node_deg[i], pFile);

        for (int j = 1; j < node_deg[i] + 1; ++j) {
            tmp1[j] = tmp[j] * UI_MAX;
        }

        if (tmp1[node_deg[i]] <= 0)
            tmp1[node_deg[i]] = UI_MAX;

        weights.push_back(tmp1);
    }
}

void Graph::writeToFile(const char *filename) {/*
	ofstream output(filename);
	for (unsigned int i = 0; i < numNodes; ++i){
		for (unsigned int j = 0; j < adjList[i].size(); ++j){
			if (adjList[i][j] > i){
				output << adjList[i][j] << " " << i << " " << weights[i][j] << endl;
			}
		}
	}
	output.close();
*/
}

Graph::Graph() {

}

Graph::~Graph() {

}
