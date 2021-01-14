#pragma once

#include <stdint.h>
#include <vector>
#include <iostream>

using namespace std;

typedef uint32_t UI;
typedef uint64_t ULL;

class Graph {
public:
    UI UI_MAX = 4294967295U;
    ULL ULL_MAX = 18446744073709551615ULL;

    // number of nodes
    unsigned int numNodes;
    // number of edges
    unsigned int numEdges;
    // adjacency list
    vector<int> testV;
    std::vector <std::vector<int>> adjList;
    std::vector<int> node_deg;
    std::vector <std::vector<UI>> weights;
    std::vector<bool> dart;

    Graph();
    ~Graph();

    // get a vector of neighbours of node u
    const std::vector<int> &operator[](int u) const;

    // return weights of neighbours of node u
    const std::vector <UI> &getWeight(int u) const;

    // get a vector of neighbours of node u
    const std::vector<int> &operator[](int u);

    // return weights of neighbours of node u
    const std::vector <UI> &getWeight(int u);

    // get degree of node u
    int getDegree(int u) const;

    // get size of the graph
    int getSize() const;

    // get number of edges
    int getEdge() const;

    // read graph from a file
    void readGraphLT(const char *filename);

    // read graph from a file
    void readGraphIC(const char *filename);

    // write the graph to file
    void writeToFile(const char *filename);
};