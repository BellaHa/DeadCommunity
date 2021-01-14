#pragma once

#include "Graph.h"
#include "SFMT.h"

class HyperGraph {
public:
    // store the edges that a node is incident to
    std::vector<std::vector<int>> node_edge;
    // store hyperedges
    std::vector<std::vector<int>> edge_node;
    unsigned int curEdge;
    unsigned int maxDegree;
    unsigned int numNodes;
    sfmt_t sfmtSeed;

    // randomly selects a node in a list according to the weight distribution
    inline int randIndex_bin(const std::vector<UI> &w, unsigned int si);

    inline int randIndex_lin(const std::vector<UI> &w, unsigned int si);

    inline int randIndex_dart(const std::vector<UI> &w, unsigned int si);

    HyperGraph(unsigned int n);

    void updateDeg();

    void updateEdge();

    void addEdge(std::vector<int> &edge);

    void addEdgeD(std::vector<int> &edge);

    int getMaxDegree();

    const std::vector<int> &getEdge(int e) const;

    const std::vector<int> &getEdge(int e);

    const std::vector<int> &getNode(int n) const;

    const std::vector<int> &getNode(int n);

    int getNumEdge() const;

    void clearEdges();

    // different polling processes for different purposes: model and returned output
    void pollingLT1(Graph &g, std::vector<bool> &visit, std::vector<int> &mark_visit);

    bool pollingLT2(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::vector<bool> &visit,
                    std::vector<int> &mark_visit);

    bool pollingLT(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::vector<bool> &visit,
                   std::vector<int> &mark_visit);

    void pollingIC1(Graph &g, std::vector<bool> &visit, std::vector<int> &visit_mark);

    bool pollingIC2(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::vector<bool> &visit,
                    std::vector<int> &visit_mark);

    bool pollingIC(Graph &g, std::vector<unsigned int> &link, unsigned int k, std::vector<bool> &visit,
                   std::vector<int> &visit_mark);
};
