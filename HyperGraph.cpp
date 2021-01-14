#include "HyperGraph.h"

inline int HyperGraph::randIndex_lin(const vector <UI> &w, unsigned int si) {
    UI ranNum = sfmt_genrand_uint32(&sfmtSeed);
    if (si <= 1 || ranNum > w[si - 1])
        return -1;

    for (unsigned int i = 1; i < si; ++i) {
        if (ranNum <= w[i])
            return i;
    }
    return -1;
}

// choose a random live edge in LT model based on binary search
inline int HyperGraph::randIndex_bin(const vector <UI> &w, unsigned int si) {
    UI ran = sfmt_genrand_uint32(&sfmtSeed);
    if (si <= 1 || ran > w[si - 1])
        return -1;
    int left = 1;
    int right = si - 1;
    int prob;
    for (unsigned int i = 0; i < si; ++i) {
        prob = (left + right) / 2;
        if (w[prob - 1] > ran) {
            right = prob - 1;
            continue;
        }
        if (w[prob] <= ran) {
            left = prob + 1;
            continue;
        }
        break;
    }
    return prob;
}

inline int HyperGraph::randIndex_dart(const vector <UI> &w, unsigned int si) {
    int prob = 0;
    while (prob == 0) {
        UI ran = sfmt_genrand_uint32(&sfmtSeed) % si;
        UI ran2 = sfmt_genrand_uint32(&sfmtSeed);
        if (w[ran] >= ran2)
            prob = ran;
    }
    return prob;
}

HyperGraph::HyperGraph(unsigned int n) {
    sfmt_init_gen_rand(&sfmtSeed, rand());
    node_edge = vector < vector < int > > (n + 1);
    maxDegree = 0;
    numNodes = n;
    curEdge = 0;
}

void HyperGraph::updateDeg() {
    unsigned int num = edge_node.size();
    for (unsigned int i = curEdge; i < num; ++i) {
        unsigned int num2 = edge_node[i].size();
        for (unsigned int j = 0; j < num2; ++j) {
            node_edge[edge_node[i][j]].push_back(i);
        }
    }
    curEdge = edge_node.size();
}

void HyperGraph::updateEdge() {
    curEdge = edge_node.size();
}

/*
* Add a hyperedge into the hypergraph
*/
void HyperGraph::addEdge(vector<int> &edge) {
    edge_node.push_back(edge);
    unsigned int ind = edge_node.size() - 1;
    for (unsigned int i = 0; i < edge.size(); ++i)
        node_edge[edge[i]].push_back(ind);
}

/*
* Add a hyperedge into the hypergraph while keeping track of the node with max degree
*/
void HyperGraph::addEdgeD(vector<int> &edge) {
    edge_node.push_back(edge);
    int ind = edge_node.size() - 1;
    for (unsigned int i = 0; i < edge.size(); ++i) {
        node_edge[edge[i]].push_back(ind);
        if (node_edge[edge[i]].size() > maxDegree)
            maxDegree = node_edge[edge[i]].size();
    }
}

/*
* get an edge from the hypergraph
*/
const vector<int> &HyperGraph::getEdge(int e) const {
    return edge_node[e];
}

const vector<int> &HyperGraph::getEdge(int e) {
    return edge_node[e];
}

/*
* get the list of hyperedges incident to node n
*/
const vector<int> &HyperGraph::getNode(int n) const {
    return node_edge[n];
}

const vector<int> &HyperGraph::getNode(int n) {
    return node_edge[n];
}

/*
* get the number of hyperedges
*/
int HyperGraph::getNumEdge() const {
    return edge_node.size();
}

/*
* get the maximum degree
*/
int HyperGraph::getMaxDegree() {
    return maxDegree;
}

/*
* remove all the hyperedges
*/
void HyperGraph::clearEdges() {
    edge_node.clear();
    node_edge.clear();
    cout << "clear edges!" << endl;
    maxDegree = 0;
}

/*
* polling process under LT model
*/
bool HyperGraph::pollingLT2(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit,
                            vector<int> &visit_mark) {
    unsigned int i;
    bool t = false;
    unsigned int gSize = g.getSize();
    unsigned int cur = sfmt_genrand_uint32(&sfmtSeed) % gSize + 1;
    unsigned int num_marked = 0;
    for (i = 0; i < gSize; ++i) {
        if (visit[cur] == true) break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;
        if (link[cur] < k)
            t = true;
        int ind;
        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
            break;

        cur = g.adjList[cur][ind - 1];
    }
    edge_node.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + num_marked));
    for (i = 0; i < num_marked; ++i) {
        visit[visit_mark[i]] = false;
    }
    return t;
}

bool HyperGraph::pollingLT(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit,
                           vector<int> &visit_mark) {
    unsigned int i;
    bool t = false;
    unsigned int gSize = g.getSize();
    unsigned int cur = sfmt_genrand_uint32(&sfmtSeed) % gSize + 1;
    unsigned int num_marked = 0;
    for (i = 0; i < gSize; ++i) {
        if (link[cur] < k) {
            t = true;
            break;
        }
        if (visit[cur] == true) break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;
        int ind;

        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
            break;

        cur = g.adjList[cur][ind - 1];
    }
    for (i = 0; i < num_marked; ++i) {
        visit[visit_mark[i]] = false;
    }
    return t;
}

void HyperGraph::pollingLT1(Graph &g, vector<bool> &visit, vector<int> &visit_mark) {
    unsigned int i;
    unsigned int gSize = g.getSize();
    unsigned int cur = sfmt_genrand_uint32(&sfmtSeed) % gSize + 1;
    unsigned int num_marked = 0;
    for (i = 0; i < gSize; ++i) {
        if (visit[cur] == true) break;
        visit[cur] = true;
        visit_mark[num_marked] = cur;
        num_marked++;
        const vector<int> &neigh = g[cur];
        int ind;

        if (g.weights[cur].size() >= 32)
            ind = randIndex_bin(g.weights[cur], g.node_deg[cur]);
        else
            ind = randIndex_lin(g.weights[cur], g.node_deg[cur]);

        if (ind == -1)
            break;

        cur = neigh[ind - 1];
    }
    edge_node.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + num_marked));
    for (i = 0; i < num_marked; ++i) {
        visit[visit_mark[i]] = false;
    }
}

void HyperGraph::pollingIC1(Graph &g, vector<bool> &visit, vector<int> &visit_mark) {
    int i;
    unsigned int cur = sfmt_genrand_uint32(&sfmtSeed) % (g.getSize()) + 1;
    int num_marked = 1;
    int curPos = 0;
    visit[cur] = true;
    visit_mark[0] = cur;
    while (curPos < num_marked) {
        cur = visit_mark[curPos];
        curPos++;
        const vector <UI> &w = g.getWeight(cur);
        const vector<int> &neigh = g[cur];
        for (i = 0; i < g.node_deg[cur]; ++i) {
            if (sfmt_genrand_uint32(&sfmtSeed) < w[i + 1]) {
                if (!visit[neigh[i]]) {
                    visit[neigh[i]] = true;
                    visit_mark[num_marked] = neigh[i];
                    num_marked++;
                }
            }
        }
    }
    edge_node.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + num_marked));
    for (i = 0; i < num_marked; ++i) {
        visit[visit_mark[i]] = false;
    }
}

/*
* polling process under IC model
*/
bool HyperGraph::pollingIC2(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit,
                            vector<int> &visit_mark) {
    int i;
    unsigned int cur = sfmt_genrand_uint32(&sfmtSeed) % (g.getSize()) + 1;
    int curPos = 0;
    int num_marked = 1;
    visit[cur] = true;
    visit_mark[0] = cur;
    bool t = false;
    while (curPos < num_marked) {
        cur = visit_mark[curPos];
        curPos++;
        if (link[cur] < k)
            t = true;
        const vector <UI> &w = g.getWeight(cur);
        const vector<int> &neigh = g[cur];
        for (i = 0; i < g.node_deg[cur]; ++i) {
            if (sfmt_genrand_uint32(&sfmtSeed) < w[i + 1]) {
                if (!visit[neigh[i]]) {
                    visit[neigh[i]] = true;
                    visit_mark[num_marked] = neigh[i];
                    num_marked++;
                }
            }
        }
    }
    edge_node.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + num_marked));

    for (i = 0; i < num_marked; ++i) {
        visit[visit_mark[i]] = false;
    }
    return t;
}

bool HyperGraph::pollingIC(Graph &g, vector<unsigned int> &link, unsigned int k, vector<bool> &visit,
                           vector<int> &visit_mark) {
    int i;
    unsigned int cur = sfmt_genrand_uint32(&sfmtSeed) % (g.getSize()) + 1;
    int curPos = 0;
    int num_marked = 1;
    visit[cur] = true;
    visit_mark[0] = cur;
    bool t = false;
    while (curPos < num_marked) {
        cur = visit_mark[curPos];
        curPos++;
        if (link[cur] < k) {
            t = true;
            break;
        }
        const vector <UI> &w = g.getWeight(cur);
        const vector<int> &neigh = g[cur];
        for (i = 0; i < g.node_deg[cur]; ++i) {
            if (sfmt_genrand_uint32(&sfmtSeed) < w[i + 1]) {
                if (!visit[neigh[i]]) {
                    visit[neigh[i]] = true;
                    visit_mark[num_marked] = neigh[i];
                    num_marked++;
                }
            }
        }
    }
    for (i = 0; i < num_marked; ++i) {
        visit[visit_mark[i]] = false;
    }
    return t;
}
