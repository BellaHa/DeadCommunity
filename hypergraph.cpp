#include "hypergraph.h"

long long addHyperedge(Graph & g, HyperGraph & hg, int t, long long num, bool lt)
{
    int numNodes = g.getSize();

    omp_set_num_threads(t);

    long long iter = 0;
    int c = 100;

#pragma omp parallel
    {
        vector<int> visit_mark(numNodes+1,0);
        vector<bool> visit(numNodes+1,false);
        vector<unsigned int> link;
        if (lt == 0){
            while (iter < num){
                for (int i = 0; i < c; ++i){
                    vector<int> he;
                    hg.pollingLT1(g,visit,visit_mark);
                }
#pragma omp atomic
                iter += c;
            }
        } else {
            while (iter < num){
                for (int i = 0; i < c; ++i){
                    vector<int> he;
                    hg.pollingIC1(g,visit,visit_mark);
                }

#pragma omp atomic
                iter += c;
            }

        }
    }
    hg.updateDeg();
    return hg.getNumEdge();
}

/*
* find seed nodes procedure using greedy algorithm
*/
void buildSeedSet(HyperGraph & hg, vector<int> & seeds, unsigned int n, int k, vector<double> &degree)
{
    long long i;
    unsigned int j,l,maxInd;
    vector<int> e, nList;

    vector<int> nodeDegree(n,0);
    vector<int> indx(n,0);
    for (j = 0; j < n; ++j){
        indx[j] = j;
        nodeDegree[j] = hg.getNode(j+1).size();
    }

    InfCost<int> hd(&nodeDegree[0]);
    MappedHeap<InfCost<int> > heap(indx,hd);
    long long numEdge = hg.getNumEdge();

    // check if an edge is removed
    vector<bool> edgeMark(numEdge, false);
    vector<bool> nodeMark(n+1, true);

    double totalCost = 0;

    i=1;
    // building each seed at a time
    while(totalCost < k && !heap.empty()){
        maxInd = heap.pop()+1;
        nodeMark[maxInd] = false;

        totalCost++;

        e = hg.getNode(maxInd);

        degree[i] = degree[i-1]+nodeDegree[maxInd-1];

        seeds.push_back(maxInd);
        for (j = 0; j < e.size(); ++j){
            if (edgeMark[e[j]]){
                continue;
            }

            nList = hg.getEdge(e[j]);
            for (l = 0; l < nList.size(); ++l){
                nodeDegree[nList[l]-1]--;
                if (nodeMark[nList[l]]){
                    heap.heapify(nList[l]-1);
                }
            }
            edgeMark[e[j]] = true;
        }
        i++;
    }

    vector<int>().swap(nodeDegree);
    vector<int>().swap(e);
    vector<int>().swap(nList);
    vector<int>().swap(indx);
    vector<bool>().swap(edgeMark);
}
