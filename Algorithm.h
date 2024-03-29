#pragma once

#include<vector>
#include"DCRgraph.h"
#include"SocialGraph.h"
#include "DCRgenerator.h"
#include "Constant.h"
#include "Common.h"
#include <algorithm>

using namespace std;

class Algorithm {
public:
    Algorithm(SocialGraph *g);

    Algorithm();

    ~Algorithm();

    virtual double getDeterministicSolution(vector<int> *sol) = 0;

    virtual double estimate(vector<int> *sol, double epsilon2, double delta, int tMax);

    virtual double getSolution(vector<int> *sol, double *est);

    virtual double getSolution2Step(vector<int> *sol, double *est);

    virtual double estimateInf(vector<int> *sol);

    double estimateInfMig(vector<int> *sol);

    virtual void generateDCRgraphs(int number);

    void generateDCRgraphsMig(int number);

    void initiate();

    virtual void initiateMig();

    double calculateCost(vector<int> sol);

    double calculateCostMig(vector<int> sol);

    map<int, double> intialGainB; // initial fraction inf of each nodes, USED in sandwich solution
    double nMax;
    double n1;
    int iMax;
    double delta1;
    double c;

protected:
    SocialGraph *g;
    vector<DCRgraph *> dcrSet;
    DCRgenerator gen;
    Common *commonInstance;

    map<int, double> intialGain; // initial fraction inf of each nodes, USED in sandwich solution
    map<int, int> initialDead; // initital number of influenced graphs, USED in greedy solution
    map<int, double> initialOtherCommunityGain;  // initial fraction inf of each nodes to samples that have source community is not its community, USED in MAF
    vector<int> indx;
    map<int, int> mapNodeIdx;

    double D, e1, e2, e3; // D - lambda, epsilon 1 2 3
    int rMax; // maximum number of ri graph

    void clear();

    void clearMig();

    bool isMaf = false;
};

