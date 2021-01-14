#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "Algorithm.h"
#include "hypergraph.h"

class SSA : public Algorithm {
public:
    string seedFile;
    string graphBinFile;

    SSA(SocialGraph *g);

    ~SSA();

    double getDeterministicSolution(vector<int> *sol);

    double getSolution(vector<int> *sol, double *est);

    bool
    calculateInfluence(HyperGraph &hg, Graph &g, vector<int> &seeds, int t, double &deg, float epsilon, float delta,
                       int m, long long int maxSamples, int iter);

private:
    string graphSSAformat;
    const char *formatCmd;
};

