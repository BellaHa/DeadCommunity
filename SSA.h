#pragma once

#include "Algorithm.h"

class SSA : public Algorithm {
public:
    string graphBinFile;
    string seedFile;

    SSA(SocialGraph *g);

    ~SSA();

    double getDeterministicSolution(vector<int> *sol);

    double getSolution(vector<int> *sol, double *est);
    double getSolutionBS(vector<int> *sol, double *est, int left, int right);

private:
    string graphSSAformat;
    const char *formatCmd;
};

