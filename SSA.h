#pragma once

#include "Algorithm.h"
#include <time.h>
#include "StopWatch.h"

class SSA : public Algorithm {
public:
    string graphBinFile;
    string seedFile;
    double bsTime = 0;

    SSA(SocialGraph *g);

    ~SSA();

    double getDeterministicSolution(vector<int> *sol);

    double getSolution(vector<int> *sol, double *est);

    double getSolutionBS(vector<int> *sol, double *est, int left, int right);

private:
    string graphSSAformat;
    const char *formatCmd;
};

