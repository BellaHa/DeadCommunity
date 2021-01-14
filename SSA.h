#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "Algorithm.h"
#include "sstream"
#include <process.h>
#include "mappedheap.hpp"
#include "HeapData.hpp"

class SSA : public Algorithm {
public:
    string seedFile;
    string graphBinFile;

    SSA(SocialGraph *g);

    ~SSA();

    double getDeterministicSolution(vector<int> *sol);

    double getSolution(vector<int> *sol, double *est);


    string intToStr(int i);

    unsigned int strToInt(string s);

    float getCurrentMemoryUsage();

private:
    string graphSSAformat;
    const char *formatCmd;
};

