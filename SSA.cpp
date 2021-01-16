#include "SSA.h"
#include <iostream>
#include <fstream>
#include <string>

SSA::SSA(SocialGraph *g) : Algorithm(g) {
    graphSSAformat = "graphSSA.txt";
    // g->generateFileIM(graphSSAformat);
    graphBinFile = "graphSSA.bin";
    //
    // string tmp = " ../SSA_release_2.0/DSSA/el2bin " + graphSSAformat + " " + graphBinFile;
    // formatCmd = tmp.c_str();
    // system(formatCmd);
    //
    seedFile = "ssa.seeds";
}

SSA::~SSA() {
}

double SSA::getDeterministicSolution(vector<int> *sol) {
    return 0.0;
}

double SSA::getSolution(vector<int> *sol, double *est) {
    sol->clear();
    initiate();

    string tmp2 = "D:/DSSA/cmake-build-release/DSSA -i " + graphBinFile + " -o "
                  + seedFile + " -k " + to_string(Constant::K) + " -epsilon "
                  + to_string(Constant::EPSILON) + " -delta " + to_string(Constant::DELTA);
    const char *runSSAcmd = tmp2.c_str();
    cout << runSSAcmd << endl;
    system(runSSAcmd);

    ifstream inputFile;
    inputFile.open(seedFile);
    if (inputFile) {
        double et;
        inputFile >> et;
        bsTime += et;
        vector<int> *listNodes = g->getListNodeIds();
        int nodeIdx;
        while (inputFile >> nodeIdx) {
            int nodeId = listNodes->at(nodeIdx - 1);
            sol->push_back(nodeId);
        }
        inputFile.close();
    }

    StopWatch sw;
    sw.start();
    *est = estimate(sol, Constant::EPSILON, Constant::DELTA, 100000000);
    sw.stop();
    bsTime += sw.getSeconds();
    clear();
    return 1;
}

double SSA::getSolutionBS(vector<int> *sol, double *est, int left, int right) {
    if (left >= right) return 1;
    vector<int> sol1;
    double est1;
    double K = (double) g->getNumberOfCommunities();
    double e = Constant::EPSILON;
    Constant::K = left + (right - left) / 2;
    getSolution(&sol1, &est1);
    cout << "K: " << Constant::K << endl;
    cout << "size: " << sol1.size() << endl;
    // est1 = Algorithm::estimate(&sol1, Constant::EPSILON, Constant::DELTA, 100000000);
    cout << "est1: " << est1 << endl;
    cout << "ratio: " << est1 * 100. / g->getNumberOfCommunities() << endl;
    cout << "-----\n";
    if (est1 >= K - e * K) {
        *sol = sol1;
        *est = est1;
        if (sol1.size() < Constant::K) {
            Constant::K = sol1.size();
        }
        return getSolutionBS(sol, est, left, Constant::K);
    }
    return getSolutionBS(sol, est, Constant::K + 1, right);

    // sol->clear();
    // initiate();
    //
    // string tmp2 = "D:/DSSA/cmake-build-release/DSSA -i " + graphBinFile + " -o "
    //               + seedFile + " -k " + to_string(Constant::K) + " -epsilon "
    //               + to_string(Constant::EPSILON) + " -delta " + to_string(Constant::DELTA);
    // const char *runSSAcmd = tmp2.c_str();
    // cout << runSSAcmd << endl;
    // system(runSSAcmd);
    //
    // ifstream inputFile;
    // inputFile.open(seedFile);
    // if (inputFile) {
    //     vector<int> *listNodes = g->getListNodeIds();
    //     int nodeIdx;
    //     while (inputFile >> nodeIdx) {
    //         int nodeId = listNodes->at(nodeIdx - 1);
    //         sol->push_back(nodeId);
    //     }
    //     inputFile.close();
    // }
    //
    // *est = estimate(sol, Constant::EPSILON, Constant::DELTA, 100000000);
    // clear();
    return 1;
}
