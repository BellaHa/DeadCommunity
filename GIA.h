#pragma once

#include "Algorithm.h"

class GIA : public Algorithm {
public:
    double nMax;
    double n1;
    int iMax;
    double delta1;

    GIA(SocialGraph *g);

    ~GIA();

    double getDeterministicSolution(vector<int> *sol);

    double getDeterministicSolutionMig(vector<int> *sol);

    //double estimate(vector<int> * sol, double epsilon2, double delta, int tMax);
    double getSolution(vector<int> *sol, double *est);

    double getSolutionMig(vector<int> *sol, double *est);

    double getSolution2Step(vector<int> *sol, double *est);

    double estimate(vector<int> *sol, double epsilon, double delta, int tMax);

    double estimateInf(vector<int> *sol);


    void generateDCRgraphs(int number);

    void initiate();

    vector<double> currentLiveB; // store current live node in each dcr graph after each iteration in greedy

private:
    double getMarginalGain(int nodeId, vector<int> *sol);

    vector<vector<int>> currentLive; // store current live node in each dcr graph after each iteration in greedy
};
