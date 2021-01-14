#pragma once
#include "Algorithm.h"
class SSA : public Algorithm
{
public:
    string seedFile;
    string graphBinFile;
	SSA(SocialGraph *g);
	~SSA();
	double getDeterministicSolution(vector<int> * sol);
	double getSolution(vector<int> * sol, double *est);

private:
    string graphSSAformat;
	const char * formatCmd;
};

