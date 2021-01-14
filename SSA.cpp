#include "SSA.h"
#include <iostream>
#include <fstream>
#include <string>

SSA::SSA(SocialGraph *g) : Algorithm(g) {
    // graphSSAformat = "graphSSA.txt";
    // // g->generateFileIM(graphSSAformat);
    // graphBinFile = "graphSSA.bin";
    //
    // string tmp = " ../SSA_release_2.0/DSSA/el2bin " + graphSSAformat + " " + graphBinFile;
    // formatCmd = tmp.c_str();
    // cout << formatCmd << endl;
    // system(formatCmd);
    //
    // seedFile = "ssa.seeds";
}

SSA::~SSA() {
}

double SSA::getDeterministicSolution(vector<int> *sol) {
    return 0.0;
}

double SSA::getSolution(vector<int> *sol, double *est) {
    // sol->clear();
    // srand(time(NULL));
    //
    // char *inFile = new char[graphBinFile.length() + 1];
    // strcpy(inFile, graphBinFile.c_str());
    //
    // char *outFile = new char[seedFile.length() + 1];
    // strcpy(outFile, seedFile.c_str());
    //
    // char *model = "IC";
    //
    // Graph *sg = new Graph();
    // if (strcmp(model, "LT") == 0) {
    //     sg->readGraphLT(inFile);
    // } else if (strcmp(model, "IC") == 0) {
    //     sg->readGraphIC(inFile);
    // } else {
    //     printf("Incorrect model option!");
    //     return -1;
    // }
    //
    // int n = sg->getSize();
    //
    // float epsilon = Constant::EPSILON;
    //
    // float delta = Constant::DELTA;
    //
    // double k = Constant::K;
    //
    // int t = Constant::NUM_THREAD;
    //
    // HyperGraph hg(n);
    // vector<double> degree(k + 1, 0);
    //
    // vector<int> seeds;
    //
    // double f =
    //         (log(6 / delta) + lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) * n / (k * log(6 * log2(n) / delta));
    //
    // double lambda = (2 + 2 * epsilon / 3) * log(3 * log2(f) / delta) / (epsilon * epsilon);
    //
    // long long int totalSamples = (long long int) lambda;
    // cout << lambda << " " << totalSamples << endl;
    //
    // int mo = 0;
    // if (strcmp(model, "IC") == 0)
    //     mo = 1;
    //
    // int iter = 1;
    //
    // addHyperedge(*sg, hg, t, totalSamples, mo);
    // double nmax = (2 + 2 * epsilon / 3) * (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + log(6 / delta)) * n /
    //               (epsilon * epsilon * k);
    //
    // clock_t start = clock();
    // cout << totalSamples << " " << nmax << " " << lgamma(n + 1) << " " << lgamma(k + 1) << " " << lgamma(n - k + 1)
    //      << endl;
    //
    // while (totalSamples < nmax) {
    //     seeds.clear();
    //     totalSamples = hg.getNumEdge();
    //     cout << "Total Samples: " << totalSamples << endl;
    //     buildSeedSet(hg, seeds, n, k, degree);
    //     if (calculateInfluence(hg, *sg, seeds, t, degree[k], epsilon, delta, mo, totalSamples, iter)) {
    //         break;
    //     }
    //     iter++;
    // }
    // cout << "Seed Nodes: ";
    // ofstream out(outFile);
    // for (unsigned int i = 0; i < seeds.size(); ++i) {
    //     cout << seeds[i] << " ";
    //     out << seeds[i] << endl;
    // }
    // out.close();
    // cout << endl;
    // printf("Influence: %0.2lf\n", (double) degree[k] * n / totalSamples);
    // cout << "Time: " << (float) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    // cout << "Memory: " << getCurrentMemoryUsage() << " MB" << endl;

    return 1.;
    //
    // string tmp2 = "../SSA_release_2.0/DSSA/DSSA -i " + graphBinFile + " -o "
    //               + seedFile + " -k " + to_string(Constant::K) + " -epsilon "
    //               + to_string(Constant::EPSILON) + " -delta " + to_string(Constant::DELTA);
    // const char *runSSAcmd = tmp2.c_str();
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
    // return 1;
}

/*
* convert from an integer to a string
*/
string SSA::intToStr(int i) {
    stringstream ss;
    ss << i;
    return ss.str();
}

/*
* convert from a strong to an integer
*/
unsigned int SSA::strToInt(string s) {
    unsigned int i;
    istringstream myStream(s);

    if (myStream >> i) {
        return i;
    } else {
        cout << "String " << s << " is not a number." << endl;
        return atoi(s.c_str());
    }
    return i;
}

/*
* measure the consumed memory
*/
float SSA::getCurrentMemoryUsage() {

    string pid = intToStr(unsigned(_getpid()));
    string outfile = "tmp_" + pid + ".txt";
    string command = "pmap " + pid + " | grep -i Total | awk '{print $2}' > " + outfile;
    system(command.c_str());

    string mem_str;
    ifstream ifs(outfile.c_str());
    std::getline(ifs, mem_str);
    ifs.close();

    mem_str = mem_str.substr(0, mem_str.size() - 1);
    float mem = (float) strToInt(mem_str);

    command = "rm " + outfile;
    system(command.c_str());

    return mem / 1024;

    return 0;
}

