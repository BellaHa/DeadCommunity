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
    sol->clear();
    initiate();
    srand(time(NULL));

    char *inFile = new char[graphBinFile.length() + 1];
    strcpy(inFile, graphBinFile.c_str());

    char *outFile = new char[seedFile.length() + 1];
    strcpy(outFile, seedFile.c_str());

    char *model = "IC";

    Graph g;
    if (strcmp(model, "LT") == 0) {
        g.readGraphLT(inFile);
    } else if (strcmp(model, "IC") == 0) {
        g.readGraphIC(inFile);
    } else {
        printf("Incorrect model option!");
        return -1;
    }

    int n = g.getSize();

    float epsilon = Constant::EPSILON;

    float delta = Constant::DELTA;

    double k = Constant::K;

    int t = Constant::NUM_THREAD;

    HyperGraph hg(n);
    vector<double> degree(k + 1, 0);

    vector<int> seeds;

    double f =
            (log(6 / delta) + lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) * n / (k * log(6 * log2(n) / delta));

    double lambda = (2 + 2 * epsilon / 3) * log(3 * log2(f) / delta) / (epsilon * epsilon);

    long long int totalSamples = (long long int) lambda;
    cout << lambda << " " << totalSamples << endl;

    int mo = 0;
    if (strcmp(model, "IC") == 0)
        mo = 1;

    int iter = 1;

    addHyperedge(g, hg, t, totalSamples, mo);
    double nmax = (2 + 2 * epsilon / 3) * (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + log(6 / delta)) * n /
                  (epsilon * epsilon * k);

    clock_t start = clock();
    cout << totalSamples << " " << nmax << " " << lgamma(n + 1) << " " << lgamma(k + 1) << " " << lgamma(n - k + 1)
         << endl;

    while (totalSamples < nmax) {
        seeds.clear();
        totalSamples = hg.getNumEdge();
        cout << "Total Samples: " << totalSamples << endl;
        buildSeedSet(hg, seeds, n, k, degree);
        if (calculateInfluence(hg, g, seeds, t, degree[k], epsilon, delta, mo, totalSamples, iter)) {
            break;
        }
        iter++;
    }
    cout << "Seed Nodes: ";
    ofstream out(outFile);
    for (unsigned int i = 0; i < seeds.size(); ++i) {
        cout << seeds[i] << " ";
        out << seeds[i] << endl;
    }
    out.close();
    cout << endl;
    printf("Influence: %0.2lf\n", (double) degree[k] * n / totalSamples);
    cout << "Time: " << (float) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    cout << "Memory: " << getCurrentMemoryUsage() << " MB" << endl;
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

bool
SSA::calculateInfluence(HyperGraph &hg, Graph &g, vector<int> &seeds, int t, double &deg, float epsilon, float delta,
                        int m, long long int maxSamples, int iter) {
    long long counter = 0;
    int n = g.getSize();
    unsigned k = seeds.size();
    vector<unsigned int> link(n + 1, seeds.size());
    double f =
            (log(6 / delta) + lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) * n /
            (k * log(6 * log2(n) / delta));
    double lambda1 = 1 + (1 + epsilon) * (2 + 2 * epsilon / 3) * log(3 * log2(f) / delta) / (epsilon * epsilon);
    double degree = 0;
    for (unsigned int i = 0; i < k; ++i) {
        link[seeds[i]] = i;
    }
    vector<bool> maxSeed(t, false);

    omp_set_num_threads(t);
#pragma omp parallel
    {
        vector<bool> visit(n + 1, false);
        vector<int> visit_mark(n, 0);
        int id = omp_get_thread_num();

        if (m == 0) {
            while (counter < maxSamples) {
                maxSeed[id] = hg.pollingLT2(g, link, k, visit, visit_mark);
#pragma omp critical
                {
                    counter += 1;
                    if (maxSeed[id]) {
                        degree++;
                    }
                }
            }
        } else {
            while (counter < maxSamples) {
                maxSeed[id] = hg.pollingIC2(g, link, k, visit, visit_mark);
#pragma omp critical
                {
                    counter += 1;
                    if (maxSeed[id]) {
                        degree++;
                    }
                }
            }
        }
    }
//	cout << "Degree: " << degree << " " << counter << endl;

    if (degree >= lambda1) {
        double epsilon_1 = (deg * n / maxSamples) / (degree * n / counter) - 1;
        cout << "Epsilon_1 = " << epsilon_1 << endl;
        double epsilon_2 = epsilon * sqrt(n * (1 + epsilon) / (degree * n * pow(2, iter - 1) / counter));
        cout << "Epsilon_2 = " << epsilon_2 << " "
             << epsilon * sqrt(n * (1 + epsilon) / (degree * n * pow(2, iter - 1) / counter)) << " "
             << pow(2, iter - 1)
             << " " << pow(3, iter - 1) << endl;
        double epsilon_3 = epsilon * sqrt(n * (1 + epsilon) * (1 - 1 / exp(1) - epsilon) /
                                          ((1 + epsilon / 3) * degree * n * pow(2, iter - 1) / counter));
        cout << "Epsilon_3 = " << epsilon_3 << endl;
        cout << "Epsilon_t = " << (epsilon_1 + epsilon_2 + epsilon_1 * epsilon_2) * (1 - 1 / exp(1) - epsilon) +
                                  epsilon_3 * (1 - 1 / exp(1)) << endl;

        if ((epsilon_1 + epsilon_2 + epsilon_1 * epsilon_2) * (1 - 1 / exp(1) - epsilon) +
            epsilon_3 * (1 - 1 / exp(1)) <= epsilon) {
            return true;
        }
    }

    hg.updateDeg();
    return false;
}
