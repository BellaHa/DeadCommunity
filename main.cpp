#include<iostream>
#include "SocialGraph.h"
#include "DCRgenerator.h"
//#include "ExactSolution.h"
#include "GreedySolution.h"
#include "SandwichSolution.h"
#include "BoundedThres.h"
#include "CompareGreedy.h"
#include "HighInfluence.h"
#include "HighTouch.h"
#include "HighBenefit.h"
#include "SSA.h"
#include <omp.h>
#include <time.h>
#include <fstream>
#include "Constant.h"
#include <stdlib.h>
#include <string>
#include "GIA.h"
#include "HighDegree.h"
#include "HyperGraph.h"
#include <string.h>


using namespace std;

#pragma warning(disable : 4996)

SocialGraph *g;
Graph *sg;
ofstream writefile;
string graphBinFile;
string seedFile;

/*
* convert from an integer to a string
*/
string intToStr(int i) {
    stringstream ss;
    ss << i;
    return ss.str();
}

/*
* convert from a strong to an integer
*/
unsigned int strToInt(string s) {
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
float getCurrentMemoryUsage() {

    string pid = intToStr(unsigned(getpid()));
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

long long addHyperedge(Graph &g, HyperGraph &hg, int t, long long int num, bool lt) {
    int numNodes = g.getSize();

    omp_set_num_threads(t);

    long long iter = 0;
    int c = 100;

#pragma omp parallel
    {
        vector<int> visit_mark(numNodes + 1, 0);
        vector<bool> visit(numNodes + 1, false);
        vector<unsigned int> link;
        if (lt == 0) {
            while (iter < num) {
                for (int i = 0; i < c; ++i) {
                    vector<int> he;
                    hg.pollingLT1(g, visit, visit_mark);
                }
#pragma omp atomic
                iter += c;
            }
        } else {
            while (iter < num) {
                for (int i = 0; i < c; ++i) {
                    vector<int> he;
                    hg.pollingIC1(g, visit, visit_mark);
                }

#pragma omp atomic
                iter += c;
            }

        }
    }
    hg.updateDeg();
    return hg.getNumEdge();
}

void buildSeedSet(HyperGraph &hg, vector<int> &seeds, unsigned int n, int k, vector<double> &degree) {
    long long i;
    unsigned int j, l, maxInd;
    vector<int> e, nList;

    vector<int> nodeDegree(n, 0);
    vector<int> indx(n, 0);
    for (j = 0; j < n; ++j) {
        indx[j] = j;
        nodeDegree[j] = hg.getNode(j + 1).size();
    }

    InfCost<int> hd(&nodeDegree[0]);
    MappedHeap<InfCost<int>> heap(indx, hd);
    long long numEdge = hg.getNumEdge();

    // check if an edge is removed
    vector<bool> edgeMark(numEdge, false);
    vector<bool> nodeMark(n + 1, true);

    double totalCost = 0;

    i = 1;
    // building each seed at a time
    while (totalCost < k && !heap.empty()) {
        maxInd = heap.pop() + 1;
        nodeMark[maxInd] = false;

        totalCost++;

        e = hg.getNode(maxInd);

        degree[i] = degree[i - 1] + nodeDegree[maxInd - 1];

        seeds.push_back(maxInd);
        for (j = 0; j < e.size(); ++j) {
            if (edgeMark[e[j]]) {
                continue;
            }

            nList = hg.getEdge(e[j]);
            for (l = 0; l < nList.size(); ++l) {
                nodeDegree[nList[l] - 1]--;
                if (nodeMark[nList[l]]) {
                    heap.heapify(nList[l] - 1);
                }
            }
            edgeMark[e[j]] = true;
        }
        i++;
    }

    vector<int>().swap(nodeDegree);
    vector<int>().swap(e);
    vector<int>().swap(nList);
    vector<int>().swap(indx);
    vector<bool>().swap(edgeMark);
}

bool
calculateInfluence(HyperGraph &hg, Graph &g, vector<int> &seeds, int t, double &deg, float epsilon, float delta, int m,
                   long long int maxSamples, int iter) {
    long long counter = 0;
    int n = g.getSize();
    unsigned k = seeds.size();
    vector<unsigned int> link(n + 1, seeds.size());
    double f =
            (log(6 / delta) + lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) * n / (k * log(6 * log2(n) / delta));
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
             << epsilon * sqrt(n * (1 + epsilon) / (degree * n * pow(2, iter - 1) / counter)) << " " << pow(2, iter - 1)
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

double getSSASolution(Graph *sg, vector<int> *sol, double *est) {
    sol->clear();
    srand(time(NULL));

    char *inFile = new char[graphBinFile.length() + 1];
    strcpy(inFile, graphBinFile.c_str());

    char *outFile = new char[seedFile.length() + 1];
    strcpy(outFile, seedFile.c_str());

    char *model = "IC";

    // Graph g;
    if (strcmp(model, "LT") == 0) {
        sg->readGraphLT(inFile);
    } else if (strcmp(model, "IC") == 0) {
        sg->readGraphIC(inFile);
    } else {
        printf("Incorrect model option!");
        return -1;
    }

    int n = sg->getSize();

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

    addHyperedge(*sg, hg, t, totalSamples, mo);
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
        if (calculateInfluence(hg, *sg, seeds, t, degree[k], epsilon, delta, mo, totalSamples, iter)) {
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
    return 0.;
}

void printResult(bool isScalable, bool isLargeFile) {
    vector<int> sol;
    sol.clear();

    // HighDegree hd(g);
    // long startHD = time(NULL);
    // double reHD = 0;
    // if (isScalable)
    //     hd.getSolution(&sol, &reHD);
    // else
    //     hd.getSolution2Step(&sol, &reHD);
    // long timeHD = time(NULL) - startHD;
    // cout << "HD: " << reHD << endl;
    // cout << "HD Time: " << timeHD << endl;
    // cout << sol.size() << endl;
    //
    // GIA gia(g);
    // long startGIA = time(NULL);
    // double reGIA = 0;
    // if (isScalable)
    //     gia.getSolution(&sol, &reGIA);
    // else
    //     gia.getSolution2Step(&sol, &reGIA);
    // long timeGIA = time(NULL) - startGIA;
    // cout << "GIA: " << reGIA << endl;
    // cout << "GIA Time: " << timeGIA << endl;
    // cout << sol.size() << endl;
    //
    // Constant::K = sol.size();
    //
    // SandwichSolution ubg(g);
    // long startUbg = time(NULL);
    // double reubg = 0;
    // if (isScalable)
    //     ubg.getSolution(&sol, &reubg);
    // else
    //     ubg.getSolution2Step(&sol, &reubg);
    // long timeUbg = time(NULL) - startUbg;
    // cout << "UBG: " << reubg << endl;
    // cout << "UBG Time: " << timeUbg << endl;

    // SSA ssa(g);
    double reSSA = 0;
    // ssa.graphBinFile = graphBinFile;
    // ssa.seedFile = seedFile;
    long startSSA = time(NULL);
    // ssa.getSolution(sg, &sol, &reSSA);
    getSSASolution(sg, &sol, &reSSA);
    long timeSSA = time(NULL) - startSSA;
    cout << "SSA: " << reSSA << endl;
    cout << "SSA Time: " << timeSSA << endl;

    // GreedySolution maf(g);
    // long startMaf = time(NULL);
    // double remaf = 0;
    // if (isScalable)
    //     maf.getSolution(&sol, &remaf);
    // else
    //     maf.getSolution2Step(&sol, &remaf);
    // long timeMaf = time(NULL) - startMaf;
    //
    // cout << "MAF: " << remaf << endl;
    // cout << "MAF Time: " << timeMaf << endl;

    //
    // CompareGreedy grd(g);
    // long startGrd = time(NULL);
    // double reGrd = 0;
    // if (isScalable)
    // grd.getSolution(&sol, &reGrd);
    // else
    // grd.getSolution2Step(&sol, &reGrd);
    // long timeGrd = time(NULL) - startGrd;
    //
    // long timeBt = 0;
    // double reBt = 0;
    // if (Constant::IS_BOUNDED_THRESHOLD && !isLargeFile) {
    //     BoundedThres bt(g);
    //     long startBt = time(NULL);
    //     reBt = 0;
    //     if (isScalable)
    //         bt.getSolution(&sol, &reBt);
    //     else
    //         bt.getSolution2Step(&sol, &reBt);
    //     timeBt = time(NULL) - startBt;
    // }
    //
    // HighBenefit hb(g);
    // long startHB = time(NULL);
    // double reHb = 0;
    // hb.getSolution(&sol, &reHb);
    // long timeHB = time(NULL) - startHB;
    //
    // cout << "HB: " << reHb << endl;

    // writefile << Constant::K << " \t " << Constant::COMMUNITY_POPULATION << "\t" << remaf << " " << timeMaf
    //           << "\t" << reubg << " " << ratio << " " << timeUbg
    //           //<< "\t" << reGrd << " " << timeGrd
    //           << (Constant::IS_BOUNDED_THRESHOLD && !isLargeFile ? "\t" + to_string(reBt) + " " + to_string(timeBt)
    //                                                              : "")
    //           //<< "\t" << reHi << " " << timeHi
    //           << "\t" << reHb << " " << timeHB
    //           << "\t" << reSSA << " " << timeSSA << endl;
}

void runExperiment(string input, string inputCommunity, int min, int max, int step,
                   bool isScalable = true, bool isBoundedThres = false, bool isLargeFile = false,
                   bool changeK = false, bool isDirected = true, bool isWeighted = false, bool isCommMM = true) {
    Constant::IS_BOUNDED_THRESHOLD = isBoundedThres;
    Constant::IS_WEIGHTED = isWeighted;
    Constant::COMMUNITY_POPULATION = 8;
    if (!isDirected) {
        if (!isLargeFile)
            g->readSocialGraphFromFile("../data/" + input);
        else
            g->readSocialGraphFromLargeFile("../data/" + input);
    } else {
        g->readSocialGraph("../data/" + input);
        g->readCommunityFile("../data/" + inputCommunity, isCommMM);
    }


    string outfilename = "../results/" + input + "_result_"
                         + (isScalable ? "scalable_" : "2step_")
                         + (Constant::IS_BOUNDED_THRESHOLD ? "boundedThres" : "freeThres")
                         + (changeK ? "_changeK" : "_changePop")
                         + (isWeighted ? "_weighted" : "_unweighted")
                         + (isCommMM ? "_isMM" : "_isClauset") + ".txt";
    writefile.open(outfilename);
    if (writefile.is_open()) {
        if (Constant::IS_BOUNDED_THRESHOLD && !isLargeFile)
            writefile << "k \t Pop \t maf \t ubg-ratio \t grd \t bt \t hb \t ssa" << endl;
        else
            writefile << "k \t Pop \t maf \t ubg-ratio \t grd \t hb \t ssa" << endl;

        if (changeK) {
            if (!isDirected)
                g->formCommunitiesFromActualCommunities();
            for (Constant::K = min; Constant::K <= max; Constant::K += step) {
                printResult(isScalable, isLargeFile);
            }
        } else {
            // for (Constant::COMMUNITY_POPULATION = min;
            //      Constant::COMMUNITY_POPULATION <= max; Constant::COMMUNITY_POPULATION += step) {
            g->formCommunitiesFromActualCommunities();
            printResult(isScalable, isLargeFile);
            // }
        }

        writefile.close();
    }
}

int main() {
    g = new SocialGraph();
    sg = new Graph();
    omp_set_num_threads(Constant::NUM_THREAD);

    vector<string> graphs{"facebook"};
    // vector<string> graphs{"facebook", "wiki","epinions", "dblp","pokec"};
    string pre = "../data/";
    for (int i = 0; i < graphs.size(); ++i) {
        string graph = pre + graphs[i];
        graphBinFile = graph + "SSA.bin";
        seedFile = graph + ".seeds";
        string input = graph + "E.txt";
        string inputCommunity = graph + "Comm.txt";
        runExperiment(input, inputCommunity, 4, 10, 2, true, false, false, false, true, false);
    }
    // runExperiment("wikiE.txt", "wikiComm.txt", 4, 10, 2, true, false, false, false, true, false);
    // runExperiment("epinionsE.txt", "epinionsComm.txt", 4, 10, 2, true, false, false, false, true, false);
    // runExperiment("dblpE.txt", "dblpComm.txt", 4, 10, 2, true, false, false, false, true, false);

    // Generate community
    // g->generateEdgeWeight("../data/facebook.txt", "../data/facebookE.txt");
    // g->readSocialGraph("../data/facebookE.txt", true);
    // g->generateFileIM("../data/facebookSSA.txt");
    // g->formCommunityModularity("../data/facebookComm.txt", "../data/facebook.adj", false);
    //
    // g->generateEdgeWeight("../data/wiki.txt", "../data/wikiE.txt");
    // g->readSocialGraph("../data/wikiE.txt", true);
    // g->generateFileIM("../data/wikiSSA.txt");
    // g->formCommunityModularity("../data/wikiComm.txt", "../data/wiki.adj", false);
    //
    // g->generateEdgeWeight("../data/epinions.txt", "../data/epinionsE.txt");
    // g->readSocialGraph("../data/epinionsE.txt", true);
    // g->generateFileIM("../data/epinionsSSA.txt");
    // g->formCommunityModularity("../data/epinionsComm.txt", "../data/epinions.adj", false);
    //
    // g->generateEdgeWeight("../data/dblp.txt", "../data/dblpE.txt");
    // g->readSocialGraph("../data/dblpE.txt", true);
    // g->generateFileIM("../data/dblpSSA.txt");
    // g->formCommunityModularity("../data/dblpComm.txt", "../data/dblp.adj", false);

    // g->generateEdgeWeight("../data/pokec.txt", "../data/pokecE.txt");
    // g->readSocialGraph("../data/pokecE.txt", true);
    // g->generateFileIM("../data/pokecSSA.txt");
    // g->formCommunityModularity("../data/pokecComm.txt", "../data/pokec.adj", false);

    delete g;
    return 0;
}