#include <iostream>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <string>
#include "SocialGraph.h"
#include "DCRgenerator.h"
#include "GreedySolution.h"
#include "SandwichSolution.h"
#include "SSA.h"
#include "GIA.h"
#include "HighDegree.h"
#include "BoundedThres.h"
#include "CompareGreedy.h"
#include "HighInfluence.h"
#include "HighTouch.h"
#include "HighBenefit.h"
#include "Constant.h"

using namespace std;

#pragma warning(disable : 4996)

SocialGraph *g;
ofstream writefile;
string graphBinFile;
string seedFile;

void printResult(bool isScalable, bool isLargeFile) {
    vector<int> sol;
    sol.clear();

    cout << "GIA..." << endl;
    GIA gia(g);
    double reGIA = 0;
    long startGIA = time(NULL);
    gia.getSolution(&sol, &reGIA);
    long timeGIA = time(NULL) - startGIA;
    double costGIA = gia.calculateCost(sol);
    cout << "GIA: " << reGIA << endl;
    cout << "GIA Cost: " << costGIA << endl;
    cout << "GIA Time: " << timeGIA << endl;

    cout << "HD..." << endl;
    HighDegree hd(g);
    double reHD = 0;
    long startHD = time(NULL);
    hd.getSolution(&sol, &reHD);
    long timeHD = time(NULL) - startHD;
    double costHD = hd.calculateCost(sol);
    cout << "HD: " << reHD << endl;
    cout << "HD Cost: " << costHD << endl;
    cout << "HD Time: " << timeHD << endl;

    cout << "MAF..." << endl;
    GreedySolution maf(g);
    double remaf = 0;
    long startMaf = time(NULL);
    maf.getSolution(&sol, &remaf);
    long timeMaf = time(NULL) - startMaf;
    double costMaf = maf.calculateCost(sol);
    cout << "MAF: " << remaf << endl;
    cout << "MAF Cost: " << costMaf << endl;
    cout << "MAF Time: " << timeMaf << endl;

    cout << "UBG..." << endl;
    SandwichSolution ubg(g);
    double reubg = 0;
    long startUbg = time(NULL);
    ubg.getSolutionFast(&sol, &reubg);
    long timeUbg = time(NULL) - startUbg;
    double costUbg = ubg.calculateCost(sol);
    cout << "UBG: " << reubg << endl;
    cout << "UBG Cost: " << costUbg << endl;
    cout << "UBG Time: " << timeUbg << endl;

    cout << "DSSA..." << endl;
    SSA ssa(g);
    double reSSA = 0;
    ssa.graphBinFile = graphBinFile;
    ssa.seedFile = seedFile;
    long startSSA = time(NULL);
    ssa.getSolution(&sol, &reSSA);
    long timeSSA = time(NULL) - startSSA;
    double costSSA = ssa.calculateCost(sol);
    cout << "DSSA: " << reSSA << endl;
    cout << "DSSA Cost: " << costSSA << endl;
    cout << "DSSA Time: " << timeSSA << endl;

    writefile << Constant::K << "," << Constant::COMMUNITY_POPULATION
              << "," << reGIA << "," << costGIA << "," << timeGIA
              << "," << reHD << "," << costHD << "," << timeHD
              << "," << remaf << "," << costMaf << "," << timeMaf
              << "," << reubg << "," << costUbg << "," << timeUbg
              << "," << reSSA << "," << costSSA << "," << timeSSA << endl;

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
    // cout << "HB: " << reHb << endl;
}

void runExperiment(string input, string inputCommunity, int min, int max, int step, bool isScalable = true,
                   bool isBoundedThres = false, bool isLargeFile = false, bool changeK = false, bool isDirected = true,
                   bool isWeighted = false, bool isCommMM = true) {
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
                         + (isCommMM ? "_isMM" : "_isClauset") + ".csv";
    writefile.open(outfilename);
    if (writefile.is_open()) {
        if (Constant::IS_BOUNDED_THRESHOLD && !isLargeFile)
            writefile << "k \t Pop \t maf \t ubg-ratio \t grd \t bt \t hb \t ssa" << endl;
        else
            writefile << "k,Pop,gia,gia-cost,gia-time,"
                      << "hd,hd-cost,hd-time,"
                      << "maf,maf-cost,maf-time,"
                      << "ubg,ubg-cost,ubg-time,"
                      << "dssa,dssa-cost,dssa-time" << endl;
        // writefile << "k \t Pop \t maf \t ubg-ratio \t grd \t hb \t ssa" << endl;

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
    omp_set_num_threads(Constant::NUM_THREAD);

    vector<string> graphs{"facebook"};
    // vector<string> graphs{"facebook", "wiki", "epinions", "dblp", "pokec"};
    for (int i = 0; i < graphs.size(); ++i) {
        string graph = graphs[i];
        graphBinFile = "D:/DeadCommunity/data/" + graph + "SSA.bin";
        seedFile = "D:/DeadCommunity/data/" + graph + ".seeds";
        string input = graph + "E.txt";
        string inputCommunity = graph + "Comm.txt";
        runExperiment(input, inputCommunity, 4, 10, 2, true, false, false, false, true, false);
    }

    // Generation
    // for (int i = 0; i < graphs.size(); ++i) {
    //     string dir = "../data/";
    //     string graph = dir + graphs[i];
    //     string txt = graph + ".txt";
    //     string etxt = graph + "E.txt";
    //     string comm = graph + "Comm.txt";
    //     string adj = graph + ".adj";
    //     g->generateEdgeWeight(txt, etxt);
    //     g->readSocialGraph(etxt, true);
    //     g->formCommunityModularity(comm, adj, false);
    // }

    delete g;
    return 0;
}