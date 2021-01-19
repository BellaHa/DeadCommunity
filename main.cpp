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
#include "StopWatch.h"
#include "BoundedThres.h"
#include "CompareGreedy.h"
#include "HighInfluence.h"
#include "HighTouch.h"
#include "HighBenefit.h"
#include "Constant.h"
#include "EIG.h"

using namespace std;

#pragma warning(disable : 4996)

SocialGraph *g;
ofstream writefile;
string graphBinFile;
string seedFile;

void printResult(bool isScalable, bool isLargeFile) {
    vector<int> sol;
    sol.clear();
    StopWatch sw;

    double reGIA, reEIG, reHD, reMAF, reUBG, reSSA;
    double timeGIA, timeEIG, timeHD, timeMAF, timeUBG, timeSSA;
    double costGIA, costEIG, costHD, costMAF, costUBG, costSSA;

    cout << "GIA..." << endl;
    GIA gia(g);
    sw.start();
    if (Constant::GCS)
        gia.getSolutionMig(&sol, &reGIA);
    else
        gia.getSolution(&sol, &reGIA);
    sw.stop();
    timeGIA = sw.getSeconds();
    costGIA = gia.calculateCost(sol);
    cout << "GIA: " << reGIA << endl;
    cout << "GIA Cost: " << costGIA << endl;
    cout << "GIA Time: " << timeGIA << endl << endl;

    cout << "EIG..." << endl;
    EIG eig(g);
    sw.start();
    if (Constant::GCS)
        eig.getSolutionMig(&sol, &reEIG);
    else
        eig.getSolution(&sol, &reEIG);
    sw.stop();
    timeEIG = sw.getSeconds();
    costEIG = eig.calculateCost(sol);
    cout << "EIG: " << reEIG << endl;
    cout << "EIG Cost: " << costEIG << endl;
    cout << "EIG Time: " << timeEIG << endl << endl;

    cout << "HD..." << endl;
    HighDegree hd(g);
    sw.start();
    if (Constant::GCS)
        hd.getSolutionMig(&sol, &reHD);
    else
        hd.getSolution(&sol, &reHD);
    sw.stop();
    timeHD = sw.getSeconds();
    costHD = hd.calculateCost(sol);
    cout << "HD: " << reHD << endl;
    cout << "HD Cost: " << costHD << endl;
    cout << "HD Time: " << timeHD << endl << endl;
    // return;

    if (!Constant::GCS) {
        cout << "MAF..." << endl;
        GreedySolution maf(g);
        sw.start();
        maf.getSolutionBS(&sol, &reMAF, 1, g->getNumberOfNodes());
        sw.stop();
        timeMAF = sw.getSeconds();
        costMAF = maf.calculateCost(sol);
        cout << "MAF: " << reMAF << endl;
        cout << "MAF Cost: " << costMAF << endl;
        cout << "MAF Time: " << timeMAF << endl << endl;

        cout << "UBG..." << endl;
        SandwichSolution ubg(g);
        sw.start();
        ubg.getSolutionFastBS(&sol, &reUBG, 1, g->getNumberOfNodes());
        sw.stop();
        timeUBG = sw.getSeconds();
        costUBG = ubg.calculateCost(sol);
        cout << "UBG: " << reUBG << endl;
        cout << "UBG Cost: " << costUBG << endl;
        cout << "UBG Time: " << timeUBG << endl << endl;

        cout << "DSSA..." << endl;
        SSA ssa(g);
        ssa.graphBinFile = graphBinFile;
        ssa.seedFile = seedFile;
        ssa.bsTime = 0;
        ssa.getSolutionBS(&sol, &reSSA, 1, g->getNumberOfNodes());
        timeSSA = ssa.bsTime;
        costSSA = ssa.calculateCost(sol);
        cout << "DSSA: " << reSSA << endl;
        cout << "DSSA Cost: " << costSSA << endl;
        cout << "DSSA Time: " << timeSSA << endl << endl;
    }

    writefile << g->getNumberOfCommunities() << "," << Constant::COMMUNITY_POPULATION
              << "," << reGIA << "," << costGIA << "," << timeGIA
              << "," << reHD << "," << costHD << "," << timeHD;
    if (!Constant::GCS) {
        writefile << "," << reMAF << "," << costMAF << "," << timeMAF
                  << "," << reUBG << "," << costUBG << "," << timeUBG
                  << "," << reSSA << "," << costSSA << "," << timeSSA;
    }
    writefile << endl;

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
        // g->readSocialGraph("../data/" + input);
        g->readSocialGraphBin("../data/" + input);
        g->readCommunityFile("../data/" + inputCommunity, isCommMM);
    }

    string outfilename = "../results/" + input + "_result_"
                         + (Constant::GCS ? "gcs_" : "ucs_")
                         + ".csv";
    ifstream rIn(outfilename);
    if (rIn.peek() == EOF) {
        rIn.close();
        ofstream out(outfilename);
        out << "k,Pop,gia,gia-cost,gia-time,"
            << "hd,hd-cost,hd-time,"
            << "maf,maf-cost,maf-time,"
            << "ubg,ubg-cost,ubg-time,"
            << "dssa,dssa-cost,dssa-time" << endl;
        out.close();
    }
    rIn.close();
    writefile.open(outfilename, ofstream::out | ofstream::app);

    if (writefile.is_open()) {
        // if (Constant::IS_BOUNDED_THRESHOLD && !isLargeFile)
        //     writefile << "k \t Pop \t maf \t ubg-ratio \t grd \t bt \t hb \t ssa" << endl;
        // else
        //     writefile << "k,Pop,gia,gia-cost,gia-time,"
        //               << "hd,hd-cost,hd-time,"
        //               << "maf,maf-cost,maf-time,"
        //               << "ubg,ubg-cost,ubg-time,"
        //               << "dssa,dssa-cost,dssa-time" << endl;
        // writefile << "k \t Pop \t maf \t ubg-ratio \t grd \t hb \t ssa" << endl;
        if (changeK) {
            if (!isDirected)
                g->formCommunitiesFromActualCommunities();
            for (Constant::K = min; Constant::K <= 20; Constant::K += step) {
                printResult(isScalable, isLargeFile);
            }
        } else {
            // for (Constant::COMMUNITY_POPULATION = min;
            //      Constant::COMMUNITY_POPULATION <= max; Constant::COMMUNITY_POPULATION += step) {
            for (Constant::NUMBER_OF_COMMS = min; Constant::NUMBER_OF_COMMS <= max; Constant::NUMBER_OF_COMMS += step) {
                g->formCommunitiesFromActualCommunities();
                printResult(isScalable, isLargeFile);
            }
            // }
        }
        writefile.close();
    }
}

int main() {
    g = new SocialGraph();
    omp_set_num_threads(Constant::NUM_THREAD);

    vector<string> graphs{"facebook"};
    vector<int> min{10, 100, 200, 500, 1000};
    vector<int> max{90, 900, 2000, 5000, 10000};
    vector<int> step{10, 100, 200, 500, 1000};
    // vector<string> graphs{"facebook", "wiki", "epinions", "dblp", "pokec"};
    for (int i = 0; i < graphs.size(); ++i) {
        string graph = graphs[i];
        graphBinFile = "D:/DeadCommunity/data/" + graph + "SSA.bin";
        seedFile = "D:/DeadCommunity/data/" + graph + ".seeds";
        // string inputTxt = graph + "E.txt";
        string inputBin = graph + ".bin";
        string inputCommunity = graph + "Comm.txt";
        runExperiment(inputBin, inputCommunity, min[i], max[i], step[i], true, false, false, false, true, false);
    }

    // Generation
    // for (int i = 0; i < graphs.size(); ++i) {
    //     string dir = "../data/";
    //     string graph = dir + graphs[i];
    //     cout << graph << endl;
    //     string txt = graph + ".txt";
    //     string bin = graph + ".bin";
    //     string etxt = graph + "E.txt";
    //     string comm = graph + "Comm.txt";
    //     string adj = graph + ".adj";
    //     string ssaFile = graph + "SSA.txt";
    //     // g->generateEdgeWeightBinFile(txt, etxt, bin);
    //     g->readSocialGraphBin(bin, true);
    //     g->generateFileIM(ssaFile);
    //     g->formCommunityModularity(comm, adj, false);
    // }

    delete g;
    return 0;
}