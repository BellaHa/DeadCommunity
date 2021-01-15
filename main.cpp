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

using namespace std;

#pragma warning(disable : 4996)

using namespace std;

SocialGraph *g;
ofstream writefile;
string graphBinFile;
string seedFile;

void printResult(bool isScalable, bool isLargeFile) {
    vector<int> sol;
    sol.clear();

    GIA gia(g);
    long startGIA = time(NULL);
    double reGIA = 0;
    gia.getSolution(&sol, &reGIA);
    long timeGIA = time(NULL) - startGIA;
    cout << "GIA: " << reGIA << endl;
    cout << "GIA Time: " << timeGIA << endl;
    cout << sol.size() << endl;

    HighDegree hd(g);
    long startHD = time(NULL);
    double reHD = 0;
    hd.getSolution(&sol, &reHD);
    long timeHD = time(NULL) - startHD;
    cout << "HD: " << reHD << endl;
    cout << "HD Time: " << timeHD << endl;
    cout << sol.size() << endl;

    GreedySolution maf(g);
    long startMaf = time(NULL);
    double remaf = 0;
    maf.getSolution(&sol, &remaf);
    long timeMaf = time(NULL) - startMaf;
    cout << "MAF: " << remaf << endl;
    cout << "MAF Time: " << timeMaf << endl;
    cout << sol.size() << endl;

    SandwichSolution ubg(g);
    long startUbg = time(NULL);
    double reubg = 0;
    ubg.getSolution(&sol, &reubg);
    long timeUbg = time(NULL) - startUbg;
    cout << "UBG: " << reubg << endl;
    cout << "UBG Time: " << timeUbg << endl;
    cout << sol.size() << endl;

    SSA ssa(g);
    double reSSA = 0;
    long startSSA = time(NULL);
    ssa.graphBinFile = graphBinFile;
    ssa.seedFile = seedFile;
    ssa.getSolution(&sol, &reSSA);
    long timeSSA = time(NULL) - startSSA;
    cout << "SSA: " << reSSA << endl;
    cout << "SSA Time: " << timeSSA << endl;
    cout << sol.size() << endl;

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

    writefile << Constant::K << "," << Constant::COMMUNITY_POPULATION << "," << reGIA << "," << timeGIA
              << "," << reHD << "," << timeHD
              << "," << remaf << "," << timeMaf
              << "," << reubg << "," << timeUbg
              << "," << reSSA << "," << timeSSA << endl;
    // //<< "\t" << reGrd << " " << timeGrd
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
            writefile << "k,Pop,gia,gia-time,hd,hd-time,maf,maf-time,ubg,ubg-time,ssa,ssa-time" << endl;
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

    // Generate community
    for (int i = 0; i < graphs.size(); ++i) {
        string dir = "../data/";
        string graph = dir + graphs[i];
        string txt = graph + ".txt";
        string etxt = graph + "E.txt";
        string comm = graph + "Comm.txt";
        string adj = graph + ".adj";
        g->generateEdgeWeight(txt, etxt);
        g->readSocialGraph(etxt, true);
        g->formCommunityModularity(comm, adj, false);
    }

    delete g;
    return 0;
}