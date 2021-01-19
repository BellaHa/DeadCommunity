#include "EIG.h"

EIG::EIG(SocialGraph *g) : Algorithm(g) {

}

EIG::~EIG() {
}

/*Test for speed up greedy*/
double EIG::getDeterministicSolution(vector<int> *sol) {
    sol->clear();
    vector<int> *listNodeIds = g->getListNodeIds();

    IloEnv env;
    try {
        IloModel model(env);

        // Variable: a node is chose as a seed
        IloNumVarArray x(env, g->getNumberOfNodes());
        for (int i = 0; i < g->getNumberOfNodes(); ++i) {
            x[i] = IloNumVar(env, 0, 1, ILOBOOL);
        }

        for (int i = 0; i < dcrSet.size(); ++i) {
            DCRgraph *dcr = dcrSet[i];
            map<int, vector<int>> *mapReachable = dcr->getMapReachable();
            vector<int> *commNodeIds = dcr->getCommunityNodeIds();
            IloNumVarArray y(env, mapReachable->size(), 0, 1, ILOBOOL);
            IloExpr deadNodes(env);
            for (int j = 0; j < commNodeIds->size(); ++j) {
                int nodeId = commNodeIds->at(j);
                vector<int> reachable = mapReachable->at(nodeId);
                IloExpr reachableSeed(env);
                for (int k = 0; k < reachable.size(); ++k) {
                    reachableSeed += x[mapNodeIdx[reachable[k]]];
                }
                model.add(IloIfThen(env, reachableSeed >= 1, y[j] == 1));
                model.add(IloIfThen(env, reachableSeed < 1, y[j] == 0));
                // model.add(reachableSeed >= y[j]);
                deadNodes += y[j];
                reachableSeed.end();
            }
            model.add(deadNodes >= dcr->getThreshold());
            deadNodes.end();
        }

        IloExpr obj(env);
        for (int i = 0; i < g->getNumberOfNodes(); ++i) {
            obj += x[i];
        }
        model.add(IloMinimize(env, obj));
        obj.end();

        IloCplex cplex(model);
        cplex.setParam(IloCplex::Threads, Constant::NUM_THREAD);
        cplex.setParam(IloCplex::TiLim, Constant::EIG_TIME);
        cplex.solve();

        for (int i = 0; i < g->getNumberOfNodes(); ++i) {
            if (cplex.getValue(x[i]) >= .5) {
                sol->push_back(listNodeIds->at(i));
            }
        }
    } catch (IloException &e) {
        cerr << "Conver exception caught: " << e << endl;
    } catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();

    return (double) (Constant::IS_WEIGHTED ? g->getNumberOfNodes() : g->getNumberOfCommunities());
}

double EIG::getDeterministicSolutionMig(vector<int> *sol) {
    sol->clear();
    vector<int> *listNodeIds = g->getListNodeIds();

    IloEnv env;
    try {
        IloModel model(env);

        // Variable: a node is chose as a seed
        IloNumVarArray x(env, g->getNumberOfNodes());
        for (int i = 0; i < g->getNumberOfNodes(); ++i) {
            x[i] = IloNumVar(env, 0, 1, ILOBOOL);
        }

        for (int i = 0; i < dcrSet.size(); ++i) {
            DCRgraph *dcr = dcrSet[i];
            map<int, vector<int>> *mapReachable = dcr->getMapReachable();
            vector<int> *commNodeIds = dcr->getCommunityNodeIds();
            IloNumVarArray y(env, mapReachable->size(), 0, 1, ILOBOOL);
            IloExpr deadNodes(env);
            for (int j = 0; j < commNodeIds->size(); ++j) {
                int nodeId = commNodeIds->at(j);
                vector<int> reachable = mapReachable->at(nodeId);
                IloExpr reachableSeed(env);
                for (int k = 0; k < reachable.size(); ++k) {
                    reachableSeed += x[mapNodeIdx[reachable[k]]];
                }
                model.add(IloIfThen(env, reachableSeed >= 1, y[j] == 1));
                model.add(IloIfThen(env, reachableSeed < 1, y[j] == 0));
                deadNodes += y[j] * g->mapNodeBenefit[nodeId];
                reachableSeed.end();
            }
            model.add(deadNodes >= dcr->thresholdB);
            deadNodes.end();
        }

        IloExpr obj(env);
        for (int i = 0; i < g->getNumberOfNodes(); ++i) {
            obj += x[i] * g->mapNodeCost[listNodeIds->at(i)];
        }
        model.add(IloMinimize(env, obj));
        obj.end();

        IloCplex cplex(model);
        cplex.setParam(IloCplex::Threads, Constant::NUM_THREAD);
        cplex.setParam(IloCplex::TiLim, Constant::EIG_TIME);
        cplex.solve();

        for (int i = 0; i < g->getNumberOfNodes(); ++i) {
            if (cplex.getValue(x[i]) >= .5) {
                sol->push_back(listNodeIds->at(i));
            }
        }
    } catch (IloException &e) {
        cerr << "Conver exception caught: " << e << endl;
    } catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();

    return (double) (Constant::IS_WEIGHTED ? g->getNumberOfNodes() : g->getNumberOfCommunities());
}

/*official running*/
// double EIG::getDeterministicSolution(vector<int> *sol) {
//     sol->clear();
//     vector<int> nodeIds(*(g->getListNodeIds()));
//     currentLive.clear();
//     for (int i = 0; i < dcrSet.size(); i++) {
//         DCRgraph *dcr = dcrSet[i];
//         vector<int> *commNodeIds = dcr->getCommunityNodeIds();
//         currentLive.push_back(vector<int>(*commNodeIds));
//     }
//
//
//     double gain = 0;
//     while (sol->size() < Constant::K) {
//         int maxIndex = 0;
//         double maxGain = 0;
//
// #pragma omp parallel for
//         for (int i = 0; i < nodeIds.size(); i++) {
//             int u = nodeIds[i];
//             double marginalGain = getMarginalGain(u, sol);
//
// #pragma omp critical
//             {
//                 if (marginalGain > maxGain) {
//                     maxIndex = i;
//                     maxGain = marginalGain;
//                 }
//             }
//         }
//
//         if (maxGain > 0) {
//             sol->push_back(nodeIds[maxIndex]);
//             gain += maxGain;
//
//             // update current live
// #pragma omp parallel for
//             for (int i = 0; i < dcrSet.size(); i++) {
//                 dcrSet[i]->getCurrentLiveAfterAddingNode(nodeIds[maxIndex], &(currentLive[i]));
//             }
//             nodeIds.erase(nodeIds.begin() + maxIndex);
//         } else break;
//     }
//     return gain * g->getNumberOfCommunities() / dcrSet.size();
// }

/*
double EIG::estimateInf(vector<int>* sol, int noDcr) {
	double re = 0.0;

	#pragma omp parallel for
	for (int i = 0; i < noDcr; i++) {
		DCRgraph * g = dcrSet[i];
		double fr = g->fractionalInf(sol);
		#pragma omp critical
		{
			re += 1.0;
		}
	}

	return re * g->getNumberOfCommunities() / dcrSet.size();
}
*/

double EIG::getSolution(vector<int> *sol, double *est) {
    sol->clear();
    initiateMig();
    omp_set_num_threads(Constant::NUM_THREAD);
    generateDCRgraphs((int) n1);
    double epsilon = Constant::EPSILON;
    double K = (double) g->getNumberOfCommunities();

    double re = 0.;
    for (int i = 0; i < iMax; ++i) {
        re = getDeterministicSolution(sol);
        *est = estimateInf(sol);
        if (*est >= (K - epsilon * K) || i == iMax - 1) {
            break;
        } else {
            generateDCRgraphs(dcrSet.size());
        }
    }
    clear();
    return *est / re;

    // while (dcrSet.size() < rMax) {
    //     double re = getDeterministicSolution(sol);
    //     int tmp = dcrSet.size();
    //     *est = estimateInf(sol);
    //     if (dcrSet.size() * (*est) / (Constant::IS_WEIGHTED ? g->getNumberOfNodes() : g->getNumberOfCommunities()) >=
    //         D) {
    //         double re2 = estimate(sol, e2, Constant::DELTA / 3, dcrSet.size());
    //         cout << re << " " << re2 << " " << time(NULL) << endl;
    //         if (re <= (1 + e1) * re2)
    //             return (*est) / re;
    //     }
    //
    //     int tmp2 = dcrSet.size();
    //     generateDCRgraphs(2 * tmp - tmp2);
    //
    // }
    // double re = getDeterministicSolution(sol);
    // *est = estimateInf(sol);
    // clear();
    // return (*est) / re;
}

double EIG::getSolutionMig(vector<int> *sol, double *est) {
    sol->clear();
    initiateMig();
    omp_set_num_threads(Constant::NUM_THREAD);
    generateDCRgraphsMig((int) n1);
    double epsilon = Constant::EPSILON;
    double K = (double) g->getNumberOfCommunities();

    double re = 0.;
    for (int i = 0; i < iMax; ++i) {
        re = getDeterministicSolutionMig(sol);
        *est = estimateInfMig(sol);
        if (*est >= (K - epsilon * K) || i == iMax - 1) {
            break;
        } else {
            generateDCRgraphsMig(dcrSet.size());
        }
    }
    clearMig();
    return *est / re;

    // while (dcrSet.size() < rMax) {
    //     double re = getDeterministicSolution(sol);
    //     int tmp = dcrSet.size();
    //     *est = estimateInf(sol);
    //     if (dcrSet.size() * (*est) / (Constant::IS_WEIGHTED ? g->getNumberOfNodes() : g->getNumberOfCommunities()) >=
    //         D) {
    //         double re2 = estimate(sol, e2, Constant::DELTA / 3, dcrSet.size());
    //         cout << re << " " << re2 << " " << time(NULL) << endl;
    //         if (re <= (1 + e1) * re2)
    //             return (*est) / re;
    //     }
    //
    //     int tmp2 = dcrSet.size();
    //     generateDCRgraphs(2 * tmp - tmp2);
    //
    // }
    // double re = getDeterministicSolution(sol);
    // *est = estimateInf(sol);
    // clear();
    // return (*est) / re;
}

double EIG::getSolution2Step(vector<int> *sol, double *est) {
    sol->clear();
    initiate();
    omp_set_num_threads(Constant::NUM_THREAD);
    generateDCRgraphs((int) rMax);
    double re = getDeterministicSolution(sol);
    *est = estimateInf(sol);
    clear();
    return (*est) / re;
}

double EIG::estimate(vector<int> *sol, double epsilon, double delta, int tMax) {
    double lamda = 0.72;
    double tmp = 4 * lamda * log(2 / delta) / (epsilon * epsilon);
    double lambda = 1 + (1 + epsilon) * tmp;

    int T = 0;
    double inf = 0.0;

#pragma omp parallel for
    for (int i = 0; i < tMax; i++) {
        DCRgraph *g = gen.generateDCRgraph();
        double fr = g->fractionalInf(sol);

#pragma omp critical
        {
            dcrSet.push_back(g);
            g->updateInitalGain(&intialGain, &initialDead);
            if (tMax > 0) {
                T++;
                inf += fr;
                if (inf >= lambda) {
                    tMax = -1;
                }
            }
        }
    }

    return (tMax == -1 ? lambda * (Constant::IS_WEIGHTED ? g->getNumberOfNodes() : g->getNumberOfCommunities()) / T
                       : -1);
}

double EIG::getMarginalGain(int nodeId, vector<int> *sol) {
    double re = 0.0;
    for (int i = 0; i < dcrSet.size(); i++) {
        DCRgraph *dcr = dcrSet[i];
        int gain = dcr->getMarginalGain(nodeId, &(currentLive[i]));
        re += ((double) gain) / dcr->getThreshold();
    }
    return re;
}

void EIG::generateDCRgraphs(int number) {
#pragma omp parallel for
    for (int i = 0; i < number; i++) {
        DCRgraph *dcr = gen.generateDCRgraphMig();
#pragma omp critical
        {
            dcrSet.push_back(dcr);
            if (!isMaf)
                dcr->updateInitalGain(&intialGain, &initialDead);
            else
                dcr->updateInitalGain(&intialGain, &initialDead, &initialOtherCommunityGain);
        }
        //cout << i << endl;
    }
    //cout << "done generating samples" << endl;
}

double EIG::estimateInf(vector<int> *sol) {
    double K = (double) g->getNumberOfCommunities();
    double T = (double) dcrSet.size();
    double nb1 = K - ((K * c) / (3. * T));

    double Xsol = 0.;
#pragma omp parallel for
    for (int i = 0; i < dcrSet.size(); i++) {
        bool kill = dcrSet[i]->isKill(sol);

        if (kill) {
#pragma omp critical
            {
                Xsol += 1.0;
            }
        }
    }

    double eSigma = (K / T) * Xsol;
    double nb2 = K + (K / T) * ((2. * c / 3) - sqrt((4. * c * c / 9.) + (2. * T * c * (eSigma / K))));
    return min(nb1, nb2);
}

