#include "DCRgenerator.h"
#include "Constant.h"
#include <algorithm>
#include <iostream>


DCRgenerator::DCRgenerator() {
    commonInstance = Common::getInstance();
}

DCRgenerator::DCRgenerator(SocialGraph *g) {
    this->g = g;
    commonInstance = Common::getInstance();
}

DCRgenerator::~DCRgenerator() {
}

void DCRgenerator::setSocialGraph(SocialGraph *g) {
    this->g = g;
}

// LT model
DCRgraph *DCRgenerator::generateDCRgraphLT() {
    int commId = g->randomSelectCommunity();
    vector<int> *nodeIds = g->getNodesOfCommunity(commId);
    int threshold = Constant::IS_BOUNDED_THRESHOLD ? 2 : (int) (Constant::PERCENTAGE_THRESHOLD * nodeIds->size());
    DCRgraph *dcr = new DCRgraph(commId, threshold > 0 ? threshold : 1, nodeIds);
    map<int, bool> ck;
    map<int, int> p;

    vector<int> queue;
    for (int i = 0; i < nodeIds->size(); i++) {
        queue.push_back(nodeIds->at(i));
        ck[nodeIds->at(i)] = true;
    }

    while (!queue.empty()) {
        int u = queue[0];
        queue.erase(queue.begin());
        vector<pair<int, double>> *incommingNeighbors = g->getIncommingNeighbors(u);
        if (incommingNeighbors != nullptr && incommingNeighbors->size() > 0) {
            int select = commonInstance->randomInThread() % (incommingNeighbors->size() + 1);
            if (select < incommingNeighbors->size()) {
                int parent = incommingNeighbors->at(select).first;
                p[u] = parent;
                if (ck.find(parent) == ck.end()) {
                    queue.push_back(parent);
                    ck[parent] = true;
                }
            }
        }
    }

    // reverse parent to get reachable set
    for (int i = 0; i < nodeIds->size(); i++) {
        int nodeId = nodeIds->at(i);
        vector<int> reachable;
        int trace = nodeId;
        while (find(reachable.begin(), reachable.end(), trace) == reachable.end()) {
            reachable.push_back(trace);
            if (p.find(trace) != p.end())
                trace = p[trace];
            else break;
        }
        dcr->addReachable(nodeId, &reachable);
    }

    return dcr;
}

DCRgraph *DCRgenerator::generateDCRgraph() {
    if (Constant::MODEL)
        return generateDCRgraphLT();
    else
        return generateDCRgraphIC();
}

// IC model
DCRgraph *DCRgenerator::generateDCRgraphIC() {
    int commId = g->randomSelectCommunity();
    vector<int> *nodeIds = g->getNodesOfCommunity(commId);
    int threshold = Constant::IS_BOUNDED_THRESHOLD ? 2 : (int) (Constant::PERCENTAGE_THRESHOLD * nodeIds->size());
    DCRgraph *dcr = new DCRgraph(commId, threshold > 0 ? threshold : 1, nodeIds);
    map<int, vector<int>> mapNeighbors;
    map<int, map<int, bool>> st;
    map<int, bool> ck;

    vector<int> queue;
    for (int i = 0; i < nodeIds->size(); i++) {
        queue.push_back(nodeIds->at(i));
        ck[nodeIds->at(i)] = true;
    }

    while (!queue.empty()) {
        int u = queue[0];
        queue.erase(queue.begin());
        vector<pair<int, double>> *incommingNeighbors = g->getIncommingNeighbors(u);
        if (incommingNeighbors != nullptr && incommingNeighbors->size() > 0) {
            for (int i = 0; i < incommingNeighbors->size(); i++) {
                pair<int, double> tmp = incommingNeighbors->at(i);
                int v = tmp.first;
                double w = tmp.second;

                map<int, bool> *mapP = &(st[u]);

                if (mapP->find(v) == mapP->end()) {
                    double coin = ((double) (commonInstance->randomInThread() % 1000)) / 1000;
                    (*mapP)[v] = (coin <= w);
                }

                if ((*mapP)[v] && (ck.find(v) == ck.end())) {
                    queue.push_back(v);
                    ck[v] = true;
                    mapNeighbors[u].push_back(v);
                }
            }
        }

    }

    // reverse dfs to get reachable set
    for (int i = 0; i < nodeIds->size(); i++) {
        int nodeId = nodeIds->at(i);
        vector<int> reachable;

        map<int, bool> ck;
        vector<int> queue;
        queue.push_back(nodeId);
        ck[nodeId] = true;
        while (!queue.empty()) {
            int u = queue[0];
            queue.erase(queue.begin());
            reachable.push_back(u);
            vector<int> p = mapNeighbors[u];
            for (int j = 0; j < p.size(); j++) {
                if (ck.find(p[j]) == ck.end()) {
                    queue.push_back(p[j]);
                    ck[p[j]] = true;
                }
            }
        }

        //dfs(nodeId, &reachable, &mapNeighbors);
        dcr->addReachable(nodeId, &reachable);
    }

    return dcr;
}

void DCRgenerator::dfs(int u, vector<int> *reachable, map<int, vector<int>> *mapNeighbors) {
    reachable->push_back(u);
    vector<int> neighbor = (*mapNeighbors)[u];
    for (int i = 0; i < neighbor.size(); i++) {
        int v = neighbor[i];
        if (find(reachable->begin(), reachable->end(), v) == reachable->end()) {
            dfs(v, reachable, mapNeighbors);
        }
    }
}

DCRgraph *DCRgenerator::generateDCRgraphICMig() {
    int commId = g->randomSelectCommunity();
    vector<int> *nodeIds = g->getNodesOfCommunity(commId);
    int threshold = Constant::IS_BOUNDED_THRESHOLD ? 2 : (int) (Constant::PERCENTAGE_THRESHOLD * nodeIds->size());
    DCRgraph *dcr = new DCRgraph(commId, threshold > 0 ? threshold : 1, nodeIds);
    map<int, vector<int>> mapNeighbors;
    map<int, map<int, bool>> st;
    // map<int, bool> ck;

    for (int i = 0; i < nodeIds->size(); ++i) {
        int nodeId = nodeIds->at(i);
        vector<int> reachable{nodeId};
        vector<int> queue;
        queue.emplace_back(nodeId);
        map<int, bool> ck;
        ck[nodeId] = true;
        while (!queue.empty()) {
            int u = queue.at(0);
            queue.erase(queue.begin());
            vector<pair<int, double>> *incommingNeighbors = g->getIncommingNeighbors(u);
            if (incommingNeighbors != nullptr && incommingNeighbors->size() > 0) {
                for (int j = 0; j < incommingNeighbors->size(); ++j) {
                    pair<int, double> tmp = incommingNeighbors->at(j);
                    int v = tmp.first;
                    double w = tmp.second;
                    map<int, bool> *mapP = &(st[u]);

                    if (mapP->find(v) == mapP->end()) {
                        double coin = ((double) (commonInstance->randomInThread() % 1000)) / 1000;
                        (*mapP)[v] = (coin <= w);
                    }

                    if ((*mapP)[v] && ck.find(v) == ck.end()) {
                        queue.push_back(v);
                        ck[v] = true;
                        reachable.push_back(v);
                    }
                }
            }
        }
        dcr->addReachable(nodeId, &reachable);
    }
    return dcr;

    // vector<int> queue;
    // for (int i = 0; i < nodeIds->size(); i++) {
    //     queue.push_back(nodeIds->at(i));
    //     ck[nodeIds->at(i)] = true;
    // }
    //
    // while (!queue.empty()) {
    //     int u = queue[0];
    //     queue.erase(queue.begin());
    //     vector<pair<int, double>> *incommingNeighbors = g->getIncommingNeighbors(u);
    //     if (incommingNeighbors != nullptr && incommingNeighbors->size() > 0) {
    //         for (int i = 0; i < incommingNeighbors->size(); i++) {
    //             pair<int, double> tmp = incommingNeighbors->at(i);
    //             int v = tmp.first;
    //             double w = tmp.second;
    //
    //             map<int, bool> *mapP = &(st[u]);
    //
    //             if (mapP->find(v) == mapP->end()) {
    //                 double coin = ((double) (commonInstance->randomInThread() % 1000)) / 1000;
    //                 (*mapP)[v] = (coin <= w);
    //             }
    //
    //             if ((*mapP)[v] && (ck.find(v) == ck.end())) {
    //                 queue.push_back(v);
    //                 ck[v] = true;
    //                 mapNeighbors[u].push_back(v);
    //             }
    //         }
    //     }
    //
    // }
    //
    // // reverse dfs to get reachable set
    // for (int i = 0; i < nodeIds->size(); i++) {
    //     int nodeId = nodeIds->at(i);
    //     vector<int> reachable;
    //
    //     map<int, bool> ck;
    //     vector<int> queue;
    //     queue.push_back(nodeId);
    //     ck[nodeId] = true;
    //     while (!queue.empty()) {
    //         int u = queue[0];
    //         queue.erase(queue.begin());
    //         reachable.push_back(u);
    //         vector<int> p = mapNeighbors[u];
    //         for (int j = 0; j < p.size(); j++) {
    //             if (ck.find(p[j]) == ck.end()) {
    //                 queue.push_back(p[j]);
    //                 ck[p[j]] = true;
    //             }
    //         }
    //     }
    //
    //     //dfs(nodeId, &reachable, &mapNeighbors);
    //     dcr->addReachable(nodeId, &reachable);
    // }
    //
    // return dcr;
}

DCRgraph *DCRgenerator::generateDCRgraphICMigB() {
    int commId = g->randomSelectCommunity();
    vector<int> *nodeIds = g->getNodesOfCommunity(commId);
    int threshold = Constant::IS_BOUNDED_THRESHOLD ? 2 : (int) (Constant::PERCENTAGE_THRESHOLD * nodeIds->size());
    double commBenefit = g->mapCommBenefit[commId];
    double thresholdB = (int) (Constant::PERCENTAGE_THRESHOLD * commBenefit);
    thresholdB = thresholdB > 0 ? thresholdB : 1;
    // cout << threshold << endl;
    // cout << thresholdB << endl;
    DCRgraph *dcr = new DCRgraph(g, commId, threshold > 0 ? threshold : 1, commBenefit, thresholdB, nodeIds);
    map<int, vector<int>> mapNeighbors;
    map<int, map<int, bool>> st;
    // map<int, bool> ck;

    for (int i = 0; i < nodeIds->size(); ++i) {
        int nodeId = nodeIds->at(i);
        vector<int> reachable{nodeId};
        vector<int> queue;
        queue.emplace_back(nodeId);
        map<int, bool> ck;
        ck[nodeId] = true;
        while (!queue.empty()) {
            int u = queue.at(0);
            queue.erase(queue.begin());
            vector<pair<int, double>> *incommingNeighbors = g->getIncommingNeighbors(u);
            if (incommingNeighbors != nullptr && incommingNeighbors->size() > 0) {
                for (int j = 0; j < incommingNeighbors->size(); ++j) {
                    pair<int, double> tmp = incommingNeighbors->at(j);
                    int v = tmp.first;
                    double w = tmp.second;
                    map<int, bool> *mapP = &(st[u]);

                    if (mapP->find(v) == mapP->end()) {
                        double coin = ((double) (commonInstance->randomInThread() % 1000)) / 1000;
                        (*mapP)[v] = (coin <= w);
                    }

                    if ((*mapP)[v] && ck.find(v) == ck.end()) {
                        queue.push_back(v);
                        ck[v] = true;
                        reachable.push_back(v);
                    }
                }
            }
        }
        dcr->addReachable(nodeId, &reachable);
    }
    return dcr;

    // vector<int> queue;
    // for (int i = 0; i < nodeIds->size(); i++) {
    //     queue.push_back(nodeIds->at(i));
    //     ck[nodeIds->at(i)] = true;
    // }
    //
    // while (!queue.empty()) {
    //     int u = queue[0];
    //     queue.erase(queue.begin());
    //     vector<pair<int, double>> *incommingNeighbors = g->getIncommingNeighbors(u);
    //     if (incommingNeighbors != nullptr && incommingNeighbors->size() > 0) {
    //         for (int i = 0; i < incommingNeighbors->size(); i++) {
    //             pair<int, double> tmp = incommingNeighbors->at(i);
    //             int v = tmp.first;
    //             double w = tmp.second;
    //
    //             map<int, bool> *mapP = &(st[u]);
    //
    //             if (mapP->find(v) == mapP->end()) {
    //                 double coin = ((double) (commonInstance->randomInThread() % 1000)) / 1000;
    //                 (*mapP)[v] = (coin <= w);
    //             }
    //
    //             if ((*mapP)[v] && (ck.find(v) == ck.end())) {
    //                 queue.push_back(v);
    //                 ck[v] = true;
    //                 mapNeighbors[u].push_back(v);
    //             }
    //         }
    //     }
    //
    // }
    //
    // // reverse dfs to get reachable set
    // for (int i = 0; i < nodeIds->size(); i++) {
    //     int nodeId = nodeIds->at(i);
    //     vector<int> reachable;
    //
    //     map<int, bool> ck;
    //     vector<int> queue;
    //     queue.push_back(nodeId);
    //     ck[nodeId] = true;
    //     while (!queue.empty()) {
    //         int u = queue[0];
    //         queue.erase(queue.begin());
    //         reachable.push_back(u);
    //         vector<int> p = mapNeighbors[u];
    //         for (int j = 0; j < p.size(); j++) {
    //             if (ck.find(p[j]) == ck.end()) {
    //                 queue.push_back(p[j]);
    //                 ck[p[j]] = true;
    //             }
    //         }
    //     }
    //
    //     //dfs(nodeId, &reachable, &mapNeighbors);
    //     dcr->addReachable(nodeId, &reachable);
    // }
    //
    // return dcr;
}

DCRgraph *DCRgenerator::generateDCRgraphMig() {
    if (Constant::MODEL)
        return generateDCRgraphLT();
    else
        return generateDCRgraphICMig();
}

DCRgraph *DCRgenerator::generateDCRgraphMigB() {
    if (Constant::MODEL)
        return generateDCRgraphLT();
    else
        return generateDCRgraphICMigB();
}

/*
DCRgraph * DCRgenerator::generateDCRgraph()
{
	int commId = g->randomSelectCommunity();
	vector<int>* nodeIds = g->getNodesOfCommunity(commId);
	DCRgraph * dcr = new DCRgraph(commId, (int)(Constant::PERCENTAGE_THRESHOLD * nodeIds->size()), nodeIds);
	map<int, vector<std::pair<int, double>>> mapIncommingNeighbor(*g->getMapIncommingNeighbors());
	
	// reverse dfs to get reachable set

	for (int i = 0; i < nodeIds->size(); i++) {
		int nodeId = nodeIds->at(i);
		vector<int> reachable;
		reverseDfs(nodeId, &reachable, &mapIncommingNeighbor);
		sort(reachable.begin(), reachable.end());
		dcr->addReachable(nodeId, &reachable);
	}

	return dcr;
}

void DCRgenerator::reverseDfs(int nodeId,
	vector<int> * reachable,
	map<int, vector<std::pair<int, double>>> * mapIncommingNeighbor)
{
	reachable->push_back(nodeId);
	vector<pair<int, double>> * incommingNeighbors = &(*mapIncommingNeighbor)[nodeId];
	for (int i = 0; i < incommingNeighbors->size(); i++) {
		int neighbor = (*incommingNeighbors)[i].first;
		if (find(reachable->begin(), reachable->end(), nodeId) == reachable->end()) {
			double w = (*incommingNeighbors)[i].second;
			if (w < 0.999999) {
				double coin = ((double)(rand() % 1000)) / 1000;
				if (coin <= w) {
					(*incommingNeighbors)[i] = pair<int, double>(neighbor, 1);
					reverseDfs(neighbor, reachable, mapIncommingNeighbor);
				}
			}
			else {
				reverseDfs(neighbor, reachable, mapIncommingNeighbor);
			}
		}
	}
}
*/