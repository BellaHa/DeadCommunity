#include "SocialGraph.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "Constant.h"
#include <string>
#include <algorithm>
#include <sstream>
#include <queue>

using namespace std;


SocialGraph::SocialGraph() {
    commonInstance = Common::getInstance();
}

SocialGraph::~SocialGraph() {
}

void SocialGraph::readSocialGraphFromFile(string file) {
    clear();
    srand(time(NULL));
    ifstream inputFile;
    inputFile.open(file);
    if (inputFile) {
        int nodeId = 0, commId = 0;
        inputFile >> nodeId >> commId;
        while (nodeId != -1) {
            if (listCommListNodeIds.size() > commId) {
                listCommListNodeIds[commId].push_back(nodeId);
            } else {
                vector<int> tmp;
                tmp.push_back(nodeId);
                listCommListNodeIds.push_back(tmp);
            }
            listNodeIds.push_back(nodeId);
            mapNodeId2CommId[nodeId] = commId;
            inputFile >> nodeId >> commId;
        }

        // setup hMax
        if (Constant::IS_BOUNDED_THRESHOLD)
            hMax = 2;
        else {
            for (int i = 0; i < listCommListNodeIds.size(); i++) {
                if (hMax < listCommListNodeIds[i].size() * Constant::PERCENTAGE_THRESHOLD)
                    hMax = (int) listCommListNodeIds[i].size() * Constant::PERCENTAGE_THRESHOLD;
            }
        }

        int startNode, endNode;
        double weight;
        while (inputFile >> startNode >> endNode >> weight) {
            if (mapIncommingNeighbors.find(endNode) != mapIncommingNeighbors.end()) {
                mapIncommingNeighbors[endNode].push_back(pair<int, double>(startNode, weight));
            } else {
                vector<pair<int, double>> tmp;
                tmp.push_back(pair<int, double>(startNode, weight));
                mapIncommingNeighbors.insert(pair<int, vector<pair<int, double>>>(endNode, tmp));
            }
            noOfEdges++;
        }
        inputFile.close();
    }
}

void SocialGraph::readSocialGraph(string inputFile, bool isDirected) {
    this->isDirected = isDirected;
    clear();
    ifstream file;
    file.open(inputFile);
    if (file) {
        int startNode, endNode;
        double weight;
        while (file >> startNode >> endNode >> weight) {
            if (mapIncommingNeighbors.find(endNode) != mapIncommingNeighbors.end()) {
                mapIncommingNeighbors[endNode].push_back(pair<int, double>(startNode, weight));
            } else {
                vector<pair<int, double>> tmp;
                tmp.emplace_back(startNode, weight);
                mapIncommingNeighbors.insert(pair<int, vector<pair<int, double>>>(endNode, tmp));
            }
            if (mapOutgoingNeighbors.find(startNode) != mapOutgoingNeighbors.end()) {
                mapOutgoingNeighbors[startNode].push_back(pair<int, double>(endNode, weight));
            } else {
                vector<pair<int, double>> tmp;
                tmp.emplace_back(endNode, weight);
                mapOutgoingNeighbors.insert(pair<int, vector<pair<int, double>>>(startNode, tmp));
            }
            if (find(listNodeIds.begin(), listNodeIds.end(), startNode) == listNodeIds.end())
                listNodeIds.push_back(startNode);

            if (find(listNodeIds.begin(), listNodeIds.end(), endNode) == listNodeIds.end())
                listNodeIds.push_back(endNode);
            noOfEdges++;
        }
        outgoingDegree.clear();
        for (int i = 0; i < getNumberOfNodes(); ++i) {
            outgoingDegree.emplace_back(mapOutgoingNeighbors[listNodeIds[i]].size());
        }
        file.close();
        // int startId, endId;
        // int count = 0;
        // while (file >> startId >> endId) {
        //     mapIncommingNeighbors[endId].push_back(pair<int, double>(startId, 0));
        //
        //     if (!isDirected)
        //         mapIncommingNeighbors[startId].push_back(pair<int, double>(endId, 0));
        //
        //     if (find(listNodeIds.begin(), listNodeIds.end(), startId) == listNodeIds.end())
        //         listNodeIds.push_back(startId);
        //
        //     if (find(listNodeIds.begin(), listNodeIds.end(), endId) == listNodeIds.end())
        //         listNodeIds.push_back(endId);
        //
        //     noOfEdges++;
        // }
    }
}

void extractFromCharArray(char *cArray, int sp, int ep) {
    std::string::size_type sz;
    int startId = stoi(cArray + sp, &sz);
    int endId = stoi(cArray + sz);
}

void SocialGraph::readSocialGraphFromLargeFile(string inputFile) {
    clear();
    ifstream is(inputFile);
    is.seekg(0, is.end);
    long bufSize = is.tellg();
    is.seekg(0, is.beg);
    int item = 0;

    char *buffer = new char[bufSize];

    is.read(buffer, bufSize);
    is.close();


    std::string::size_type sz = 0;
    long sp = 0;
    int startId, endId;
    bool isStart = true;

    map<int, vector<int>> mapTmp;
    map<int, bool> mapNode;

    while (sp < bufSize) {
        char c = buffer[sp];
        item = item * 10 + c - 48;
        sp++;
        if (sp == bufSize || (sp < bufSize && (buffer[sp] == '\t' || buffer[sp] == '\n'))) {
            sp++;

            if (isStart) {
                startId = item;
                isStart = false;
            } else {
                endId = item;
                isStart = true;
                mapTmp[endId].push_back(startId);
                noOfEdges++;
            }

            mapNode[item] = true;

            item = 0;
        }
    }

    for (map<int, bool>::iterator it = mapNode.begin(); it != mapNode.end(); ++it) {
        listNodeIds.push_back(it->first);
    }

    for (map<int, vector<int>>::iterator it = mapTmp.begin(); it != mapTmp.end(); ++it) {
        double w = 1.0 / it->second.size();
        vector<pair<int, double>> tmp;
        for (int i = 0; i < it->second.size(); ++i) {
            tmp.push_back(pair<int, double>(it->second[i], w));
        }
        mapIncommingNeighbors[it->first] = tmp;
    }

    // assign community

    // assign community
    int numberOfNodes = listNodeIds.size();
    int numberOfCommunities = ceil(((double) numberOfNodes) / Constant::COMMUNITY_POPULATION);

    vector<int> index;
    for (int i = 0; i < numberOfNodes; i++)
        index.push_back(i);


    std::random_shuffle(index.begin(), index.end());
    listCommListNodeIds = vector<vector<int>>(numberOfCommunities, vector<int>());
    for (int i = 0; i < numberOfNodes; i++) {
        int commId = i / Constant::COMMUNITY_POPULATION;
        listCommListNodeIds[commId].push_back(listNodeIds[index[i]]);
        mapNodeId2CommId[listNodeIds[index[i]]] = commId;
    }

    /*int numberOfNodes = listNodeIds.size();
    int numberOfCommunities = ceil(((double)numberOfNodes) / Constant::COMMUNITY_POPULATION);

    srand(time(NULL));
    for (int i = 0; i < numberOfNodes; i++) {
        int nodeId = listNodeIds[i];
        int comm = rand() % numberOfCommunities;
        listCommListNodeIds[comm].push_back(nodeId);
    }*/

    // setup hMax
    hMax = Constant::IS_BOUNDED_THRESHOLD ? 2 : 0;
    if (!Constant::IS_BOUNDED_THRESHOLD) {
        for (int i = 0; i < listCommListNodeIds.size(); i++) {
            if (hMax < listCommListNodeIds[i].size() * Constant::PERCENTAGE_THRESHOLD)
                hMax = (int) listCommListNodeIds[i].size() * Constant::PERCENTAGE_THRESHOLD;
            if (bMin > listCommListNodeIds[i].size())
                bMin = listCommListNodeIds[i].size();
        }
    }

    std::cout << "done reading file";
    delete[] buffer;
}

void SocialGraph::generateFile(string inputFile) {
    ifstream file;
    file.open(inputFile);
    if (file) {
        map<int, vector<int>> mapIncommingNeighbors;
        vector<int> listNodeIds;
        int startId, endId;
        int count = 0;
        while (file >> startId >> endId) {
            /*cout << count << endl;
            count++;*/
            mapIncommingNeighbors[endId].push_back(startId);
            /*
            if (mapIncommingNeighbors.find(endId) != mapIncommingNeighbors.end()) {
                mapIncommingNeighbors[endId].push_back(startId);
            }
            else {
                vector<int> tmp;
                tmp.push_back(startId);
                mapIncommingNeighbors.insert(std::pair<int, vector<int>>(endId, tmp));
            }
            */
            if (find(listNodeIds.begin(), listNodeIds.end(), startId) == listNodeIds.end())
                listNodeIds.push_back(startId);

            if (find(listNodeIds.begin(), listNodeIds.end(), endId) == listNodeIds.end())
                listNodeIds.push_back(endId);
        }

        // assign community
        int numberOfNodes = listNodeIds.size();
        int numberOfCommunities = ceil(((double) numberOfNodes) / Constant::COMMUNITY_POPULATION);

        vector<int> index;
        for (int i = 0; i < numberOfNodes; i++)
            index.push_back(i);


        std::random_shuffle(index.begin(), index.end());
        map<int, vector<int>> mapCommunities; // map community Id -> list of nodes in community
        for (int i = 0; i < numberOfNodes; i++) {
            int commId = i / Constant::COMMUNITY_POPULATION;
            mapCommunities[commId].push_back(listNodeIds[index[i]]);
        }


        //srand(time(NULL));
        //for (int i = 0; i < numberOfNodes; i++) {
        //	int nodeId = listNodeIds[i];
        //	int comm = rand() % numberOfCommunities;
        //	if (mapCommunities.find(comm) != mapCommunities.end()) {
        //		mapCommunities[comm].push_back(nodeId);
        //	}
        //	else {
        //		vector<int> tmp;
        //		tmp.push_back(nodeId);
        //		mapCommunities.insert(std::pair<int, vector<int>>(comm, tmp));
        //	}
        //}

        // write to file
        string outFileName = "data_" + to_string(numberOfNodes) + ".txt";
        ofstream writeFile(outFileName);
        if (writeFile.is_open()) {
            // write list of communities first. Format: node id - comm id
            for (map<int, vector<int>>::iterator it = mapCommunities.begin(); it != mapCommunities.end(); ++it) {
                for (int i = 0; i < it->second.size(); i++) {
                    writeFile << it->second[i] << " " << it->first << "\n";
                }
            }
            // separate list of comm and list of edges by -1
            writeFile << "-1 -1\n";

            // write list of edges. Format: startid endid weight
            for (map<int, vector<int>>::iterator it = mapIncommingNeighbors.begin();
                 it != mapIncommingNeighbors.end(); ++it) {
                for (int i = 0; i < it->second.size(); i++) {
                    writeFile << it->first << " " << it->second[i] << " ";
                    writeFile << ((double) 1) / it->second.size();
                    writeFile << "\n";
                }
            }
            writeFile.close();
        }
        file.close();
    }
}

void SocialGraph::generateFileIM(string outputFile) {
    map<int, int> mapNodeIdx;
    for (int i = 0; i < listNodeIds.size(); i++) {
        mapNodeIdx[listNodeIds[i]] = i + 1;
    }
    ofstream writeFile(outputFile);
    if (writeFile.is_open()) {
        // first line: no of nodes - no of edges
        writeFile << listNodeIds.size() << " " << noOfEdges << endl;

        // next - list of edge: start node - end node - weight
        for (int i = 0; i < listNodeIds.size(); i++) {
            int nodeId = listNodeIds[i];
            vector<std::pair<int, double>> neighbors = mapIncommingNeighbors[nodeId];
            if (neighbors.size() > 0) {
                double w = Constant::MODEL ? 1.0 / (neighbors.size() + 1) : 1. /
                                                                            neighbors.size(); // weight is different between LT and IC
                for (int j = 0; j < neighbors.size(); j++) {
                    int tmp = neighbors[j].first;
                    writeFile << mapNodeIdx[tmp] << " " << mapNodeIdx[nodeId] << " " << w << endl;
                }
            }
        }
        writeFile.close();
    }
}

void SocialGraph::standardize(string file, bool header) {
    map<int, int> mapNodeId2Index;
    for (int i = 0; i < listNodeIds.size(); i++) {
        mapNodeId2Index[listNodeIds[i]] = i;
    }
    ofstream writeFile(file);
    if (writeFile.is_open()) {
        if (header)
            writeFile << listNodeIds.size() << " " << noOfEdges << endl;

        for (map<int, vector<std::pair<int, double>>>::iterator it = mapIncommingNeighbors.begin();
             it != mapIncommingNeighbors.end(); ++it) {
            int nodeId = it->first;
            vector<std::pair<int, double>> listNei = it->second;
            for (int i = 0; i < listNei.size(); i++) {
                if ((!isDirected && nodeId < listNei[i].first) || isDirected)
                    writeFile << mapNodeId2Index[nodeId] << " " << mapNodeId2Index[listNei[i].first] << endl;
            }
        }

        writeFile.close();
    }
}

void SocialGraph::formCommunityModularity(string output, string tmpfile, bool directed) {
    string file = tmpfile;
    standardize(file);
    string tmp2 = "./ldf -i " + file + " -o " + output;
    if (directed)
        tmp2 = tmp2 + " -d";
    const char *runMMcmd = tmp2.c_str();
    cout << runMMcmd << endl;
    // system(runMMcmd);
}

void SocialGraph::formCommunityClauset(string output) {
    string file = "tmp2.adj";
    standardize(file, false);
    string tmp2 = "../snap-master/examples/community/community -i:" + file + " -o:" + output;
    const char *runCcmd = tmp2.c_str();
    system(runCcmd);
}

void SocialGraph::readCommunityFile(string file, bool isMM) {
    actualCommNodeIds.clear();
    ifstream inputFile;
    inputFile.open(file);
    if (inputFile) {
        if (isMM) { // community file generated by ldf
            int commId = 0;
            string tmp;
            vector<int> commNodes;
            while (getline(inputFile, tmp)) {
                istringstream iss(tmp);
                while (iss) {
                    int nodeIndex;
                    iss >> nodeIndex;
                    commNodes.push_back(listNodeIds[nodeIndex]);
                }
                commNodes.pop_back();
                actualCommNodeIds.push_back(commNodes);
                commNodes.clear();
            }
        } else { // community file generated by Girvan algorithms
            string tmp;

            // first 6 lines are useless
            getline(inputFile, tmp);
            getline(inputFile, tmp);
            getline(inputFile, tmp);
            getline(inputFile, tmp);
            getline(inputFile, tmp);
            getline(inputFile, tmp);

            int nodeIndex, actualCommId, prevCommmId = 0, commId = 0;

            vector<int> commNodes;
            while (getline(inputFile, tmp)) {
                istringstream iss(tmp);
                iss >> nodeIndex >> actualCommId;
                if (actualCommId == prevCommmId) {
                    commNodes.push_back(listNodeIds[nodeIndex]);
                } else {
                    actualCommNodeIds.push_back(commNodes);
                    commNodes.clear();
                    commNodes.push_back(listNodeIds[nodeIndex]);
                    prevCommmId = actualCommId;
                }

            }
            if (!commNodes.empty()) {
                actualCommNodeIds.push_back(commNodes);
            }
        }
    }
}

void SocialGraph::formCommunitiesFromActualCommunities() {
    listCommListNodeIds.clear();
    mapNodeId2CommId.clear();
    hMax = 0;
    bMin = 10000;
    vector<int> commNodes;
    int commId = 0;
    int numberOfComms = Constant::NUMBER_OF_COMMS - 1;
    for (int i = 0; i < actualCommNodeIds.size(); i++) {
        if (actualCommNodeIds[i].size() <= Constant::COMMUNITY_POPULATION) {
            addCommunity(&actualCommNodeIds[i]);
            for (int j = 0; j < actualCommNodeIds[i].size(); j++)
                mapNodeId2CommId[actualCommNodeIds[i][j]] = commId;
            commId++;
            if (commId > numberOfComms) break;
        } else {
            commNodes.clear();
            for (int j = 0; j < actualCommNodeIds[i].size(); j++) {
                int nodeId = actualCommNodeIds[i][j];
                commNodes.push_back(nodeId);
                mapNodeId2CommId[nodeId] = commId;
                if (commNodes.size() >= Constant::COMMUNITY_POPULATION) {
                    addCommunity(&commNodes);
                    commNodes.clear();
                    commId++;
                    if (commId > numberOfComms) break;
                }
            }
            if (commId > numberOfComms) break;
            if (!commNodes.empty()) {
                addCommunity(&commNodes);
                commNodes.clear();
                commId++;
                if (commId > numberOfComms) break;
            }
        }
    }
    cout << "pause" << endl;
}

int SocialGraph::randomSelectCommunity() {
    if (!Constant::IS_WEIGHTED)
        return commonInstance->randomInThread() % listCommListNodeIds.size();
    else {
        int tmp = commonInstance->randomInThread() % listNodeIds.size();
        return mapNodeId2CommId[listNodeIds[tmp]];
    }
}

vector<int> *SocialGraph::getNodesOfCommunity(int commId) {
    return &listCommListNodeIds[commId];
}

vector<std::pair<int, double>> *SocialGraph::getIncommingNeighbors(int nodeId) {
    if (mapIncommingNeighbors.find(nodeId) != mapIncommingNeighbors.end())
        return &mapIncommingNeighbors[nodeId];
    else
        return nullptr;
}

map<int, vector<std::pair<int, double>>> *SocialGraph::getMapIncommingNeighbors() {
    return &mapIncommingNeighbors;
}

vector<int> *SocialGraph::getListNodeIds() {
    return &listNodeIds;
}

int SocialGraph::getMaxThreshold() {
    return Constant::IS_BOUNDED_THRESHOLD ? 2 : hMax;
}

int SocialGraph::getNumberOfNodes() {
    return listNodeIds.size();
}

int SocialGraph::getNumberOfCommunities() {
    return listCommListNodeIds.size();
}

int SocialGraph::getCommunityThreshold(int commId) {
    if (Constant::IS_BOUNDED_THRESHOLD)
        return 2;
    else
        return listCommListNodeIds[commId].size() * Constant::PERCENTAGE_THRESHOLD;
}

int SocialGraph::getCommunityId(int nodeId) {
    return mapNodeId2CommId[nodeId];
}

int SocialGraph::getCommunitySize(int commId) {
    return listCommListNodeIds[commId].size();
}

int SocialGraph::getMinBenefit() {
    return bMin;
}

void SocialGraph::clear() {
    listNodeIds.clear();
    listCommListNodeIds.clear();
    mapIncommingNeighbors.clear();
    mapNodeId2CommId.clear();
    noOfEdges = 0;
    hMax = 0;
    bMin = 10000;

    outgoingDegree.clear();
    incomingDegree.clear();
    mapNodeWeight.clear();
    mapNodeCost.clear();
    mapNodeBenefit.clear();
    mapIncomingNodes.clear();
    mapOutgoingNodes.clear();
}

void SocialGraph::addCommunity(vector<int> *commNodes) {
    listCommListNodeIds.push_back(vector<int>(*commNodes));
    if (hMax < commNodes->size() * Constant::PERCENTAGE_THRESHOLD)
        hMax = (int) commNodes->size() * Constant::PERCENTAGE_THRESHOLD;
    if (bMin > commNodes->size())
        bMin = commNodes->size();
    double commBenefit = 0.;
    for (int i = 0; i < commNodes->size(); ++i) {
        commBenefit += mapNodeBenefit[commNodes->at(i)];
    }
    mapCommBenefit[listCommListNodeIds.size() - 1] = commBenefit;
}

void SocialGraph::generateEdgeWeightBinFile(string file, string outfile, string binfile) {
    clear();
    ifstream in(file);
    if (in) {
        int src, dst;
        int numEdges = 0;
        // int numNodes = 0;
        vector<int> srcNodes;
        vector<int> dstNodes;
        map<int, bool> nodeIds;
        map<int, vector<int>> mapIncomingNodes;
        map<int, vector<int>> mapOutgoingNodes;
        map<int, int> mapInDegree;
        map<int, int> mapOutDegree;
        map<int, double> mapDstWeight;

        // in >> numNodes >> numEdges;

        while (in >> src >> dst) {
            numEdges++;
            srcNodes.push_back(src);
            dstNodes.push_back(dst);

            if (nodeIds.find(src) == nodeIds.end()) {
                nodeIds[src] = true;
            }
            if (nodeIds.find(dst) == nodeIds.end()) {
                nodeIds[dst] = true;
            }

            if (mapIncomingNodes.find(dst) != mapIncomingNodes.end()) {
                mapIncomingNodes[dst].push_back(src);
                mapInDegree[dst]++;
            } else {
                mapIncomingNodes.emplace(dst, vector<int>{src});
                mapInDegree.emplace(dst, 1);
            }
            if (mapOutgoingNodes.find(dst) != mapOutgoingNodes.end()) {
                mapOutgoingNodes[src].push_back(dst);
                mapOutDegree[src]++;
            } else {
                mapOutgoingNodes.emplace(src, vector<int>{dst});
                mapOutDegree.emplace(src, 1);
            }
        }
        in.close();

        for (auto &incomingNodes: mapIncomingNodes) {
            double w = (double) 1 / incomingNodes.second.size();
            mapDstWeight[incomingNodes.first] = w;
        }

        int numNodes = nodeIds.size();
        ofstream out(outfile);
        out << nodeIds.size() << " " << numEdges << endl;
        for (int i = 0; i < numEdges; i++) {
            double w = mapDstWeight[dstNodes[i]];
            out << srcNodes[i] << " " << dstNodes[i] << " " << w << endl;
        }
        out.close();

        ofstream outBin(binfile, fstream::out | fstream::binary);
        // outputs number of nodes
        outBin.write((char *) &numNodes, sizeof(int));
        // outputs number of edges
        outBin.write((char *) &numEdges, sizeof(int));
        // outputs node id list
        for (auto nodeId : nodeIds) {
            int id = nodeId.first;
            outBin.write((char *) &id, sizeof(int));
        }
        // outputs weights of incoming for dstNode
        for (auto nodeId : nodeIds) {
            double weight = 0;
            if (mapInDegree.find(nodeId.first) != mapInDegree.end()) {
                weight = mapDstWeight[nodeId.first];
            }
            outBin.write((char *) &weight, sizeof(double));
        }
        // outputs node's cost
        map<int, double> mapNodeCost;
        for (auto nodeId : nodeIds) {
            double cost = 1;
            if (mapOutDegree.find(nodeId.first) != mapOutDegree.end()) {
                cost = ((double) numNodes * mapOutDegree[nodeId.first]) / numEdges;
            }
            mapNodeCost[nodeId.first] = cost;
            outBin.write((char *) &cost, sizeof(double));
        }
        // outputs node's benefit
        for (auto nodeId : nodeIds) {
            double benefit = ((double) (commonInstance->randomInThread() % 1000)) / 1000;
            // double benefit = sfmt_genrand_real1(&sfmtSeed) + mapNodeCost[nodeId.first];
            outBin.write((char *) &benefit, sizeof(double));
        }
        // outputs cumulative indegree sequence
        for (auto nodeId : nodeIds) {
            int degree = 0;
            if (mapInDegree.find(nodeId.first) != mapInDegree.end()) {
                degree = mapInDegree[nodeId.first];
            }
            outBin.write((char *) &degree, sizeof(int));
        }
        // outputs incoming of dstNode
        for (auto nodeId : nodeIds) {
            int degree = 0;
            vector<int> inNodes(degree);
            if (mapInDegree.find(nodeId.first) != mapInDegree.end()) {
                degree = mapInDegree[nodeId.first];
                inNodes = mapIncomingNodes[nodeId.first];
            }
            if (degree > 0) {
                outBin.write((char *) &inNodes[0], degree * sizeof(int));
            }
        }
        // outputs cumulative outdegree sequence
        for (auto nodeId : nodeIds) {
            int degree = 0;
            if (mapOutDegree.find(nodeId.first) != mapOutDegree.end()) {
                degree = mapOutDegree[nodeId.first];
            }
            outBin.write((char *) &degree, sizeof(int));
        }
        // outputs outgoing of dstNode
        for (auto nodeId : nodeIds) {
            int degree = 0;
            vector<int> outNodes(degree);
            if (mapOutDegree.find(nodeId.first) != mapOutDegree.end()) {
                degree = mapOutDegree[nodeId.first];
                outNodes = mapOutgoingNodes[nodeId.first];
            }
            if (degree > 0) {
                outBin.write((char *) &outNodes[0], degree * sizeof(int));
            }
        }
        outBin.close();
        cout << "Generate \"" << outfile << "\" successfully" << endl;
    } else {
        cerr << "Can't open file: " << file << endl;
    }
}

void SocialGraph::readSocialGraphBin(string file, bool isDirected) {
    this->isDirected = isDirected;
    clear();
    srand(time(NULL));
    ifstream inputFile;
    inputFile.open(file, fstream::in | fstream::binary);
    if (inputFile) {
        int noOfNodes;
        inputFile.read((char *) &noOfNodes, sizeof(int));
        if (inputFile.rdstate() != ios::goodbit) {
            cerr << "The file " << file << " is not a valid graph" << endl;
            exit(EXIT_FAILURE);
        }

        // reads number of edges
        inputFile.read((char *) &noOfEdges, sizeof(int));

        // read unique node id list
        // vector<int> listNodeIds(noOfNodes);
        listNodeIds.resize(noOfNodes);
        // mapNodeIdx.clear();
        inputFile.read((char *) &listNodeIds[0], noOfNodes * sizeof(int));
        // for (int i = 0; i < noOfNodes; ++i) {
        //     mapNodeIdx.insert(make_pair(listNodeIds[i], i));
        // }

        // reads weights of incoming for dstNode
        for (int i = 0; i < noOfNodes; ++i) {
            inputFile.read((char *) &mapNodeWeight[listNodeIds[i]], sizeof(double));
        }

        // reads costs
        vector<double> nodeCosts(noOfNodes);
        inputFile.read((char *) &nodeCosts[0], noOfNodes * sizeof(double));
        if (!Constant::GCS) {
            nodeCosts.clear();
            nodeCosts.resize(noOfNodes, 1);
        }
        for (int i = 0; i < noOfNodes; ++i) {
            mapNodeCost[listNodeIds[i]] = nodeCosts[i];
        }
        // for (int i = 0; i < noOfNodes; ++i) {
        //     inputFile.read((char *) &mapNodeCost[listNodeIds[i]], sizeof(double));
        // }

        // reads benefits
        vector<double> nodeBenefits(noOfNodes);
        inputFile.read((char *) &nodeBenefits[0], noOfNodes * sizeof(double));
        if (!Constant::GCS) {
            nodeBenefits.clear();
            nodeBenefits.resize(noOfNodes, 1);
        }
        for (int i = 0; i < noOfNodes; ++i) {
            mapNodeBenefit[listNodeIds[i]] = nodeBenefits[i];
        }
        // for (int i = 0; i < noOfNodes; ++i) {
        //     inputFile.read((char *) &mapNodeBenefit[listNodeIds[i]], sizeof(double));
        // }

        // reads cumulative indegree sequence
        incomingDegree.resize(noOfNodes);
        inputFile.read((char *) &incomingDegree[0], noOfNodes * sizeof(int));

        // reads incoming of dstNode
        for (int i = 0; i < noOfNodes; ++i) {
            mapIncomingNodes[listNodeIds[i]].resize(incomingDegree[i]);
            if (incomingDegree[i] > 0) {
                inputFile.read((char *) &mapIncomingNodes[listNodeIds[i]][0], incomingDegree[i] * sizeof(int));
            }
        }

        // reads cumulative outdegree sequence
        outgoingDegree.resize(noOfNodes);
        inputFile.read((char *) &outgoingDegree[0], noOfNodes * sizeof(int));

        // reads outgoing of dstNode
        for (int i = 0; i < noOfNodes; ++i) {
            mapOutgoingNodes[listNodeIds[i]].resize(outgoingDegree[i]);
            if (outgoingDegree[i] > 0) {
                inputFile.read((char *) &mapOutgoingNodes[listNodeIds[i]][0], outgoingDegree[i] * sizeof(int));
            }
        }

        for (map<int, vector<int>>::iterator it = mapIncomingNodes.begin(); it != mapIncomingNodes.end(); ++it) {
            int dstNode = it->first;
            vector<int> srcNodes = it->second;
            double weight = mapNodeWeight[dstNode];
            for (int i = 0; i < srcNodes.size(); ++i) {
                int srcNode = srcNodes[i];
                if (mapIncommingNeighbors.find(dstNode) != mapIncommingNeighbors.end()) {
                    mapIncommingNeighbors[dstNode].push_back(make_pair(srcNode, weight));
                } else {
                    vector<pair<int, double>> tmp;
                    tmp.push_back(make_pair(srcNode, weight));
                    mapIncommingNeighbors.insert(make_pair(dstNode, tmp));
                }
            }
        }
        inputFile.close();
    }
}
