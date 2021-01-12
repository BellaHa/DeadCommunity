/*
 * clustering.h
 *
 *  Created on: May 9, 2011
 *      Author: thang
 */

#ifndef CLUSTERING_H_
#define CLUSTERING_H_

#include <vector>
#include <algorithm>
#include <string>


struct  Triple {
            int i, j, k;
            double dif;
            Triple(int x, int y, int z, double di): i(x), j(y), k(z), dif(di) {}
            bool operator<(const Triple &b ) const {
                return this->dif < b.dif;
            }
};

typedef	 std::vector<int> Membership;	 // CS membership array
typedef	 std::vector< std::vector<int> > DisjoinCS;	// Disjoin community structure



/**
 * @brief Write pajek cluster given by a Membership array to a file
 */
bool writePajekCluster(std::string file_name, const Membership mbs);

/**
 * @brief Read a cluster from a file to a membership array
 * @param file_name
 * @return Membership array
 */
Membership loadPajekCluster(std::string file_name);

/**
 * Convert a disjoint community structure into a membership array
 * @param cs	A disjoint community structure
 * @return		The membership array
 */
Membership	cs_convert(const DisjoinCS& cs);

/**
 * Convert a membership array into a disjoint community structure
 * @param cs	A membership array
 * @return		The disjoint community structure
 */
DisjoinCS	cs_convert(const Membership& mbs);

#endif /* CLUSTERING_H_ */
