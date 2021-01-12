/*
 * clustering.cpp
 *
 *  Created on: May 9, 2011
 *      Author: thang
 */

#include "clustering.h"
#include "common.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <string>

using namespace std;


/**
 * @brief Write pajek cluster given by a Membership array to a file
 */
bool writePajekCluster(std::string file_name, const Membership mbs) {
	FILE* f= fopen(file_name.c_str(),"w");
	if (!f) {
		cerr <<"Error: Cannot create file " << file_name  << endl;
		return false;
	}
	fprintf(f, "*Vertices %d\n",mbs.size());
	FOR(i, mbs.size())
		 fprintf(f,"%d\n", mbs[i]);
	fclose(f);
	return true;
}

/**
 * @brief Read a cluster from a file to a membership array
 * @param file_name
 * @return Membership array
 */
Membership loadPajekCluster(std::string file_name) {
	Membership mbs;
	FILE* f= fopen(file_name.c_str(),"r");
	if (!f) {
		cerr <<"Error: Cannot access file " << file_name  << endl;
		return mbs;
	}
	char	buf[1000];
	int n = 0;
	fscanf(f, "%s%d\n",buf, &n);
	mbs.resize(n);
	int u;
	FOR(i, n) {
         u = -1;
		 fscanf(f,"%d\n", &u);
		 if (u == -1 && feof(f)) {
			 cerr <<"Error: Invalid format Pajek cluster file " << file_name  << endl;
			 mbs.clear();
			 return mbs;
		 }
		 mbs[i] = u;
	}
	fclose(f);
	return mbs;
}


/**
 * Convert a disjoint community structure into a membership array
 * @param cs	A disjoint community structure
 * @return		The membership array
 */
Membership	cs_convert(const DisjoinCS& cs) {
	Membership mbs;	
	int n = 0;
	for(int i =0; i < cs.size(); ++i)
		n += cs[i].size();
	mbs.resize( n, -1);
	for(int i = 0; i < cs.size(); ++i)
		 for(int j=0; j < cs[i].size(); ++j)  {
			 DCHECK( (0 <=cs[i][j]) && (cs[i][j]<n), "Error: Invalid DisjoinCS "<<cs[i][j]);
			 mbs[ cs[i][j] ] = i;
		 }
#ifdef	_DEBUG	
		for(int i = 0; i <n;++i)
				DCHECK(mbs[i] >= 0, "Error: Invalid DisjoinCS, not found "<<i);
#endif	
	return mbs;
}

/**
 * Convert a membership array into a disjoint community structure
 * @param cs	A membership array
 * @return		The disjoint community structure
 */
DisjoinCS	cs_convert(const Membership& mbs) {	
	int n = mbs.size(),
		ncs = *max_element( mbs.begin(), mbs.end() );
	DisjoinCS cs(ncs+1);
	for(int i =0; i < mbs.size();++i) {
#ifdef _DEBUG
		if (mbs[i] <0) {
			cerr <<"Error: Invalid Membership "<<i<<endl;
			exit(1);
		}
#endif
		cs[ mbs[i] ].push_back(i);
	}
#ifdef _DEBUG
	for(int i =0; i < cs.size();++i)
		if (cs[i].size() <= 0) {
			cerr << "Error: CS with empty cluster !"<<i<<endl;
			exit(1);
		}
#endif
	return cs;
}
