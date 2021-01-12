/*
 * ldf_class.h
 *
 *  Created on: May 10, 2012
 *      Author: Thang N. Dinh
 */

#include "clustering.h"
#include "generalGraph.hpp"
#include "graphAlgorithm.hpp"
#include "edgeMap.hpp"
#include "mappedHeap.hpp"
#include "config.hpp"
#include "clustering.h"
#include "common.h"
#include "Timer.h"
#include <ctime>
#include <queue>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <fstream>
#include <functional>
#include <utility>

using namespace std;

#ifndef LDF_CLASS_H_
#define LDF_CLASS_H_


#define		SQR(x) ( (x)*(x))	// Square macro

/**
 * Class to implement LDF-based algorithms
 *
 */
template<class T>
class Ldf {
	GeneralGraph<T> *g; // Pointer to graph G
	int n; // Network size
	/**
	 * Abstract representation of G. Each vertex corresponds to a community
	 * (a subset of vertices) of G. The nodes' ids are given by the ids of the
	 * cp vector.
	 */
	GeneralGraph<T> ga; // The abstract graph
	vector<double> volc; // Volume of each community
	double inside_edge, m_out;

	/* Communities (or vertex partition) of the network
	 * cp[i]: The list of nodes in the ith community	 *
	 */
	DisjoinCS cp;
	Membership mbs; // Membership array, which community a node belongs to
public:	
	DisjoinCS best_cp;
	Membership best_mbs;
	double best_q, current_q; // Modularity values

	/**
	 * Decide how a node choose one neighbor as
	 * its leader.
	 * JP_RANDOM: The leader is chosen randomly among all neighbors
	 * JP_BEST_Q: The leader is chosen so that the modularity increases
	 * the most.
	 */
	enum JT_ALGO {
		JTA_LDF_RD,			// Ldf
		JTA_LDF_GR,			// Ldf best modularity increase
		JTA_LS,				// Local search
		JTA_LDF_RD_LS,		// Ldf_random + local search
		JTA_LDF_GR_LS,		// Ldf greedy + local search
		JTA_PRESET,			// Planned CS with preset_com() +localsearh
		JTA_NONE			// Doing nothing
	};
	enum	Vt_label  {VT_LEADER, VT_MEMBER, VT_ORBITER, VT_UNLABELED};

	/**
	 * @return The modularity value of the current clustering
	 * specified by the membership array "mbs"
	 */
	double get_modularity_check() {
		// TODO: revise for directed network
		double m = 0, mi = 0;
		// fraction of edges inside communities
		ALL_EDGES(*g, h, t)
			{
				m += g->edw(id);
				if (mbs[h] == mbs[t])
					mi += g->edw(id);
			}
		// volume of each community
		vector<double> vol(cp.size(), 0);
		MY_FOR(u, n)
			vol[mbs[u]] += g->wdeg_out(u);
		double vol2 = 0;
		MY_FOR(i, cp.size() )
			vol2 += (vol[i] * vol[i]);
		if (fabs(inside_edge - mi) > 1e-8) {
			cout <<"Error: Mismatched internal edges : "<< inside_edge
				<< mi << endl;
		}
		if ( fabs(m-m_out) > 1e-8) {
			cout <<"Error: Mismatched total edges "<< m <<" " << mi << endl;
		}
		return (mi / m - vol2 / (m * m));
	}

	/**
	 * @return The modularity value of the current clustering
	 * specified by the membership array "mbs"
	 */
	double get_modularity() {
		// TODO: revise for directed network
		double vol2 = 0;
		MY_FOR(i, volc.size() )
			vol2 += (volc[i] * volc[i]);
		return (inside_edge - vol2 / m_out) / m_out;
	}

	/**
	 * Provided the algorithm with a given CS
	 */
	double preset_com (const Membership &mb) {
		mbs = mb;
		cp = cs_convert(mbs);
		m_out = inside_edge = 0;
		ALL_EDGES(*g, h, t)
		{
			m_out += g->edw(id);
			if (mbs[h] == mbs[t])
			inside_edge += g->edw(id);
		}
		volc = vector<double> (n, 0);
		FOR(u, n)
		volc[mbs[u]] += g->wdeg_out(u);

		best_cp = cp;
		best_mbs = mbs;
		best_q = current_q = get_modularity();
		return current_q;
	}

	/**
	 * Initiate a clustering in which each vertex is in
	 * its own community. This initial solution will be
	 * assigned as the currently best community structure.
	 * @return The modularity value of the corresponding clustering
	 */
	double singleton_cluster() {
		Membership mb  = vector<int> (n, -1),
				   perm= vector<int> (n);

		// Put all isolated community into the first community
		int c = 0;
		FOR(i, n) {
			perm[i] = i;
			if (g->degree(i) <= 0) {
				mb[i] = 0;
				c = 1;
			}
		}
		random_shuffle(perm.begin(), perm.end());
		FOR(i, n)
			if (mb[ perm[i] ] < 0)
				mb[ perm[i] ] = c++;
		return 	preset_com(mb);
	}

	/**
	 * Default constructor
	 * @param ig	Pointer to a weighted graph for communities detection
	 */
	Ldf(GeneralGraph<T> *ig) :
		g(ig) {
		n = g->size();
	}

	/**
	 * Create the abstract graph GA from G
	 */
	void build_abstract_graph() {
		remove_empty_cluster();
		// Make a list of edges
		vector<int> head, tail;
		vector<T> w;
		ALL_EDGES(*g, h, t)	{
			head.push_back(mbs[h]);
			tail.push_back(mbs[t]);
			w.push_back(g->edw(id));
		}
		ga	= GeneralGraph<T> (cp.size(), head.size(), &head[0],
								&tail[0], &w[0]);
		// Aggregate the weight + allow loops
		ga.simplify(true, true);
	}

	void	reset_membership() {
		MY_FOR(c, cp.size() )
				MY_FOR(i, cp[c].size() )
					mbs[ cp[c][i] ] = c;
	}


	bool check_modularity() {
#ifdef	_DEBUG
		if ( fabs(current_q - get_modularity())>1e-8
			|| fabs(current_q - get_modularity_check()) >1e-8) {
				cout  <<"Error: Q "<< current_q <<"\t"<<get_modularity()<<"\t"
					<<get_modularity_check()<<endl;
				cout <<"CP: "<<cp.size()<<endl;				
				MY_FOR(i, n)
					cout << i<<" \t";
				cout << endl;
				MY_FOR(i, n)
					cout << mbs[i] <<" \t";
				cout << endl;
				cout.flush();				
				get_modularity();
				get_modularity_check();
				return false;				
		}
#endif
		return true;
	}

	/**
	 * Move nodes from one component to an another if doing so
	 * increases the modularity
	 */
	void local_search() {
		bool improved;				
		vector<T> wc(2*n, 0);
		vector<int> mask(2*n, -1);
		vector<int> moved(2*n, 0);
		vector<int>	nb_com;
		nb_com.reserve(n);
		do {
			std::fill(mask.begin(), mask.end(), -1);
			std::fill(moved.begin(), moved.end(), 0);
			improved = false;			
			FOR(c, cp.size() ) {
				for (int i = cp[c].size() - 1; i >= 0; --i) {
					int u = cp[c][i];
					if (moved[u])
						continue;
					moved[u] = 1;					
					// Compute the total weights of edges from u to
					// each communities					
					
					/*
					NEIGHBORS(*g, u, v) {
						if (u==v) continue; // Skip self-loops
						int idv = mbs[v];
						if (mask[idv] != u) {
							mask[idv] = u;
							wc[idv] = g->edw(id);
						} else
							wc[idv] += g->edw(id);
					} // Neighbor
					*/
					nb_com.clear();
					int *nu=(*g)[u], *ru=(*g)[u+1];
					T*	we = g->Weight(u);
					wc[c] = 0;
					mask[c] = u;
					for( ;nu < ru; ++nu, ++we)
						 if ( (u != *nu) ) {
							int idv = mbs[*nu];
							if (mask[ idv] != u ) {
								mask[ idv] = u;
								wc[idv] = 0;
								nb_com.push_back(idv);
							} 
							wc[idv] += *we;
						}					
					// Find the best community to move u to
					int nc = c;					
					// Avoid moving to the original community
					//int u2n = u + 2*n;
					//mask[c] = u2n;
					double max_dq = 0.0;					
					double du = g->wdeg_out(u);
					/*
					NEIGHBORS(*g, u, v) {
						int idv = mbs[v];
						if (mask[idv] != u + 2 * n) {
							mask[idv] = u + 2 * n;
							// considered moving u to idv component
							//double dq = 2*(wc[idv] - wc[c])/m_out
							//-  ( (SQR(volc[idv]+du) + SQR(volc[c]-du))
							//   - (SQR(volc[idv]) + SQR(volc[c])))/SQR(m_out);
							double dq =2*( (wc[idv] - wc[c])
							-   du*(volc[idv]+du-volc[c])/m_out)/m_out;

							if (dq > max_dq) {
								nc = idv;
								max_dq = dq;
							}
						}
					} // Neighbor
					*/					
					/*
					for(nu=(*g)[u], we = g->Weight(u); nu < ru; ++nu, ++we) {
						int idv = mbs[*nu];
						if (mask[idv] != u2n) {
							mask[idv] = u2n;
							// considered moving u to idv component
							//double dq = 2*(wc[idv] - wc[c])/m_out
							//-  ( (SQR(volc[idv]+du) + SQR(volc[c]-du))
							//   - (SQR(volc[idv]) + SQR(volc[c])))/SQR(m_out);
							double dq =2*( (wc[idv] - wc[c])
							-   du*(volc[idv]+du-volc[c])/m_out)/m_out;

							if (dq > max_dq) {
								nc = idv;
								max_dq = dq;
							}
						}
					} // Neighbor
					*/
					for(int j = 0, ns = nb_com.size(); j <ns;++j) {
							int idv = nb_com[j];							
							double dq =2*( (wc[idv] - wc[c])
							-   du*(volc[idv]+du-volc[c])/m_out)/m_out;
							if (dq > max_dq) {
								nc = idv;
								max_dq = dq;
							}										
					} 
					if (nc != c) {
						improved = true;
						volc[c]  -= du;
						volc[nc] += du;
						inside_edge += 2*(wc[nc] - wc[c]);
						current_q += max_dq;
						cp[c][i] = cp[c].back();
						cp[c].pop_back();
						cp[nc].push_back(u);
						mbs[u] = nc;
					}
				} // for
			} // FOR
			remove_empty_cluster();
			//decompose();
			//build_abstract_graph();
		} while (improved);
		if (current_q > best_q) {
			best_q = current_q;
			best_mbs = mbs;
			best_cp = cp;
		}

		DCHECK(check_modularity(),"Error  local_search")			
}

	void remove_empty_cluster() {
		vector<vector<int> > ncp(0);
		FOR(i, cp.size() )
			if (!cp[i].empty())
				ncp.push_back(cp[i]);
		cp = ncp;
		volc = vector< double> (cp.size(), 0);		
		FOR(i, cp.size() ) {
			FOR(j, cp[i].size() ) {
				int u = cp[i][j];
				mbs[u] = i;
				volc[i] += g->wdeg_out(u);
			}		
		}		
		DCHECK(check_modularity(), "Error: remove empty cluster")		
	}



	void ldf_best_increase() {
		// Sort vertices in non-decreasing order of the (weighted degree).
		typedef	pair<double, int> D_id;
		vector<D_id>	d_id(n);
		MY_FOR(i, n)
			d_id[i]= make_pair( (double) g->wdeg_out(i), (int) i);
		sort(d_id.begin(), d_id.end() );
		// Assign label for each vertex until no gain in modularity
		std::vector <Vt_label> lb(n, VT_UNLABELED);
		//vector<Vt_label>  lb;
		vector	<int>		parent(n, -1);
		singleton_cluster(); // Initiate the singleton communities
		vector<T> 	wc(2 * n, 0);
		vector<int> mask(2 * n, -1);
		MY_FOR(i, n) {
			int u = d_id[i].second;
			if ((lb[u] == VT_UNLABELED) && (g->degree(u)>0)) {
				int c = mbs[u];
				// Find the best neighbor to join u to
				// Compute the total weights of edges from u to
				// each communities
				wc[c] = 0;
				NEIGHBORS(*g, u, v) {
					if (u == v)
						continue; // Skip self-loops
					int idv = mbs[v];
					if (mask[idv] != u) {
						mask[idv] = u;
						wc[idv] = g->edw(id);
					} else
						wc[idv] += g->edw(id);
				} // Neighbor
				// Find the best community to move u to
				int nc = c, u_parent = -1;
				double max_dq = -m_out, du = g->wdeg_out(u);
				// Avoid moving to the original community
				mask[c] = u + 2 * n;
				NEIGHBORS(*g, u, v) {
					int idv = mbs[v];
					if ( (mask[idv] != u + 2 * n) && ( (lb[v] == VT_UNLABELED)
							|| lb[v]== VT_LEADER) ) {
						mask[idv] = u + 2 * n;
						// Considered moving u to idv component
						/*
						double dq = 2 * (wc[idv] - wc[c]) / m_out
								- ((SQR(volc[idv]+du) + SQR(volc[c]-du))
										- (SQR(volc[idv]) + SQR(volc[c])))
										/SQR(m_out);
										*/
						double dq =2*( (wc[idv] - wc[c])
							-  (volc[idv]+du-volc[c])*du/m_out)/m_out;

						if (dq > max_dq) {
							nc = idv;
							max_dq = dq;
							u_parent = v;
						}
					}
				} // Neighbor
				if (nc == c) { // Cannot select a neighbor to be the leader
					NEIGHBORS(*g, u, v) {
						int idv = mbs[v];
						if ( (mask[idv] != u + 2 * n) &&(lb[v]==VT_MEMBER) ) {
							mask[idv] = u + 2 * n;
							// Considered moving u to idv component
							double dq = 2 * (wc[idv] - wc[c]) / m_out
									- ((SQR(volc[idv]+du) + SQR(volc[c]-du))
											- (SQR(volc[idv]) + SQR(volc[c])))
											/SQR(m_out);

							if (dq > max_dq) {
								nc = idv;
								max_dq = dq;
								u_parent = v;
							}
						}
					} // Neighbor
				}
				if ( (nc != c) && (max_dq > 0) ) {
					volc[c] -= du;
					volc[nc] += du;
					inside_edge += 2 * (wc[nc] - wc[c]);
					current_q += max_dq;					
					cp[c].pop_back();
					cp[nc].push_back(u);
					if (lb[u_parent] == VT_UNLABELED)
						lb[u_parent] = VT_LEADER;
					if (lb[u_parent] == VT_LEADER)
						lb[u] = VT_MEMBER;
					else lb[u] = VT_ORBITER;
					mbs[u] = nc;
				}
			}
		}
		DCHECK(check_modularity(),"Error: ldf best increase")
		remove_empty_cluster();
		if (current_q > best_q) {
				best_q = current_q;
				best_mbs = mbs;
				best_cp = cp;
		}
	}

	void ldf_random() {
		// TODO:
	}

	/**
	 * Find clustering with the maximum modularity
	 * using multi-level ldf algorithm
	 */
	void find_community(int level=0, JT_ALGO first_phase=JTA_LDF_GR_LS,
			JT_ALGO compact_phase = JTA_LS, bool verbose=false) {		
		// Select the joining type		
		JT_ALGO jta = first_phase;
		if (level >0) 
				jta = compact_phase;
		switch (jta) {
		 case JTA_LDF_RD: ldf_random(); break;
		 case JTA_LDF_GR: ldf_best_increase(); break;
		 case JTA_LS: singleton_cluster(); local_search(); break;
		 case JTA_LDF_RD_LS: ldf_random(); local_search(); break;
		 case JTA_LDF_GR_LS: ldf_best_increase(); local_search();break;
		 case JTA_PRESET: local_search(); break;
		 case JTA_NONE:	singleton_cluster(); return; // Do nothing
		 default: singleton_cluster(); break;
		} 		
		if (verbose)
			cout << setprecision(4)<< "Level "<< level<<":  Q="
				 << current_q <<"\t Size=("<<n<<", "<<g->Arcs()<<")"<<endl;
		if (cp.size() >= n) // Stop when no improvement is possible
			return;
		// Condense the graph and go to next phase
		build_abstract_graph();
		Ldf<T> dga(&ga);
		dga.find_community(level+1, first_phase, compact_phase, verbose);
		DCHECK(check_modularity(),"Error: find_community after dga.find_com")
		// Expand the condensed graph
		vector<vector<int> > ncp;
		inside_edge = dga.inside_edge;
		current_q = dga.current_q;
		FOR(i, dga.cp.size() ) {
			vector<int> cc(0);
			FOR(j, dga.cp[i].size() ) {
				int ic = dga.cp[i][j];
				cc.insert(cc.end(), cp[ic].begin(), cp[ic].end());
			}
			ncp.push_back(cc);
		}
		cp = ncp;
		remove_empty_cluster();
		switch (jta) {
			case JTA_LS:
			case JTA_LDF_RD_LS:
			case JTA_LDF_GR_LS:	local_search(); break;
		}
	}
};

#endif /* LDF_CLASS_H_ */
