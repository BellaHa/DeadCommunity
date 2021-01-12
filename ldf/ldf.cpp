/* Copyright June 2012, Thang N. Dinh */
/* Last modified: Aug. 2016 */

/*
 * INPUT: Graph G=(V,E).
 * OUTPUT: Community structure by maximizing modularity with the LDF algorithm
 * in the paper
 * 
 *    T. N. Dinh and M. T. Thai, Community Detection in Scale-free Networks: Approximation Algorithms for
 *    Maximizing Modularity, IEEE J. on Selected Areas in Com.: Spec. Iss. on Network Science, 2013
 */
#pragma warning(1: 4710)

#include "clustering.h"
#include "generalGraph.hpp"
#include "graphAlgorithm.hpp"
#include "edgeMap.hpp"
#include "config.hpp"
#include "Timer.h"
#include "modularity.hpp"

#include <ctime>
#include <queue>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>
#include <fstream>

using namespace std;


template<class T>
void	statistic(vector<T> v, int prc=4) {
	if (v.empty() )
		cout <<" N/A "<<endl;
	else if (v.size() == 1)
		cout <<v[0]<<" n=1"<<endl;
	else {
		double mean = accumulate(ALL(v), 0.0)/(1.0*v.size());
		T	min_e = *min_element(ALL(v) ),
			max_e = *max_element(ALL(v) );
		double std = 0.0;
		FOR(i, v.size() )
			std += (double)(v[i]-mean)*(v[i]-mean);	
		cout <<setprecision(prc)<<" avg="<<mean<<" min="<<min_e <<" max="<<max_e
			 <<" std="<<sqrt(std/(v.size()-1.0))<<" n="<<v.size()<<endl;
	}
}



Configuration    *par;

typedef		float	EdgeType;

template <class T>
typename Ldf<T>::JT_ALGO  algo_type(string key) {
  if (par->match(key,"ldf"))
	  return Ldf<T>::JTA_LDF_RD;
  if (par->match(key,"ldf_gr"))
	  return Ldf<T>::JTA_LDF_GR;
  if (par->match(key,"local"))
	  return Ldf<T>::JTA_LS;  
  if (par->match(key,"ldf_lc"))
	  return Ldf<T>::JTA_LDF_RD_LS;
  if (par->match(key,"ldf_gr_lc"))
	  return Ldf<T>::JTA_LDF_GR_LS;
  if (par->match(key,"preset"))
	  return Ldf<T>::JTA_PRESET;
  if (par->match(key,"none"))
	  return Ldf<T>::JTA_NONE;
  return Ldf<T>::JTA_LS;
}

template	<class T>
T	max_vector(const vector<T> &v) {
	return *max_element( ALL(v) );
}


void write_cs_to_file(Membership cl, string fileName, vector<bool> isolated) {
	ofstream	fo(fileName.c_str() );
	CHECK(fo.good(), "Error: Cannot create file "<< fileName)
	DisjoinCS best_cs = cs_convert(cl);
	FOR(i, best_cs.size() )   {
		int count = 0;
		FOR(j, best_cs[i].size())
			if (!isolated[ best_cs[i][j] ] )
				count++;
		if (count > 0) {
			FOR(j, best_cs[i].size())
				if (!isolated[ best_cs[i][j]] )
					fo << best_cs[i][j]<<" ";	
			fo << endl;
		}
	}
	FOR(i, isolated.size() )
		if (isolated[i])
			fo << i  << endl;
	fo.close();
}

template <class T>
inline string toStr(T x){
	ostringstream ostr;
	ostr << x;
	return ostr.str();
}


int       main(int   argc,     char *argv[])
{
  Configuration    cfg(argc, argv);
  par = &cfg;
  if (cfg.exists("help")) {
	cout <<"Clustering via maximizing modularity with Low-Degree Following Algorithm."<<endl
		 <<"Usage: "<<argv[0]<<" -i <input file>  [options]"<<endl
		 <<" -i <input file>         A file contains the adjacency list of the graph. See"<<endl
		 <<"                      the file format section for more information"<<endl
     <<" -o <output file>        The output file "<<endl
		 <<" -d                      The graph is directed"<<endl		 	 
		 <<" -we                     The graph is weighted"<<endl
		 <<" -clu        		     Produce the cluster file in the Pajek format."<< endl
		 <<" -config <config file>   Load the configuration from the specified file."<<endl
		 <<" -v						 Verbose mode. Print out more details"<<endl    	      
     <<" -repeat <#time>       Repeat the algorithm multiple times"<<endl
     <<" AND IF YOU KNOW WHAT YOU'RE DOING" << endl     
     <<" -l0 <ldf| ldf_gr | local | ldf_lc | ldf_gr_lc  > Specify the algorithm for the first phase."<<endl		 
		 <<"					  The default value is 'ldf_gr_lc'."<<endl
		 <<" -lm <ldf| ldf_gr | local | ldf_lc | ldf_gr_lc | none > Specify the algorithm "<<endl
		 <<"					  to work on the compact representation. The default is 'local'. "<<endl;		 
      return 0;
  }
  // Setting the aliases
  string kn[] = {"i",     "o",    "d",    	  "t" };
  string as[] = {"input","output","directed","type"};
  for(int i=0; i < 4;++i)
      cfg.alias(kn[i], as[i]);
  // Check for required parameters
  if ( !cfg.exists("input") ) {
      cerr <<"Error: Invalid syntax! See: "
              << argv[0]<<" -help"<<endl;
      return 0;
  }
  /* Setting default parameters if necessary */
  if (!cfg.exists("output"))
      cfg.add("output", cfg.value("input")+".com");
  // Get the graph from file
  GeneralGraph<EdgeType> g( cfg.value("input"), cfg.exists("we") );
  if (!cfg.exists("we"))
        g.setWeight(1);

  if (!cfg.exists("directed") )
      g.toUndirected();
  else {
	  //TODO: directed
	  cout<<"Error: Directed networks are not currently  supported!"<<endl;
	  return 1;
  }

  // We may repeat the same algorithm multiple times and repeat the min/max/avg.
  int rp=1;
  if (cfg.exists("repeat"))
	  rp = cfg.cast<int>("repeat");

  Ldf<EdgeType>	dg(&g);
  
  vector<double>	vq, vt, vc;  
  DisjoinCS		best_cs;
  double		  best_q=-1;
  FOR(i, rp) {
    Timer tm;	
		srand(i);
    if (!cfg.exists("l0") || !cfg.exists("lm") )
      dg.find_community();
    else
      dg.find_community(0, algo_type <EdgeType>("l0"), algo_type<EdgeType>("lm"));
		vq.push_back( dg.best_q );		
		vt.push_back( tm.stop() );
		vc.push_back( dg.best_cp.size() );
		if (cfg.exists("v")) {
			cout <<"#"<<i<<"\t Q="<< dg.best_q<<" t="
			 <<  vt.back() <<" second(s) #com="<<vc.back()<<endl;
		}
		if (dg.best_q > best_q) {
			best_q = dg.best_q;
			best_cs = dg.best_cp;
		}
  }
  if (rp > 1) {	 
    cout <<"Statistic info. for multiple runs: "<<endl;
	  cout <<"Q: ";
	  statistic(vq);
	  cout <<"T: ";
	  statistic(vt);
	  cout <<"#com: ";
	  statistic(vc);
  } else {
    cout <<"Modularity Q: " << best_q <<"\t"	  
	       <<"#community: "<< best_cs.size()<<"\t"
         <<"time: "<<vt[0] <<endl;	  
   }
  ofstream	fo(cfg.value("output").c_str());
  cout <<"Writing output to file..."<<cfg.value("output")<<endl;
  CHECK(fo.good(), "Error: Cannot create file "<< cfg.value("output"))
  FOR(i, best_cs.size() ) 
  {
    int count = 0;
    FOR(j, best_cs[i].size())
	if (g.degree( best_cs[i][j]) > 0 )
		count++;
    if (count > 0) {
	    FOR(j, best_cs[i].size())
		if (g.degree( best_cs[i][j]) > 0 )
			fo << best_cs[i][j]<<" ";	
	    fo << endl;
    }
  }
  FOR(i, g.size() )
    if (g.degree(i)<= 0)
	fo << i << endl;
  fo.close();  
  return 0;
}

