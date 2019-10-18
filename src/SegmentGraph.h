/*
Part of SQUID transcriptomic structural variation detector
(c) 2017 by  Cong Ma, Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __SEGMENTGRAPH_H__
#define __SEGMENTGRAPH_H__

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <limits>
#include <ctime>
#include <cmath>
#include <cassert>
#include <map>
#include <set>
#include "boost/algorithm/string.hpp"
#include "boost/icl/interval_set.hpp"
#include "boost/icl/interval_map.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/one_bit_color_map.hpp"
#include "boost/graph/stoer_wagner_min_cut.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/typeof/typeof.hpp"
#include "glpk.h"
#include "BPEdge.h"
#include "BPNode.h"
#include "ReadRec.h"
#include "SingleBamRec.h"
#include "Config.h"

using namespace std;
using namespace BamTools;

extern bool UsingSTAR;
extern int Concord_Dist_Pos;
extern int Concord_Dist_Idx;
extern int Min_Edge_Weight;
extern uint16_t ReadLen;
extern double DiscordantRatio;
extern int MaxAllowedDegree;

const int TMLIM = 3600000; //1 hour

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, boost::property<boost::edge_weight_t, int> > undirected_graph;
typedef boost::property_map<undirected_graph, boost::edge_weight_t>::type weight_map_type;
typedef boost::property_traits<weight_map_type>::value_type weight_type;

static std::map<char,char> Nucleotide={{'A','T'},{'C','G'},{'G','C'},{'T','A'},{'R','Y'},{'Y','R'},{'S','W'},{'W','S'},{'K','M'},{'M','K'},{'B','V'},{'V','B'},{'D','H'},{'H','D'}, {'N','N'}, {'.','.'},{'-','-'}};

void ReverseComplement(string::iterator itbegin, string::iterator itend);
pair<int,int> ExtremeValue(vector<int>::iterator itbegin, vector<int>::iterator itend);
void CountTop(vector< pair<int,int> >& x);

struct MinCutEdge_t
{
	unsigned long first;
	unsigned long second;
};

class SegmentGraph_t{
public:
	vector<Node_t> vNodes;
	vector<Edge_t> vEdges;
	vector<int> Label;
public:
	SegmentGraph_t(){};
	SegmentGraph_t(const vector<int>& RefLength, SBamrecord_t& Chimrecord, string ConcordBamfile);
	SegmentGraph_t(string graphfile);

	bool IsDiscordant(int edgeidx);
	bool IsDiscordant(Edge_t* edge);
	bool IsDiscordant(Edge_t edge);

	void BuildNode_STAR(const vector<int>& RefLength, SBamrecord_t& Chimrecord, string ConcordBamfile);
	void BuildNode_BWA(const vector<int>& RefLength, string bamfile);
	void BuildEdges(SBamrecord_t& Chimrecord, string ConcordBamfile);
	void FilterbyWeight();
	void FilterbyInterleaving(vector<bool>& KeepEdge);

	vector<int> LocateRead(int initialguess, ReadRec_t& ReadRec);
	vector<int> LocateRead(vector<int>& singleRead_Node, ReadRec_t& ReadRec);

	void RawEdgesChim(SBamrecord_t& Chimrecord);
	void RawEdgesOther(SBamrecord_t& Chimrecord, string ConcordBamfile);
	void RawEdges(SBamrecord_t& Chimrecord, string bamfile);
	
	void FilterEdges(const vector<bool>& KeepEdge);
	void UpdateNodeLink();
	void CompressNode();
	void CompressNode(vector< vector<int> >& Read_Node);
	void FurtherCompressNode();
	void OutputDegree(string outputfile);

	int DFS(int node, int curlabelid, vector<int>& Label);
	void ConnectedComponent(int & maxcomponentsize);
	void ConnectedComponent();
	void MultiplyDisEdges();
	void DeMultiplyDisEdges();

	void ExactBreakpoint(SBamrecord_t& Chimrecord, map<Edge_t, vector< pair<int,int> > >& ExactBP);
	void ExactBPConcordantSupport(string Input_BAM, SBamrecord_t& Chimrecord, const map<Edge_t, vector< pair<int,int> > >& ExactBP, map<Edge_t, vector< pair<int,int> > >& ExactBP_concord_support);
	void OutputGraph(string outputfile);

	vector< vector< vector<int> > > Ordering(int Mode,string Output_Prefix);
	vector< vector<int> > MincutRecursion(std::map<int,int> CompNodes, vector<Edge_t> CompEdges, int Mode, string Output_Prefix);
	void GenerateILP(std::map<int,int>& CompNodes, vector<Edge_t>& CompEdges, vector< vector<int> >& Z1,vector< vector<int> >& Z2, vector<int>& X1,vector<int>& X2,string Output_Prefix);
	void GenerateSqueezedILP(std::map<int,int>& CompNodes, vector<Edge_t>& CompEdges, vector< vector<int> >& Z, vector<int>& X);
	vector<Edge_t> GenerateSquidILP(std::map<int,int>& CompNodes, vector<Edge_t>& CompEdges, vector< vector<int> >& Z, vector<int>& X);
	void IterativeSquidILP(std::map<int,int>& CompNodes, vector<Edge_t>& CompEdges, vector< vector<int> >& Z1, vector< vector<int> >& Z2, vector<int>& X1,  vector<int>& X2);

	vector< vector<int> > SortComponents(vector< vector<int> >& Components);
	vector< vector<int> > MergeSingleton(vector< vector<int> >& Components, const vector<int>& RefLength, int LenCutOff=500000);
	/*bool MergeSingleton_Insert(int singleton, vector< vector< vector<int> > >& NewComponents, int LenCutOff); //return whether inserted or not
	bool MergeSingleton_Insert(vector<int> consecutive, vector< vector< vector<int> > >& NewComponents, int LenCutOff);*/
	bool MergeSingleton_Insert(vector<int> SingletonComponent, vector< vector<int> >& NewComponents);
	bool MergeSingleton_Insert(vector< vector<int> >  Consecutive, vector< vector<int> >& NewComponents);
	vector< vector<int> > MergeComponents(vector< vector<int> >& Components, int cutoff=5);

	// small functions
	int GroupConnection(int node, vector<Edge_t*>& Edges, int sumweight, vector<int>& Connection, vector<int>& Label);
	void GroupSelect(int node, vector<Edge_t*>& Edges, int sumweight, int count, vector<int>& Connection, vector<int>& Label, vector<Edge_t>& ToDelete);

	vector<Edge_t> findEdges(int ind1, int ind2);
	int testConcordant(vector< vector< vector<int> > >& Components);
	vector<Edge_t> getDiscordant();

	// Approximation
	int is_Discordant_comp(map<int, int> order, map<int, int> orientation, Edge_t edge, int u, int pos, int ori_u);
	int is_Discordant_comp(map<int,int> order, map<int, int> orientation, Edge_t edge);
	void change_order_comp(map<int, int>& order, int u, int pos);
	vector<Edge_t> SCAP_approx(map<int,int>& CompNodes, vector<Edge_t>& CompEdges, vector<int>& ordering, int Mode);
	vector<vector<int>> DCAP_approx(map<int,int>& CompNodes, vector<Edge_t>& CompEdges, int Mode);

};

#endif
