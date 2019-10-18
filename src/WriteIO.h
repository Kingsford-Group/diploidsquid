/*
Part of SQUID transcriptomic structural variation detector
(c) 2017 by  Cong Ma, Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __WRITEIO_H__
#define __WRITEIO_H__

#include "SegmentGraph.h"

using namespace std;

extern int Concord_Dist_Pos;
extern int Concord_Dist_Idx;

vector< vector<int> > ReadComponents(string file);
void WriteComponents(string outputfile, vector< vector< vector<int> > > Components);
void WriteBEDPE(string outputfile,SegmentGraph_t& SegmentGraph, vector< vector< vector<int> > >& Components, vector< vector< pair<int, int> > >& Node_NewChr, vector<string>& RefName, map<Edge_t, vector< pair<int,int> > >& ExactBP, map<Edge_t, vector< pair<int,int> > >& ExactBP_concord_support);
void OutputNewGenome(SegmentGraph_t& SegmentGraph, vector< vector< vector<int> > >& Components, const vector<string>& RefSequence, const vector<string>& RefName, vector<string> fnames);
void TmpWriteBEDPE(string outputfile, SegmentGraph_t& SegmentGraph, vector<string>& RefName);

#endif
