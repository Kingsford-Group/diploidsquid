/*
Part of Diploid-SQUID transcriptomic structural variation detector
(c) 2019 by Yutong Qiu, Cong Ma, Han Xie, Carl Kingsford and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "SingleBamRec.h"
#include "ReadRec.h"
#include "BPNode.h"
#include "BPEdge.h"
#include "SegmentGraph.h"
#include "WriteIO.h"
#include "Config.h"

using namespace std;

int main(int argc, char* argv[]){
	time_t CurrentTime;
	string CurrentTimeStr;
	
	bool success=parse_arguments(argc, argv);

	if(success){
		map<string, int> RefTable;
		vector<string> RefName;
		vector<int> RefLength;

		BuildRefName(Input_BAM, RefName, RefTable, RefLength);
		for(map<string,int>::iterator it=RefTable.begin(); it!=RefTable.end(); it++)
			cout<<"Reference name "<<it->first<<"\t-->\t"<<it->second<<endl;

		SBamrecord_t Chimrecord;
		if(Input_Chim_BAM.size()!=0){
			BuildChimericSBamRecord(Chimrecord, RefName, Input_Chim_BAM);
		}

		SegmentGraph_t SegmentGraph(RefLength, Chimrecord, Input_BAM);

		if(Print_Graph)
			SegmentGraph.OutputGraph(Output_Prefix+"_graph.txt");

		cout << Mode << endl;

		vector< vector< vector<int> > > Components=SegmentGraph.Ordering(Mode, Output_Prefix);
		cout << "We get ordering" << endl;

		cout << "Start sorting" << endl;

		// reorder the components so that the shape is (# alleles, # components, # nodes in each components)
		vector<vector<int>> Components1, Components2;
		vector<vector<vector<int>>> reorderComponents;
		for (int i = 0 ;i<Components.size();i++){
			Components1.push_back(Components[i][0]);
			Components2.push_back(Components[i][1]);
		}
		reorderComponents.push_back(Components1);
		reorderComponents.push_back(Components2);

		if(Print_Components_Ordering)
			WriteComponents(Output_Prefix+"_component_pri.txt", reorderComponents);
		
		for (int i = 0;i<reorderComponents.size();i++){
			reorderComponents[i]=SegmentGraph.SortComponents(reorderComponents[i]);
			reorderComponents[i]=SegmentGraph.MergeSingleton(reorderComponents[i], RefLength);
			reorderComponents[i]=SegmentGraph.SortComponents(reorderComponents[i]);
			reorderComponents[i]=SegmentGraph.MergeComponents(reorderComponents[i]);
		}

		// For each node, store its component and index in component at its index
		vector< vector< pair<int, int> > > Node_NewChr; Node_NewChr.resize(2);
		Node_NewChr[0].resize(SegmentGraph.vNodes.size());
		Node_NewChr[1].resize(SegmentGraph.vNodes.size());
		for(unsigned int i=0; i<reorderComponents.size(); i++)
			for (unsigned int k=0; k<reorderComponents[i].size();k++)
				for(unsigned int j=0; j<reorderComponents[i][k].size(); j++)
					Node_NewChr[i][abs(reorderComponents[i][k][j])-1]=make_pair(k, j); 

		if(Print_Total_Ordering)
			WriteComponents(Output_Prefix+"_component.txt", reorderComponents);

		map<Edge_t, vector< pair<int,int> > > ExactBP;
		SegmentGraph.ExactBreakpoint(Chimrecord, ExactBP);
		map<Edge_t, vector< pair<int,int> > > ExactBP_concord_support;
		SegmentGraph.ExactBPConcordantSupport(Input_BAM, Chimrecord, ExactBP, ExactBP_concord_support);
		SegmentGraph.DeMultiplyDisEdges();
		WriteBEDPE(Output_Prefix+"_sv.txt", SegmentGraph, reorderComponents, Node_NewChr, RefName, ExactBP, ExactBP_concord_support);

		if(Print_Rearranged_Genome){
			vector<string> RefSequence;
			bool canbuild=BuildRefSeq(Input_FASTA, RefTable, RefLength, RefSequence);
			if(canbuild){
				vector<string> fnames{ Output_Prefix+"_1_genome.fa",  Output_Prefix+"_2_genome.fa"};
				OutputNewGenome(SegmentGraph, reorderComponents, RefSequence, RefName, fnames);
			}
		}

		time(&CurrentTime);
		CurrentTimeStr=ctime(&CurrentTime);
		cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] Done."<<endl;

	}
}
