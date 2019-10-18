/*
Part of Diploid-SQUID transcriptomic structural variation detector
(c) 2019 by Yutong Qiu, Cong Ma, Han Xie, Carl Kingsford and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "WriteIO.h"

using namespace std;

vector< vector<int> > ReadComponents(string file){
	ifstream input(file);
	string line;
	vector< vector<int> > Components;
	while(getline(input,line)){
		if(line[0]=='#')
			continue;
		else{
			size_t start=line.find_first_of('\t');
			line=line.substr(start+1);
			vector<string> strs;
			boost::split(strs, line, boost::is_any_of(","));
			vector<int> tmp;
			for(int i=0; i<strs.size(); i++)
				tmp.push_back(stoi(strs[i]));
			Components.push_back(tmp);
		}
	}
	input.close();
	return Components;
};

void WriteComponents(string outputfile, vector< vector< vector<int> > > Components){
	ofstream output(outputfile, ios::out);
	output<<"# component_id\tnodes_1\tnodes_2\n";
	for(int i=0; i<Components[0].size(); i++){
		output<<i<<'\t';
		for (int k=0; k<Components.size(); k++){
			for(int j=0; j<Components[k][i].size()-1; j++)
				output<<Components[k][i][j]<<",";
			output<<Components[k][i][Components[k][i].size()-1]<<"\t";
		}
		output<<endl;
	}
	output.close();
};

void WriteBEDPE(string outputfile,  SegmentGraph_t& SegmentGraph, vector< vector< vector<int> > >& Components, vector< vector< pair<int, int> > >& Node_NewChr, 
	vector<string>& RefName, map<Edge_t, vector< pair<int,int> > >& ExactBP, map<Edge_t, vector< pair<int,int> > >& ExactBP_concord_support)
{
	sort(SegmentGraph.vEdges.begin(), SegmentGraph.vEdges.end(),  [](Edge_t a, Edge_t b){return a.Weight>b.Weight;});
	ofstream output(outputfile, ios::out);
	output<<"# chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tnum_concordantfrag_bp1\tnum_concordantfrag_bp2\n";
	for(int i=0; i<SegmentGraph.vEdges.size(); i++){
		int ind1=SegmentGraph.vEdges[i].Ind1, ind2=SegmentGraph.vEdges[i].Ind2;
		bool flag_chr=(SegmentGraph.vNodes[ind1].Chr==SegmentGraph.vNodes[ind2].Chr);
		bool flag_ori=(SegmentGraph.vEdges[i].Head1==false && SegmentGraph.vEdges[i].Head2==true);
		bool flag_dist=(SegmentGraph.vNodes[ind2].Position-SegmentGraph.vNodes[ind1].Position-SegmentGraph.vNodes[ind1].Length<=Concord_Dist_Pos || ind2-ind1<=Concord_Dist_Idx);
		if(!flag_chr || !flag_ori || !flag_dist){ //if discordant
			pair<int,int> pos11=Node_NewChr[0][SegmentGraph.vEdges[i].Ind1]; //pos -> (comp, node)
			pair<int,int> pos12=Node_NewChr[0][SegmentGraph.vEdges[i].Ind2];

			pair<int,int> pos21=Node_NewChr[1][SegmentGraph.vEdges[i].Ind1];
			pair<int,int> pos22=Node_NewChr[1][SegmentGraph.vEdges[i].Ind2];

			// Record BP if edge is made concordant after either arrangement
			bool flag1=false;
			bool flag2=false;  
			if(pos11.first==pos12.first && pos11.second<pos12.second && SegmentGraph.vEdges[i].Head1==(Components[0][pos11.first][pos11.second]<0) && SegmentGraph.vEdges[i].Head2==(Components[0][pos12.first][pos12.second]>0))
				flag1=true;
			else if(pos11.first==pos12.first && pos11.second>pos12.second && SegmentGraph.vEdges[i].Head2==(Components[0][pos12.first][pos12.second]<0) && SegmentGraph.vEdges[i].Head1==(Components[0][pos11.first][pos11.second]>0))
				flag1=true;

			if(pos21.first==pos22.first && pos21.second<pos22.second && SegmentGraph.vEdges[i].Head1==(Components[1][pos21.first][pos21.second]<0) && SegmentGraph.vEdges[i].Head2==(Components[1][pos22.first][pos22.second]>0))
				flag2=true;
			else if(pos21.first==pos22.first && pos21.second>pos22.second && SegmentGraph.vEdges[i].Head2==(Components[1][pos22.first][pos22.second]<0) && SegmentGraph.vEdges[i].Head1==(Components[1][pos21.first][pos21.second]>0))
				flag2=true;

			if(flag1 || flag2){
				// check whether prediction is autosomal or on X Y. remove contig/mitochrondia predictions
				string chrname1=RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr];
				string chrname2=RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr];
				// if(chrname1.substr(0,3)=="chr")
				// 	chrname1=chrname1.substr(3);
				// if(chrname2.substr(0,3)=="chr")
				// 	chrname2=chrname2.substr(3);
				// if((chrname1[0]<'0' || chrname1[0]>'9') && chrname1[0]!='X' && chrname1[0]!='Y')
				// 	continue;
				// else if((chrname2[0]<'0' || chrname2[0]>'9') && chrname2[0]!='X' && chrname2[0]!='Y')
				// 	continue;
				map<Edge_t, vector< pair<int,int> > >::iterator itmap=ExactBP.find(SegmentGraph.vEdges[i]);
				vector< pair<int,int> > BP;
				map<Edge_t, vector< pair<int,int> > >::const_iterator itsup = ExactBP_concord_support.find(SegmentGraph.vEdges[i]);
				if(itsup == ExactBP_concord_support.cend())
					cout<<"Error: "<<i<<"\t"<<(itmap->first).Ind1<<"\t"<<(itmap->first).Ind2<<endl;
				assert(itsup != ExactBP_concord_support.cend());
				const vector< pair<int,int> >& Support = itsup->second;
				if(itmap==ExactBP.end() || itmap->second.size()==0){
					int bp1,bp2;
					bp1=(SegmentGraph.vEdges[i].Head1?SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position:(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length));
					bp2=(SegmentGraph.vEdges[i].Head2?SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position:(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length));
					BP.push_back(make_pair(bp1,bp2));
				}
				else{
					BP=itmap->second;
				}
				if(BP.size() != Support.size())
					cout<<"Error: number of breakpoints in BPs is different from number of breakpoints in Support.\n";
				assert(BP.size() == Support.size());
				for(int k=0; k<BP.size(); k++){
					if(SegmentGraph.vEdges[i].Head1){
						output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr]<<'\t';
						output<<BP[k].first<<'\t';
						output<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length)<<'\t';
					}
					else{
						output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr]<<'\t';
						output<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position<<'\t';
						output<<BP[k].first<<'\t';
					}
					if(SegmentGraph.vEdges[i].Head2){
						output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr]<<'\t';
						output<<BP[k].second<<'\t';
						output<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length)<<'\t';
					}
					else{
						output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr]<<'\t';
						output<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position<<'\t';
						output<<BP[k].second<<'\t';
					}
					output<<".\t"<<SegmentGraph.vEdges[i].Weight<<"\t";
					output<<(SegmentGraph.vEdges[i].Head1?"-\t":"+\t");
					output<<(SegmentGraph.vEdges[i].Head2?"-\t":"+\t");
					output<<SegmentGraph.vEdges[i].Ind1<<"\t"<<SegmentGraph.vEdges[i].Ind2<<"\t";
					output<<(Support[k].first)<<"\t"<<(Support[k].second)<<endl;
				}
			}
		}	
	}
	output.close();

	// ofstream output2(outputfile2, ios::out);
	// output2<<"# chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tnum_concordantfrag_bp1\tnum_concordantfrag_bp2\n";
	// for(int i=0; i<SegmentGraph.vEdges.size(); i++){
	// 	int ind1=SegmentGraph.vEdges[i].Ind1, ind2=SegmentGraph.vEdges[i].Ind2;
	// 	bool flag_chr=(SegmentGraph.vNodes[ind1].Chr==SegmentGraph.vNodes[ind2].Chr);
	// 	bool flag_ori=(SegmentGraph.vEdges[i].Head1==false && SegmentGraph.vEdges[i].Head2==true);
	// 	bool flag_dist=(SegmentGraph.vNodes[ind2].Position-SegmentGraph.vNodes[ind1].Position-SegmentGraph.vNodes[ind1].Length<=Concord_Dist_Pos || ind2-ind1<=Concord_Dist_Idx);
	// 	if(!flag_chr || !flag_ori || !flag_dist){
	// 		pair<int,int> pos1=Node_NewChr[1][SegmentGraph.vEdges[i].Ind1];
	// 		pair<int,int> pos22=Node_NewChr[1][SegmentGraph.vEdges[i].Ind2];
	// 		bool flag=false;  // if rearrangement made new edges concordant with the original genome
	// 		if(pos1.first==pos12.first && pos1.second<pos12.second && SegmentGraph.vEdges[i].Head1==(Components[pos1.first][1][pos1.second]<0) && SegmentGraph.vEdges[i].Head2==(Components[pos12.first][1][pos12.second]>0))
	// 			flag=true;
	// 		else if(pos1.first==pos12.first && pos1.second>pos12.second && SegmentGraph.vEdges[i].Head2==(Components[pos12.first][1][pos12.second]<0) && SegmentGraph.vEdges[i].Head1==(Components[pos1.first][1][pos1.second]>0))
	// 			flag=true;
	// 		if(flag){
	// 			// chech whether prediction is autosomal or on X Y. remove contig/mitochrondia predictions
	// 			string chrname1=RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr];
	// 			string chrname2=RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr];
	// 			// if(chrname1.substr(0,3)=="chr")
	// 			// 	chrname1=chrname1.substr(3);
	// 			// if(chrname2.substr(0,3)=="chr")
	// 			// 	chrname2=chrname2.substr(3);
	// 			// if((chrname1[0]<'0' || chrname1[0]>'9') && chrname1[0]!='X' && chrname1[0]!='Y')
	// 			// 	continue;
	// 			// else if((chrname2[0]<'0' || chrname2[0]>'9') && chrname2[0]!='X' && chrname2[0]!='Y')
	// 			// 	continue;
	// 			map<Edge_t, vector< pair<int,int> > >::iterator itmap=ExactBP.find(SegmentGraph.vEdges[i]);
	// 			vector< pair<int,int> > BP;
	// 			map<Edge_t, vector< pair<int,int> > >::const_iterator itsup = ExactBP_concord_support.find(SegmentGraph.vEdges[i]);
	// 			if(itsup == ExactBP_concord_support.cend())
	// 				cout<<"Error: "<<i<<"\t"<<(itmap->first).Ind1<<"\t"<<(itmap->first).Ind2<<endl;
	// 			assert(itsup != ExactBP_concord_support.cend());
	// 			const vector< pair<int,int> >& Support = itsup->second;
	// 			if(itmap==ExactBP.end() || itmap->second.size()==0){
	// 				int bp1,bp2;
	// 				bp1=(SegmentGraph.vEdges[i].Head1?SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position:(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length));
	// 				bp2=(SegmentGraph.vEdges[i].Head2?SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position:(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length));
	// 				BP.push_back(make_pair(bp1,bp2));
	// 			}
	// 			else
	// 				BP=itmap->second;
	// 			if(BP.size() != Support.size())
	// 				cout<<"Error: number of breakpoints in BPs is different from number of breakpoints in Support.\n";
	// 			assert(BP.size() == Support.size());
	// 			for(int k=0; k<BP.size(); k++){
	// 				if(SegmentGraph.vEdges[i].Head1){
	// 					output2<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr]<<'\t';
	// 					output2<<BP[k].first<<'\t';
	// 					output2<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length)<<'\t';
	// 				}
	// 				else{
	// 					output2<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr]<<'\t';
	// 					output2<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position<<'\t';
	// 					output2<<BP[k].first<<'\t';
	// 				}
	// 				if(SegmentGraph.vEdges[i].Head2){
	// 					output2<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr]<<'\t';
	// 					output2<<BP[k].second<<'\t';
	// 					output2<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length)<<'\t';
	// 				}
	// 				else{
	// 					output2<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr]<<'\t';
	// 					output2<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position<<'\t';
	// 					output2<<BP[k].second<<'\t';
	// 				}
	// 				output2<<".\t"<<SegmentGraph.vEdges[i].Weight<<"\t";
	// 				output2<<(SegmentGraph.vEdges[i].Head1?"-\t":"+\t");
	// 				output2<<(SegmentGraph.vEdges[i].Head2?"-\t":"+\t");
	// 				output2<<(Support[k].first)<<"\t"<<(Support[k].second)<<endl;
	// 			}
	// 		}
	// 	}
	// }
	// output2.close();
};

void TmpWriteBEDPE(string outputfile, SegmentGraph_t& SegmentGraph, vector<string>& RefName){
	ofstream output(outputfile, ios::out);
	output<<"# chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\n";
	for(int i=0; i<SegmentGraph.vEdges.size(); i++){
		int ind1=SegmentGraph.vEdges[i].Ind1, ind2=SegmentGraph.vEdges[i].Ind2;
		bool flag_chr=(SegmentGraph.vNodes[ind1].Chr==SegmentGraph.vNodes[ind2].Chr);
		bool flag_ori=(SegmentGraph.vEdges[i].Head1==false && SegmentGraph.vEdges[i].Head2==true);
		bool flag_dist=(SegmentGraph.vNodes[ind2].Position-SegmentGraph.vNodes[ind1].Position-SegmentGraph.vNodes[ind1].Length<=Concord_Dist_Pos || ind2-ind1<=Concord_Dist_Idx);
		if(!flag_chr || !flag_ori || !flag_dist){
			// chech whether prediction is autosomal or on X Y. remove contig/mitochrondia predictions
			vector< pair<int,int> > BP;
			int bp1,bp2;
			bp1=(SegmentGraph.vEdges[i].Head1?SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position:(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length));
			bp2=(SegmentGraph.vEdges[i].Head2?SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position:(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length));
			BP.push_back(make_pair(bp1,bp2));
			for(int k=0; k<BP.size(); k++){
				if(SegmentGraph.vEdges[i].Head1){
					output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr]<<'\t';
					output<<BP[k].first<<'\t';
					output<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Length)<<'\t';
				}
				else{
					output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Chr]<<'\t';
					output<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind1].Position<<'\t';
					output<<BP[k].first<<'\t';
				}
				if(SegmentGraph.vEdges[i].Head2){
					output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr]<<'\t';
					output<<BP[k].second<<'\t';
					output<<(SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position+SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Length)<<'\t';
				}
				else{
					output<<RefName[SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Chr]<<'\t';
					output<<SegmentGraph.vNodes[SegmentGraph.vEdges[i].Ind2].Position<<'\t';
					output<<BP[k].second<<'\t';
				}
				output<<".\t"<<SegmentGraph.vEdges[i].Weight<<"\t";
				output<<(SegmentGraph.vEdges[i].Head1?"-\t":"+\t");
				output<<(SegmentGraph.vEdges[i].Head2?"-\t":"+\t");
				output<<endl;
			}
		}
	}
	output.close();
};

void OutputNewGenome(SegmentGraph_t& SegmentGraph, vector< vector< vector<int> > >& Components, const vector<string>& RefSequence, const vector<string>& RefName, vector<string> fnames){
	// ofstream output1(outputfile1, ios::out);
	// ofstream output2(outputfile2, ios::out);
	// vector<ofstream> outputVec;
	// outputVec.push_back(output1);
	// outputVec.push_back(output2);

	for (int m = 0;m<fnames.size();m++){
		ofstream output(fnames[m], ios::out);
		for(int i=0; i<Components[m].size(); i++){
			string info="PA:", seq, tmpseq;
			for(int j=0; j<Components[m][i].size(); j++){
				int k;
				for(k=j+1; k<Components[m][i].size() && Components[m][i][k]-Components[m][i][k-1]==1 && SegmentGraph.vNodes[abs(Components[m][i][j])-1].Chr==SegmentGraph.vNodes[abs(Components[m][i][k])-1].Chr; k++){}
				if(Components[m][i][j]>0){
					int curChr=SegmentGraph.vNodes[abs(Components[m][i][j])-1].Chr;
					int curStart=SegmentGraph.vNodes[abs(Components[m][i][j])-1].Position;
					int curLen=SegmentGraph.vNodes[abs(Components[m][i][k-1])-1].Position+SegmentGraph.vNodes[abs(Components[m][i][k-1])-1].Length-SegmentGraph.vNodes[abs(Components[m][i][j])-1].Position;
					tmpseq=RefSequence[curChr].substr(curStart, curLen);
					info+="{"+RefName[curChr]+","+to_string(curStart)+","+to_string(curLen)+"}";
				}
				else{
					int curChr=SegmentGraph.vNodes[abs(Components[m][i][k-1])-1].Chr;
					int curStart=SegmentGraph.vNodes[abs(Components[m][i][k-1])-1].Position;
					int curLen=SegmentGraph.vNodes[abs(Components[m][i][j])-1].Position+SegmentGraph.vNodes[abs(Components[m][i][j])-1].Length-SegmentGraph.vNodes[abs(Components[m][i][k-1])-1].Position;
					tmpseq=RefSequence[curChr].substr(curStart, curLen);
					info+="{"+RefName[curChr]+","+to_string(curStart)+","+to_string(curLen)+"}";
				}
				if(Components[m][i][j]<0)
					ReverseComplement(tmpseq.begin(), tmpseq.end());
				seq+=tmpseq;
				info+=((Components[m][i][j]<0)?"R-":"F-");
				j=k-1;
			}
			info=info.substr(0, info.size()-1);
			output<<">chr"<<(i+1)<<'\t'<<"LN:"<<seq.size()<<'\t'<<info<<endl;
			int idx=0;
			while(idx<seq.size()){
				int nextidx=min(idx+80, (int)seq.size());
				output<<seq.substr(idx, nextidx-idx)<<endl;
				idx=nextidx;
			}
		}
		output.close();
	}
	
	
	// ofstream output1(outputfile1, ios::out);
	// ofstream output2(outputfile2, ios::out);

	// for(int i=0; i<Components[0].size(); i++){
	// 	string info="PA:", seq, tmpseq;
	// 	for(int j=0; j<Components[i][0].size(); j++){
	// 		int k;
	// 		for(k=j+1; k<Components[i][0].size() && Components[i][0][k]-Components[i][0][k-1]==1 && SegmentGraph.vNodes[abs(Components[i][0][j])-1].Chr==SegmentGraph.vNodes[abs(Components[i][0][k])-1].Chr; k++){}
	// 		if(Components[i][0][j]>0){
	// 			int curChr=SegmentGraph.vNodes[abs(Components[i][0][j])-1].Chr;
	// 			int curStart=SegmentGraph.vNodes[abs(Components[i][0][j])-1].Position;
	// 			int curLen=SegmentGraph.vNodes[abs(Components[i][0][k-1])-1].Position+SegmentGraph.vNodes[abs(Components[i][0][k-1])-1].Length-SegmentGraph.vNodes[abs(Components[i][0][j])-1].Position;
	// 			tmpseq=RefSequence[curChr].substr(curStart, curLen);
	// 			info+="{"+RefName[curChr]+","+to_string(curStart)+","+to_string(curLen)+"}";
	// 		}
	// 		else{
	// 			int curChr=SegmentGraph.vNodes[abs(Components[i][0][k-1])-1].Chr;
	// 			int curStart=SegmentGraph.vNodes[abs(Components[i][0][k-1])-1].Position;
	// 			int curLen=SegmentGraph.vNodes[abs(Components[i][0][j])-1].Position+SegmentGraph.vNodes[abs(Components[i][0][j])-1].Length-SegmentGraph.vNodes[abs(Components[i][0][k-1])-1].Position;
	// 			tmpseq=RefSequence[curChr].substr(curStart, curLen);
	// 			info+="{"+RefName[curChr]+","+to_string(curStart)+","+to_string(curLen)+"}";
	// 		}
	// 		if(Components[i][0][j]<0)
	// 			ReverseComplement(tmpseq.begin(), tmpseq.end());
	// 		seq+=tmpseq;
	// 		info+=((Components[i][0][j]<0)?"R-":"F-");
	// 		j=k-1;
	// 	}
	// 	info=info.substr(0, info.size()-1);
	// 	output1<<">chr"<<(i+1)<<'\t'<<"LN:"<<seq.size()<<'\t'<<info<<endl;
	// 	int idx=0;
	// 	while(idx<seq.size()){
	// 		int nextidx=min(idx+80, (int)seq.size());
	// 		output1<<seq.substr(idx, nextidx-idx)<<endl;
	// 		idx=nextidx;
	// 	}

	// 	for(int j=0; j<Components[i][1].size(); j++){
	// 		int k;
	// 		for(k=j+1; k<Components[i][1].size() && Components[i][1][k]-Components[i][1][k-1]==1 && SegmentGraph.vNodes[abs(Components[i][1][j])-1].Chr==SegmentGraph.vNodes[abs(Components[i][1][k])-1].Chr; k++){}
	// 		if(Components[i][1][j]>0){
	// 			int curChr=SegmentGraph.vNodes[abs(Components[i][1][j])-1].Chr;
	// 			int curStart=SegmentGraph.vNodes[abs(Components[i][1][j])-1].Position;
	// 			int curLen=SegmentGraph.vNodes[abs(Components[i][1][k-1])-1].Position+SegmentGraph.vNodes[abs(Components[i][1][k-1])-1].Length-SegmentGraph.vNodes[abs(Components[i][1][j])-1].Position;
	// 			tmpseq=RefSequence[curChr].substr(curStart, curLen);
	// 			info+="{"+RefName[curChr]+","+to_string(curStart)+","+to_string(curLen)+"}";
	// 		}
	// 		else{
	// 			int curChr=SegmentGraph.vNodes[abs(Components[i][1][k-1])-1].Chr;
	// 			int curStart=SegmentGraph.vNodes[abs(Components[i][1][k-1])-1].Position;
	// 			int curLen=SegmentGraph.vNodes[abs(Components[i][1][j])-1].Position+SegmentGraph.vNodes[abs(Components[i][1][j])-1].Length-SegmentGraph.vNodes[abs(Components[i][1][k-1])-1].Position;
	// 			tmpseq=RefSequence[curChr].substr(curStart, curLen);
	// 			info+="{"+RefName[curChr]+","+to_string(curStart)+","+to_string(curLen)+"}";
	// 		}
	// 		if(Components[i][1][j]<0)
	// 			ReverseComplement(tmpseq.begin(), tmpseq.end());
	// 		seq+=tmpseq;
	// 		info+=((Components[i][1][j]<0)?"R-":"F-");
	// 		j=k-1;
	// 	}
	// 	info=info.substr(0, info.size()-1);
	// 	output2<<">chr"<<(i+1)<<'\t'<<"LN:"<<seq.size()<<'\t'<<info<<endl;
	// 	idx=0;
	// 	while(idx<seq.size()){
	// 		int nextidx=min(idx+80, (int)seq.size());
	// 		output2<<seq.substr(idx, nextidx-idx)<<endl;
	// 		idx=nextidx;
	// 	}
	// }
	// output1.close();
	// output2.close();
};
