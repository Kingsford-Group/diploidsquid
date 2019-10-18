#include "WriteIO.h"
#include "SegmentGraph.h"
#include "Config.h"
#include <stdlib.h>

using namespace std;

void testWriteComponent(){
    vector< vector< vector<int> > > components; 
    components.resize(2);

    for (int k=0; k<2; k++)
        for (int i=0; i<3; i++){
            vector<int> comp;
            for (int j=0;j<5;j++){
                int a = j*i;
                if (k==1)
                    a = -a;
                comp.push_back(a);
            }
            components[k].push_back(comp);
        }
    WriteComponents("testWriteComp.txt", components);
    SegmentGraph_t graph;
    components[0] = graph.SortComponents(components[0]);
    components[1] = graph.SortComponents(components[1]);
    WriteComponents("testWriteComp_sorted.txt", components);
}

void testGraph(string filename){
    SegmentGraph_t graph(filename) ;

    graph.OutputGraph("testGraph1_out.txt");
}

// void testOrdering(string filename){
//     SegmentGraph_t graph(filename);
//     graph.OutputGraph("testGraph1_out.txt");
//     vector< vector < vector<int> > > Components = graph.Ordering();
//     WriteComponents("testOrdering1_out.txt", Components);

// }

void testFindEdge(SegmentGraph_t SegmentGraph, vector< vector < vector<int> > > components){
    for(int i=0; i < components.size(); i ++){
        for(int j=0; j<components[i].size(); j++){
            for(int k=0; k<components[i][j].size()-1;k++){
                // cout << "========================" << endl;
                // cout << abs(components[i][j][k])-1 << "\t" << abs(components[i][j][k+1])-1 << endl;
                for (Edge_t e : SegmentGraph.findEdges(abs(components[i][j][k])-1,abs(components[i][j][k+1])-1)){
                    cout << e.Ind1 << "\t" <<  e.Ind2 << endl;
                }
            }
        }
    }
}

void testTestConcordant(SegmentGraph_t graph, vector<vector<vector<int>>> Components, string prefix, string outdir){
    int c = 0;
    int all_edge = 0;
    int c_D = 0;
    vector<Edge_t> discordantEdges = graph.getDiscordant();

    ofstream output1(outdir+"/ConcordantEdges_diploid_"+prefix+".txt", ios::out);
	output1<<"ind1\tind2\n";

    ofstream output2(outdir+"/C_DiscordantEdges_diploid_"+prefix+".txt", ios::out);
	output2<<"ind1\tind2\n";

	// For each node, store its component and index in component at its index
	vector< vector< pair<int, int> > > Node_NewChr; Node_NewChr.resize(2);
	Node_NewChr[0].resize(graph.vNodes.size());
	Node_NewChr[1].resize(graph.vNodes.size());
	for(unsigned int i=0; i<Components.size(); i++)
		for (unsigned int k=0; k<Components[i].size();k++)
			for(unsigned int j=0; j<Components[i][k].size(); j++)
				Node_NewChr[i][abs(Components[i][k][j])-1]=make_pair(k, j); 

    for (int i=0; i<graph.vEdges.size(); i++){
		int ind1=graph.vEdges[i].Ind1, ind2=graph.vEdges[i].Ind2;
        pair<int,int> pos11=Node_NewChr[0][graph.vEdges[i].Ind1]; //pos -> (comp, position)
        pair<int,int> pos12=Node_NewChr[0][graph.vEdges[i].Ind2];

        pair<int,int> pos21=Node_NewChr[1][graph.vEdges[i].Ind1];
        pair<int,int> pos22=Node_NewChr[1][graph.vEdges[i].Ind2];


        // cout << ind1 << "\t" << ind2 << endl;
        // cout << pos11.second << "\t" << pos12.second <<endl;
        // cout << "-------------" <<endl;

        // Record BP if edge is made concordant after either arrangement
        bool flag1=false;
        bool flag2=false;  
        if(pos11.first==pos12.first && pos11.second<pos12.second && graph.vEdges[i].Head1==(Components[0][pos11.first][pos11.second]<0) && graph.vEdges[i].Head2==(Components[0][pos12.first][pos12.second]>0))
            flag1=true;
        else if(pos11.first==pos12.first && pos11.second>pos12.second && graph.vEdges[i].Head2==(Components[0][pos12.first][pos12.second]<0) && graph.vEdges[i].Head1==(Components[0][pos11.first][pos11.second]>0))
            flag1=true;

        if(pos21.first==pos22.first && pos21.second<pos22.second && graph.vEdges[i].Head1==(Components[1][pos21.first][pos21.second]<0) && graph.vEdges[i].Head2==(Components[1][pos22.first][pos22.second]>0))
            flag2=true;
        else if(pos21.first==pos22.first && pos21.second>pos22.second && graph.vEdges[i].Head2==(Components[1][pos22.first][pos22.second]<0) && graph.vEdges[i].Head1==(Components[1][pos21.first][pos21.second]>0))
            flag2=true;

        if (flag1 || flag2){
            output1 << Components[0][pos11.first][pos11.second] << "\t" << Components[0][pos12.first][pos12.second] <<endl;
            c++;
        }


        bool flag_chr=(graph.vNodes[ind1].Chr==graph.vNodes[ind2].Chr);
		bool flag_ori=(graph.vEdges[i].Head1==false && graph.vEdges[i].Head2==true);
        // cout << i << "\t" << graph.vEdges[i].Head1 << graph.vEdges[i].Head2 << endl;
		// bool flag_dist=(graph.vNodes[ind2].Position-graph.vNodes[ind1].Position-graph.vNodes[ind1].Length<=100 || ind2-ind1<=Concord_Dist_Idx);
		// if((!flag_chr || !flag_ori ) && (flag1)){
        //     c_D ++;
        //     output2 << Components[0][pos11.first][pos11.second]<< "\t" << Components[0][pos12.first][pos12.second] <<endl;
        // }
        // if((!flag_chr || !flag_ori ) && (flag2)){
        //     c_D ++;
        //     output2 << Components[0][pos21.first][pos21.second] << "\t" << Components[0][pos22.first][pos22.second] <<endl;
        // }
        if(!flag_chr || !flag_ori ){
            if (flag1 && flag2){
                c_D ++;
                output2 << Components[0][pos11.first][pos11.second] << "\t" << Components[0][pos12.first][pos12.second] <<endl;
            }
            else if (flag1){
                c_D ++;
                output2 << Components[0][pos11.first][pos11.second] << "\t" << Components[0][pos12.first][pos12.second] <<endl;
            }
            else if (flag2){
                c_D ++;
                output2 << Components[1][pos21.first][pos21.second] << "\t" << Components[1][pos22.first][pos22.second] <<endl;
            }
        }
    }

    cout << "All discordant edges: " << discordantEdges.size() << endl;
    // cout << "All edges by diploid_squid: " << all_edge << endl;
    cout << "Conc edges: " << c << endl;
    cout << "Conc_D edges: " << c_D << endl;
    output1.close();
    output2.close();
}

int main(int argc, char * argv[]){
    // cout << "Testing WriteComponent..." <<endl;
    // testWriteComponent();
    // cout << "Testing read graph.." << endl;
    // testGraph("testGraph1.txt");
    // cout << "Testing ordering on graph1" << endl;
    // testOrdering(argv[1]);
    
    int Mode = DSQUID;
    if (argc != 5){
        cout << "This is a test function that allows you to use a GSG as input. This function outputs the concordant edges in the arranged GSG and concordant edges that were discordant in the arranged GSG."
        cout << "USAGE: test <-I/-D/-A> <graph_file> <prefix> <outdir>\n";
        return -1;
    }
    if(string(argv[1]) == "-I"){
        Mode=ISQUID;
    }
    else if (string(argv[1]) == "-D"){
        Mode=DSQUID;
    }
    else if (string(argv[1]) == "-A"){
        Mode=APPROX;
    }
    else if (string(argv[1]) == "-A2"){
        Mode=APPROX_2;
    }
    else{
        cout << "USAGE: test <-I/-D/-A> <graph_file> <prefix> <outdir>\n";
        return -1;
    }
    SegmentGraph_t graph(argv[2]);
    string prefix=argv[3];
    string outdir=argv[4];
    // graph.OutputGraph("testGraph1_out.txt");
    cout << "ordering" <<endl;
    vector< vector < vector<int> > > Components = graph.Ordering(Mode, outdir+"/"+prefix);

    vector<vector<int>> Components1, Components2;
    vector<vector<vector<int>>> reorderComponents;
    for (int i = 0 ;i<Components.size();i++){
        Components1.push_back(Components[i][0]);
        Components2.push_back(Components[i][1]);
    }
    reorderComponents.push_back(Components1);
    reorderComponents.push_back(Components2);

    // WriteComponents(outdir+"/"+prefix+"_component_pri.txt", reorderComponents);
    
    for (int i = 0;i<reorderComponents.size();i++){
    //     reorderComponents[i]=graph.SortComponents(reorderComponents[i]);
    //     reorderComponents[i]=graph.MergeSingleton(reorderComponents[i], RefLength);
         reorderComponents[i]=graph.SortComponents(reorderComponents[i]);
    //     reorderComponents[i]=graph.MergeComponents(reorderComponents[i]);
    }

    // WriteComponents(outdir+"/"+prefix+"_component_sorted.txt", reorderComponents);

    // WriteComponents("testOrdering_unsort.txt",Components);
    // Components = graph.SortComponents(Components);
    // WriteComponents("testOrdering_out.txt", Components);
    // testFindEdge(graph, Components);
    testTestConcordant(graph, reorderComponents, prefix, outdir);
    
    // cout << "Testing ordering on graph2" << endl;
    // testOrdering("testGraph2.txt");
}
