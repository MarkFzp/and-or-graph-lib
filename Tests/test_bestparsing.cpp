//
// Created by Hao Wang on 18/4/11.
//

#include <vector>
#include <string>

#include "../T-AOG/T_AOG.h"
#include "../Learner/Online_Learner.h"

using namespace std;
using namespace AOG_LIB;

int SAMPLING = 200;

void OutputToFile(string filename, const T_AOG<string>& t_aog,VertexId sample_root)
{
    //open file
    ofstream outFile;
    outFile.open(filename);   
    if(outFile.is_open())
    {
        cerr<<"inside opening file!\n";

        while (SAMPLING)
        {
            cerr<<"inside sampling!\n";

            vector<Symbolic_State<string> > res;
            cerr<<"before sampling!\n";
            
            t_aog.Sample(sample_root, res);
            SAMPLING--;
            for(auto state : res)
            {
                cerr<<"inside getcontent!\n";
                
                outFile << state.GetContent()<<" ";
            
            }
            outFile<<"\n";
        }
    }
    else
    {
        cerr<<"Error opening file!\n";
        throw exception();
    }
    outFile.close();
}

int main()
{
    string one_hundred = "100";
    string two = "2";
    string fifty = "50";
    string four = "4";
    string twenty_five = "25";
    string ten = "10";
    string five = "5";
    Symbolic_State<string> state_100(one_hundred, false);    
    Symbolic_State<string> state_2(two, true);    
    Symbolic_State<string> state_50(fifty, false);    
    Symbolic_State<string> state_4(four, true);    
    Symbolic_State<string> state_25(twenty_five, true);    
    Symbolic_State<string> state_10(ten, false);    
    Symbolic_State<string> state_5(five, true);   
    SequenceType<string> state_2_5 = {state_2, state_5};     
    SequenceType<string> state_2_50 = {state_2, state_50};
    SequenceType<string> state_10_10 = {state_10, state_10};
    SequenceType<string> state_4_25 = {state_4, state_25};
    SequenceType<string> state_5_10 = {state_5, state_10};
    SequenceType<string> state_2_25 = {state_2, state_25};
    Symbolic_Rule<string> hundred_2_50 = {state_100, state_2_50};
    Symbolic_Rule<string> hundred_4_25 = {state_100, state_4_25};
    Symbolic_Rule<string> hundred_10_10 = {state_100, state_10_10};
    Symbolic_Rule<string> fifty_5_10 = {state_50, state_5_10};
    Symbolic_Rule<string> fifty_2_25 = {state_50, state_2_25};
    Symbolic_Rule<string> ten_2_5 = {state_10, state_2_5};
    
    vector<Symbolic_Rule<string> > rules = {hundred_2_50, hundred_4_25, hundred_10_10, fifty_5_10, fifty_2_25, ten_2_5};


    shared_ptr<T_AOG<string> > related_graph = make_shared<T_AOG<string> >(rules);
    VertexId root_id;
    try
    {
        root_id = related_graph->GetRoot();
    }
    catch(...)
    {
        Symbolic_State<string> source_state;
        std::shared_ptr<AOG_Vertex<string> > source_ptr(new AOG_Vertex<string>(source_state, false, true));
        related_graph->AddVertex(source_ptr);
    }

    // change weight
    VertexId hundred_id = related_graph->GetVertexIdByState(state_100);
    unordered_map<VertexId, double> weights = related_graph->GetOutEdgeWeights(hundred_id, false);
    for (auto &weight: weights)
    {
        vector<VertexId> children = related_graph->ChildrenVertices(weight.first);
        // change branch 100 -> 2, 50 to 1.5
        if (children.size() == 2 && related_graph->GetStateByVertexId(children[0]) == state_2 && related_graph->GetStateByVertexId(children[1]) == state_50)
        {
            weight.second += 0.5;
        }
        // change branch 100 -> 4, 25 to 2        
        if (children.size() == 2 && related_graph->GetStateByVertexId(children[0]) == state_4 && related_graph->GetStateByVertexId(children[1]) == state_25)
        {
            weight.second += 1;
        }
    }
    related_graph->SetOutEdgeWeights(hundred_id, weights);

    VertexId fifty_id = related_graph->GetVertexIdByState(state_50);
    weights = related_graph->GetOutEdgeWeights(fifty_id, false);
    for (auto &weight: weights)
    {
        vector<VertexId> children = related_graph->ChildrenVertices(weight.first);
        // change branch 50 -> 5, 10 to 0.5        
        if (children.size() == 2 && related_graph->GetStateByVertexId(children[0]) == state_5 && related_graph->GetStateByVertexId(children[1]) == state_10)
        {
            weight.second -= 0.5;
        }
        // change branch 50 -> 2, 25 to 0.8        
        if (children.size() == 2 && related_graph->GetStateByVertexId(children[0]) == state_2 && related_graph->GetStateByVertexId(children[1]) == state_25)
        {
            weight.second -= 0.2;
        }

    }
    related_graph->SetOutEdgeWeights(fifty_id, weights);


    AOG_LIB::Online_Learner<string> online_Learner(related_graph);

    vector<Symbolic_State<string> > sentence = {state_2, state_5, state_2, state_5};

    cout << "Root is: " << related_graph->GetRoot() << endl;

    OutputToFile("./Grammar_Example/Output/test_bestparsing_output.txt",*related_graph, related_graph->GetVertexIdByState(state_100));

    cerr << "----------------DEBUGGING STARTS HERE------------------" << endl;
    online_Learner.MergeNewParsingWithGraph(sentence);

	return 0;
}