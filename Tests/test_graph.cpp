//
// Created by Luyao Yuan on 17/10/24.
//

#include <string>
#include "../T-AOG/T_AOG.h"
#include "../Learner/Online_Learner.h"

using namespace std;
using namespace AOG_LIB;


int main()
{
    //states for rules of square
    string square_str = "square";
    string square_adjcnt_edges_str = "square adjacent edges";
    string square_sprt_edges_str = "square separate edges";
    string atom_edge_str = "edge";

    Symbolic_State<string> square(square_str,false);
    Symbolic_State<string> sqr_adjcnt_edges(square_adjcnt_edges_str,false);
    Symbolic_State<string> sqr_sprt_edges(square_sprt_edges_str,false);
    Symbolic_State<string> edge(atom_edge_str,true);
    vector<Symbolic_State<string> > sqr_rs_1={sqr_adjcnt_edges,sqr_adjcnt_edges};
    vector<Symbolic_State<string> > sqr_rs_2={sqr_sprt_edges,sqr_sprt_edges};
    vector<Symbolic_State<string> > edge_rs={edge,edge};
    
    //rules of square
    Symbolic_Rule<string> sqr_to_edges1(square,sqr_rs_1);
    Symbolic_Rule<string> sqr_to_edges2(square,sqr_rs_2);
    Symbolic_Rule<string> sqr_edges_to_edge1(sqr_adjcnt_edges,edge_rs);
    Symbolic_Rule<string> sqr_edges_to_edge2(sqr_sprt_edges,edge_rs);


    //states for rules of triangle
    string tri_str = "triangle";
    string tri_adjcnt_edges_str = "triangle adjacent edges";
    Symbolic_State<string> triangle(tri_str,false);
    Symbolic_State<string> tri_adjcnt_edges(tri_adjcnt_edges_str,false);
    vector<Symbolic_State<string> > tri_rs_1 = {tri_adjcnt_edges,edge};

    //rules of triangle
    Symbolic_Rule<string> tri_to_parts(triangle,tri_rs_1);
    Symbolic_Rule<string> tri_edges_to_edge(tri_adjcnt_edges,edge_rs);

    //states for parallelogram
    string parallelogram_str = "parallelogram";
    Symbolic_State<string> parallelogram(parallelogram_str,false);
    vector<Symbolic_State<string> > prll_rs = {tri_adjcnt_edges,tri_adjcnt_edges};    
    //rules for parallelogram
    Symbolic_Rule<string> parallel_to_tri(parallelogram,prll_rs);

    //construct T_AOG using rules defined above
    vector<Symbolic_Rule<string> > rules;
    
    rules.push_back(sqr_to_edges1);
    rules.push_back(sqr_to_edges2);
    rules.push_back(sqr_edges_to_edge1);
    rules.push_back(sqr_edges_to_edge2);
    rules.push_back(tri_to_parts);
    rules.push_back(tri_edges_to_edge);
    rules.push_back(parallel_to_tri);
    
    shared_ptr<T_AOG<string> > related_graph = make_shared<T_AOG<string> >();
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
    
    Symbolic_State<string> root = related_graph->GetStateByVertexId(related_graph->GetRoot());
    vector<Symbolic_State<string> > next_level = {triangle, parallelogram, square};
    Symbolic_Rule<string> root_to_next = {root, next_level};

    related_graph->AddRule(sqr_to_edges1);
    related_graph->AddRule(sqr_to_edges2);
    related_graph->AddRule(sqr_edges_to_edge1);
    related_graph->AddRule(sqr_edges_to_edge2);
    related_graph->AddRule(tri_to_parts);
    related_graph->AddRule(tri_edges_to_edge);
    related_graph->AddRule(parallel_to_tri);
    related_graph->AddRule(root_to_next);
    
    
    related_graph->OutputGraph("visualize","/Users/jingyue/Desktop",false);
    
    std::shared_ptr<T_AOG<string> > new_graph = std::make_shared<T_AOG<string> >(*related_graph);
    auto test_root_id = new_graph->GetRoot();
    auto test_root_state = new_graph->GetStateByVertexId(test_root_id);
    Symbolic_State<string> test_1("hhhhhhhhhh",true);
    Symbolic_State<string> test_2("llllllllll",true);
    vector<Symbolic_State<string> > test_targets{test_1,test_2};
    Symbolic_Rule<string> test_rule(test_root_state,test_targets);
    new_graph->AddRule(test_rule);
   
    cerr<<"new graph visualization\n";
    new_graph->OutputGraph("visualize_new","/Users/jingyue/Desktop",false);
    cerr<<"original graph visualization\n";
    related_graph->OutputGraph("visualize_old","/Users/jingyue/Desktop",false);
    
    
    // cout<<"The number of rules in this T_AOG: "<<t_aog.NumOfRules()<<endl
    // <<"The number of states in this T_AOG: "<<t_aog.NumOfStates()<<endl
    // 	<< "The number of leaf states in this T_AOG: " << t_aog.NumOfLeafStates()<<endl
	// <<"Number of Vertex in this T_AOG: " <<t_aog.NumberOfVertices()<<endl<<endl<<endl;
    
    
    
    // //sample a parse tree from the T_AOG
    // VertexId sample_root = 5;
    // cout<<"Sampling a parse tree from node "<<sample_root<<" of the T_AOG Graph..."<<endl;
    // vector<Symbolic_State<string>> data_seq;
    // auto parse_tree = t_aog.Sample(sample_root,data_seq);

    // for(auto iter : parse_tree)
    // {
    
    //     cout<< "parent node: "<<iter.first<<endl;
    //     cout<<"children nodes: ";
    //     for(auto child : iter.second)
    //         cout<< child<<" ";
    //     cout<<endl;
    // }
    // cout<<endl<<endl;
   // //states for rules of square
   //  string square_str = "square";
   //  string square_adjcnt_edges_str = "square adjacent edges";
   //  string square_sprt_edges_str = "square separate edges";
   //  string atom_edge_str = "edge";

   //  Symbolic_State<string> square(square_str,false);
   //  Symbolic_State<string> sqr_adjcnt_edges(square_adjcnt_edges_str,false);
   //  Symbolic_State<string> sqr_sprt_edges(square_sprt_edges_str,false);
   //  Symbolic_State<string> edge(atom_edge_str,true);

   //  //rules of square
   //  Symbolic_Rule<string> sqr_to_edges1(square,sqr_adjcnt_edges,sqr_adjcnt_edges);
   //  Symbolic_Rule<string> sqr_to_edges2(square,sqr_sprt_edges,sqr_sprt_edges);
   //  Symbolic_Rule<string> sqr_edges_to_edge1(sqr_adjcnt_edges,edge,edge);
   //  Symbolic_Rule<string> sqr_edges_to_edge2(sqr_sprt_edges,edge,edge);

   //  //construct T_AOG using rules defined above
   //  vector<Symbolic_Rule<string> > rules;
   //  rules.push_back(sqr_to_edges1);
   //  rules.push_back(sqr_to_edges2);
   //  rules.push_back(sqr_edges_to_edge1);
   //  rules.push_back(sqr_edges_to_edge2);
    
    // T_AOG<string> t_aog(rules);
    // cout<<"The number of vertices: "<<t_aog.NumberOfVertices()<<endl;
    // cout<<"The number of rules in this T_AOG: "<<t_aog.NumOfRules()<<endl
    // 	<<"The number of states in this T_AOG: "<<t_aog.NumOfStates()<<endl
    // 	<< "The number of leaf states in this T_AOG: " << t_aog.NumOfLeafStates()<<endl;

    // cout<<"T_AOG graph:"<<endl;
    // for(int i = 0; i< t_aog.NumberOfVertices();i++)
    //   {
   	// cout<<"parent node: " << i << endl<< "children: ";
   	// for(auto j : t_aog.ChildrenVertices(i))
   	//   cout<<j<<" ";
   	// cout<<endl;
    //   }
    // cout<<endl<<endl;
    // for(int i = 0; i<  t_aog.NumberOfVertices();i++)      
	// cout<<"the content of vertex "<< i<<" is: "<< t_aog.GetStateByVertexId(i).GetContent()<< endl;
      


   //  //sample a parse tree from the T_AOG
   //  cout<<"Sampling a parse tree from the T_AOG Graph..."<<endl;
   //  auto parse_tree = t_aog.Sample(0);

   //  for(auto iter : parse_tree)
   //  {
    
   //      cout<< "parent node: "<<iter.first<<endl;
   //      cout<<"children nodes: ";
   //      for(auto child : iter.second)
   //          cout<< child<<" ";
   //      cout<<endl;
   //  }

    return 0;
}
