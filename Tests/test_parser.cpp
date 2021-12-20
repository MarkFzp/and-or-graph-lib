//
// Created by Hao Wang on 18/3/31.
//

// #include <string>
// #include "./T-AOG/T_AOG.h"
#include "../Learner/Online_Learner.h"

using namespace std;
using namespace AOG_LIB;

int main()
{
    string S = "S";
    string NP = "NP";
    string VP = "VP";
    string PP = "PP";
    string Noun = "Noun";
    string Prep = "Prep";
    string John = "John";
    string Mary = "Mary";
    string Denver = "Denver";
    string Verb = "Verb";
    string called = "called";
    string from = "from";
    string NP_down = "NP_down";
    string VP_down = "VP_down";
    Symbolic_State<string> state_S(S, false);    
    Symbolic_State<string> state_NP(NP, false);
    Symbolic_State<string> state_NP_down(NP_down, false);    
    Symbolic_State<string> state_VP(VP, false);
    Symbolic_State<string> state_VP_down(VP_down, false);
    Symbolic_State<string> state_PP(PP, false);
    Symbolic_State<string> state_Noun(Noun, false);
    Symbolic_State<string> state_Prep(Prep, false);    
    Symbolic_State<string> state_John(John, false);
    Symbolic_State<string> state_Mary(Mary, true);
    Symbolic_State<string> state_Denver(Denver, true);
    Symbolic_State<string> state_Verb(Verb, false);
    Symbolic_State<string> state_called(called, true);
    Symbolic_State<string> state_from(from, true);
    vector<Symbolic_State<string> > NPVP = {state_NP, state_VP};
    vector<Symbolic_State<string> > NPPP = {state_NP_down, state_PP};
    vector<Symbolic_State<string> > Noun_ = {state_Noun};
    vector<Symbolic_State<string> > VerbNP = {state_Verb, state_NP};
    vector<Symbolic_State<string> > VPPP = {state_VP_down, state_PP};
    vector<Symbolic_State<string> > PrepNP = {state_Prep, state_NP};
    vector<Symbolic_State<string> > John_ = {state_John};
    vector<Symbolic_State<string> > Mary_ = {state_Mary};
    vector<Symbolic_State<string> > Denver_ = {state_Denver};
    vector<Symbolic_State<string> > called_ = {state_called};
    vector<Symbolic_State<string> > from_ = {state_from};
    
    Symbolic_Rule<string> S_to_NPVP(state_S, NPVP);
    Symbolic_Rule<string> NP_to_NPPP(state_NP, NPPP);
    Symbolic_Rule<string> NP_to_Noun(state_NP, Noun_);
    Symbolic_Rule<string> VP_to_VerbNP(state_VP, VerbNP);
    Symbolic_Rule<string> VP_to_VPPP(state_VP, VPPP);
    Symbolic_Rule<string> PP_to_PrepNP(state_PP, PrepNP);
    Symbolic_Rule<string> Noun_to_John(state_Noun, John_);
    Symbolic_Rule<string> Noun_to_Mary(state_Noun, Mary_);
    Symbolic_Rule<string> Noun_to_Denver(state_Noun, Denver_);
    Symbolic_Rule<string> Verb_to_called(state_Verb, called_);
    Symbolic_Rule<string> Prep_to_from(state_Prep, from_);
    
    vector<Symbolic_Rule<string> > rules;
    rules.push_back(S_to_NPVP);
    rules.push_back(NP_to_NPPP);
    rules.push_back(NP_to_Noun);    
    rules.push_back(VP_to_VerbNP);    
    rules.push_back(VP_to_VPPP);    
    rules.push_back(PP_to_PrepNP);    
    rules.push_back(Noun_to_John);    
    rules.push_back(Noun_to_Mary);    
    rules.push_back(Noun_to_Denver);    
    rules.push_back(Verb_to_called);    
    rules.push_back(Prep_to_from);
   

    shared_ptr<T_AOG<string> > related_graph = make_shared<T_AOG<string> >(rules);
    VertexId root_id;
    try
    {
        cerr << "inside main: " << rules[0].GetSource().GetContent() << endl; 
        
        root_id = related_graph->GetRoot();
    }
    catch(...)
    {
        std::cerr << "inside catch\n"; 
        Symbolic_State<string> source_state;
        std::shared_ptr<AOG_Vertex<string> > source_ptr(new AOG_Vertex<string>(source_state, false, true));
        related_graph->AddVertex(source_ptr);
    }
    AOG_LIB::Online_Learner<string> online_Learner(related_graph);

    vector<Symbolic_State<string> > sentence = {state_John, state_called, state_Mary, state_from, state_Denver};
    cerr << "----------------DEBUGGING STARTS HERE------------------" << endl;
    // online_Learner.GetBestParsing();
    online_Learner.MergeNewParsingWithGraph(sentence);
    
    
    
    return 0;
}