#include <vector>
#include <string>
#include <cassert>
#include <unordered_map>

#include "../T-AOG/T_AOG.h"
#include "../Learner/Online_Learner.h"

using namespace std;
using namespace AOG_LIB;


//output the sample result to a file
void OutputToFile(string filename, const T_AOG<string>& t_aog,VertexId sample_root)
{
    //open file
    ofstream outFile;
    outFile.open(filename);   
    if(outFile.is_open())
    {
        while (SAMPLING)
        {
			std::vector<VertexId> res;
			double prob;
			t_aog.Sample(sample_root, res, prob);
			SAMPLING--;
			for (VertexId id : res)
			{
				outFile << t_aog.GetStateByVertexId(id).GetContent() << " ";
			}
        }
    }
    else
    {
        std::cerr<<"Error opening file!\n";
        throw exception();
    }
    outFile.close();
}

int main(int argc, char *argv[])
{
    if(argc != 4)
    {
        cerr << "Wrong usage, ./cpp grammar_filename mapping_filename grammar_index\n";
        exit(1);
    }
    string filename = string(argv[1]);
    ifstream f(filename, ifstream::in);
    string mapping_file = string(argv[2]);
    ifstream mapping_f(mapping_file, ifstream::in);

    string s = "";
    // 114, c1a1
    unordered_map<string, string> terminal_content_mapping;
    cout << "Reading mappings...\n";
    while(getline(mapping_f, s))
    {
        stringstream ss(s);
        string item;
        vector<string> tokens;
        while (getline(ss, item, ',')) 
            tokens.push_back(item);

        if (tokens.size() != 2)
        {
            cerr << "Error: mapping file read error\n";
            exit(1);
        }

        terminal_content_mapping.insert(make_pair(tokens[1], tokens[0]));
        cout << "(" << tokens[1] << "," << tokens[0] << ")\n";
    }


    for (int i = 0; i < 4; i++)
        getline(f ,s);
    
    stringstream start_symbol(s);
    string root_content = "";
    getline(start_symbol, root_content, ' ');
    getline(start_symbol, root_content, ' ');    
    int root_id = stoi(root_content);
    cout << "Root id: " << root_id << "\n";

    while(getline(f, s))
    {
        if (s == "Terminals")
            break;
    }

    cout << "Reading terminals...\n";    
    // 1 114
    unordered_map<string, string> terminal_id_mapping;
    unordered_map<string, string> id_terminal_mapping;
    while (getline(f, s))
    {
        if (s == "AndNodes")
            break;

        stringstream ss(s);
        string item;
        vector<string> tokens;
        while (getline(ss, item, '\t')) 
            tokens.push_back(item);

        if (tokens.size() != 2)
        {
            cerr << "Error: reading terminals error\n";
            exit(1);
        }

        cout << "(" << tokens[0] << "," << tokens[1] << ")\n"; 
        terminal_id_mapping.insert(make_pair(tokens[0], tokens[1]));
        id_terminal_mapping.insert(make_pair(tokens[1], tokens[0]));
    }

    cout << "Reading AndNodes...\n";
    // 1 [1,2,3]
    unordered_map<string, vector<string> > andNode_children_mapping;
    while (getline(f, s))
    {
        if (s == "OrNodes")
            break;
        stringstream ss(s);
        string item;
        vector<string> tokens;
        while (getline(ss, item, '\t')) 
            tokens.push_back(item);

        if (tokens.size() != 2)
        {
            cerr << "Error: reading terminals error\n";
            exit(1);
        }       
        // get children
        vector<string> children;
        stringstream sss(tokens[1]);

        // [2 3 4 5 6 ]
        getline(sss, item, ']');
        // 2 3 4 5 6
        stringstream ssss(item.substr(1, item.length()-2));
        while (getline(ssss, item, ' '))
            children.push_back(item);
        
        andNode_children_mapping.insert(make_pair(tokens[0], children));
        cout << tokens[0] << "->\n\t";
        for (string child : children)
            cout << child << " ";
        cout << "\n";
    }

    cout << "Reading OrNodes...\n";
    unordered_map<string, vector<string> > orNode_children_mapping;
    unordered_map<string, vector<double> > orNode_weight_mapping;
    while (getline(f, s))
    {
        stringstream ss(s);
        string item;
        vector<string> tokens;
        while (getline(ss, item, '\t')) 
            tokens.push_back(item);

        if (tokens.size() != 2)
        {
            cerr << "Error: reading terminals error\n";
            exit(1);
        }     
        vector<string> children;
        stringstream sss(tokens[1]);

        // [2 3 4 5 6 ]
        getline(sss, item, ']');
        // 2 3 4 5 6
        stringstream ssss(item.substr(1, item.length()-2));
        while (getline(ssss, item, ' '))
            children.push_back(item);
        
        if (children.size() > 1)
        {
            orNode_children_mapping.insert(make_pair(tokens[0], children));
            cout << tokens[0] << "->\n\t";
            for (string child : children)
                cout << child << " ";
            cout << "\n";

            getline(sss, item, ']');
            vector<double> weights;
            stringstream ssss1(item.substr(2, item.length()-2));
            while (getline(ssss1, item, ' '))
                weights.push_back(stod(item));
            
            orNode_weight_mapping.insert(make_pair(tokens[0], weights));

            cout << tokens[0] << "->\n\t";
            for (string child : children)
                cout << child << " ";
            cout << "\n";
            cout << "\t";
            for (double weight: weights)
                cout << weight << " ";
            cout << "\n";
        }
        else
        {
            andNode_children_mapping.insert(make_pair(tokens[0], children));
            cout << "AndNode in Ornode:\n";
            cout << tokens[0] << "->\n\t";
            for (string child : children)
                cout << child << " ";
            cout << "\n";
        }
    }

    shared_ptr<T_AOG<string>> graph_ptr(new T_AOG<string>());
    Symbolic_State<string> root(root_id);
    std::shared_ptr<AOG_Vertex<string> > source_ptr(new AOG_Vertex<string>(root, true, true));
	graph_ptr->AddVertex(source_ptr);

    for (auto it = andNode_children_mapping.begin(); it != andNode_children_mapping.end(); it++)
    {
        Symbolic_State<string> source(stoi(it->first));
        vector<Symbolic_State<string> > results;
        for (string child : it->second)
        {
            // if the child is a terminal
            if (terminal_id_mapping.find(child) != terminal_id_mapping.end())
                results.push_back(Symbolic_State<string>(terminal_content_mapping[terminal_id_mapping[child]], true));
            else
                results.push_back(Symbolic_State<string>(stoi(child)));
        }
        graph_ptr->AddRule(Symbolic_Rule<string>(source, results));
    }

    for (auto it = orNode_children_mapping.begin(); it != orNode_children_mapping.end(); it++)
    {
        Symbolic_State<string> source(stoi(it->first));

        for (string child : it->second)
        {
            vector<Symbolic_State<string> > results;
            if (terminal_id_mapping.find(child) != terminal_id_mapping.end())
                results.push_back(Symbolic_State<string>(terminal_content_mapping[terminal_id_mapping.at(child)], true));
            else
                results.push_back(Symbolic_State<string>(stoi(child)));
            graph_ptr->AddRule(Symbolic_Rule<string>(source, results));
        }

        VertexId source_id = graph_ptr->GetVertexIdByState(source);
        vector<VertexId> dummy_id = graph_ptr->ChildrenVertices(source_id);
        // assert(dummy_id.size() == it->second.size());
        unordered_map<VertexId, double> weight_map;  
        for (VertexId id : dummy_id)
        {
            vector<VertexId> children_id = graph_ptr->ChildrenVertices(id);
            assert(children_id.size() == 1);
            double weight;
            const auto& children = orNode_children_mapping[it->first];
            int state_id = graph_ptr->GetStateByVertexId(children_id[0]).GetId();
            if (state_id == -1){
                string state_content = graph_ptr->GetStateByVertexId(children_id[0]).GetContent();
                for (int i = 0; i < children.size(); i++)
                {
                    if(terminal_id_mapping.find(children[i]) == terminal_id_mapping.end()){
                        continue;
                    }
                    else if(terminal_content_mapping[terminal_id_mapping.at(children[i])] == state_content){
                        weight = orNode_weight_mapping[it->first][i];
                        break;
                    }
                }
            }
            else{
                for (int i = 0; i < children.size(); i++)
                {
                    if(terminal_id_mapping.find(children[i]) != terminal_id_mapping.end()){
                        continue;
                    }
                    else if (stoi(children[i]) == state_id){
                        weight = orNode_weight_mapping[it->first][i];
                        break;
                    }
                }
            }
            weight_map.insert(make_pair(id, weight));
        }
        graph_ptr->SetOutEdgeWeights(source_id, weight_map);
    }



    //output learned_tree.txt for benchmark
    auto rules = graph_ptr -> GetRules();

    ofstream file;
    string learned_tree_path = (string)argv[3];

	file.open( learned_tree_path, ofstream::out | ofstream::trunc);
	if (file.is_open())
	{
        std::cerr << "File opened. \n";
        for (const Symbolic_Rule<std::string>& rule : rules)
        {
            // size of result, souce id, source content, [child id, child content]*
            file << rule.GetResults().size();
            std::vector<Symbolic_State<std::string> > states = rule.GetResults();
            VertexId source_vtx_id = graph_ptr->GetVertexIdByState(rule.GetSource());
            bool isAnd = graph_ptr->GetVertexContent(source_vtx_id)->IsAnd();
            // std::cerr << "state_id: " << rule.GetSource().GetId() << ", vertex_id: " << source_vtx_id << ", isAnd: " << isAnd << std::endl;
            if(!isAnd){
                std::vector<VertexId> children_vector;
                for(const Symbolic_State<std::string>& state : states){
                    children_vector.push_back(graph_ptr->GetVertexIdByState(state));
                }
                VertexId tmpDummy;
                bool found = false;
                for (VertexId dummy : graph_ptr->ChildrenVertices(source_vtx_id))
                {
                    if (graph_ptr->ChildrenVertices(dummy) == children_vector)
                    {
                        found = true;
                        tmpDummy = dummy;
                        break;
                    }
                }
                if (!found){
                    std::cerr << "rule not found" << std::endl;
                    throw std::exception();
                }
                double weight = graph_ptr->GetOutEdgeWeights(source_vtx_id, false)[tmpDummy];
                if (weight != 0)
                    file << "," << weight;
                else{
                    std::cerr << "0 weight is found !!" << std::endl;
                    throw std::exception();
                }
            }
            else{
                file << "," << 0;
            }
            file << "," << rule.GetSource().GetId() << "," << rule.GetSource().GetContent();
            for (int i = 0; i < states.size(); i++)
            {
                file << "," << states[i].GetId() << "," << states[i].GetContent();
            }
            file << "\n";
        }
	}
	else
		cout << "Unable to open "
			 << "learned_tree.txt" << endl;

	file.close();

    // benchmark visualize file
    // graph_ptr->OutputGraph("visualize_benchmark_"+(string)argv[3]);

    // sample benchmark for KL_divergence
    // OutputToFile("benchmark_sample.txt", *graph_ptr, graph_ptr->GetRoot());
}
