#include <vector>
#include <string>

#include "../T-AOG/T_AOG.h"
#include "../Learner/Online_Learner.h"

using namespace std;
using namespace AOG_LIB;

int SAMPLING = 1000;

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
            outFile<<"\n";
        }
    }
    else
    {
        std::cerr<<"Error opening file!\n";
        throw exception();
    }
    outFile.close();
}

void GetOutput(shared_ptr<T_AOG<string> > graph)
{
    // Record results
	graph->Visualize();
    OutputToFile("sample_aog_KL_output.txt",*graph,graph->GetRoot());	
	vector<Symbolic_Rule<string> > rules = graph->GetRules();

    ofstream file;
	
	file.open("learned_tree.txt", ofstream::out|ofstream::trunc);
	if (file.is_open())
	{
		for (Symbolic_Rule<string> rule : rules)
		{
			// size of result, souce id, source content, [child id, child content]*
			file << rule.GetResults().size() << "," 
				 << rule.GetSource().GetId() << "," 
				 << rule.GetSource().GetContent();
			vector<Symbolic_State<string> > states = rule.GetResults();
			for (int i = 0; i < states.size(); i++)
			{
				file << "," << states[i].GetId() << "," << states[i].GetContent();
			}
			file << "\n";
		}
	}
	else
        cout << "Unable to open " << "learned_tree.txt" << endl;
}


int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cerr << "Wrong usage, ./cpp filename\n";
        exit(1);
    }
    string filename = string(argv[1]);
    ifstream f(filename, std::ifstream::in);
    
    vector<Symbolic_Rule<string> > rules;

    string s = "";
    int line_num = 0;
    while (getline(f, s))
    {
        line_num++;
        stringstream ss(s);
        string item;
        vector<string> tokens;
        while (getline(ss, item, ',')) 
            tokens.push_back(item);

        if (tokens.size() < 3)
        {
            cerr << "Error: row does not have source id and content\n";
            exit(1);
        }
        int num_res = stoi(tokens[0]);
        int src_id = stoi(tokens[1]);
        string src_ct = tokens[2];
        Symbolic_State<string> src(src_id);
        // if (tokens.size() != (1 + (1+num_res)*2))
        // {
        //     cerr << "Error: size doesn't match at line "<< line_num << "\n";
        //     cerr << "Token size is:" << tokens.size() << ", Expected size is:" << 1 + (1+num_res)*2 << "\n";
        //     exit(1);
        // }
        vector<Symbolic_State<string> > results;
        for (int i = 3; i < tokens.size(); i+=2)
        {
            Symbolic_State<string> state;
            int state_id = stoi(tokens[i]);
            if (state_id == -1)
            {
                string state_ct = tokens[i+1];            
                Symbolic_State<string> temp(state_ct, true);
                state = temp;
            }
            else
            {
                Symbolic_State<string> temp(state_id);
                state = temp;
            }
            results.push_back(state);
        }

        Symbolic_Rule<string> rule(src, results);
        rules.push_back(rule);
    }

    // Find root
    unordered_set<Symbolic_State<string> > top_level_rules;
    vector<Symbolic_State<string> > sources;
    unordered_set<Symbolic_State<string> > results;
    for (Symbolic_Rule<string> rule : rules)
    {
        sources.push_back(rule.GetSource());
        results.insert(rule.GetResults().begin(), rule.GetResults().end());
    }
    
    // check which source is not other sources' result
    for (Symbolic_State<string> source : sources)
    {
        if (results.find(source) == results.end())
            top_level_rules.insert(source);
    }
            
    if (top_level_rules.size() != 1)
    {
        cerr << "Dangling top level rules beside root\n";
        exit(1);
    }

    shared_ptr<T_AOG<string> > related_graph = make_shared<T_AOG<string> >();
    shared_ptr<AOG_Vertex<string> > source_ptr(new AOG_Vertex<string>(*top_level_rules.begin(), true, true));
	related_graph->AddVertex(source_ptr);

    for (auto rule : rules)
        related_graph->AddRule(rule);


    
    related_graph->OutputGraph("visualize_with_truncate");
}