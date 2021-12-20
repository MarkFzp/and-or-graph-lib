#ifndef METRICS_H
#define METRICS_H

#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<random>
#include<algorithm>
#include<cstdlib>
#include<string>
#include<vector>
#include<memory>
#include<cmath>
#include<utility>
#include<functional>
#include<unordered_map>
#include<unordered_set>
#include <boost/filesystem.hpp>
#include "../Utils/Earley_Parser_Evaluation.h"
#include "../Learner/Online_Learner.h"
#include "../T-AOG/T_AOG.h"
#include <boost/functional/hash.hpp>
#include <ctime>

using namespace std;
using namespace AOG_LIB;

#ifndef CONTAINER_HASH
#define CONTAINER_HASH
template <typename Container> // (hasher for containers) we make this generic for any container
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};
#endif


//Reconstruct the learned T-AOG from learned_tree.txt
shared_ptr<T_AOG<string> > ReconstructTree(string filename);

shared_ptr<T_AOG<string> > ReconstructTree_UniformOrNode(string filename);

//Generate unique sampled data from all sampled data
void Generate_unique_data(string file_path, string file_name);

//Calculate Precision given the learned T-AOG and the test set
double Precision( string learned_tree_path, string ground_truth_path, int sample_time = 1000);

//given outputs generated from learned grammar
double Precision2( string sample_result_path, string ground_truth_path );

//Calculate Recall given the learned T-AOG and the test set
double Recall( string learned_tree_path, string test_set_path);

void Sample_Viterbi( string learned_tree_path, unsigned sample_time, bool unique);

//Parse the input sequence with the given parser
// int Parse(const SequenceType<string> &input_seq, std::vector<Symbolic_Rule<string> > all_rules, shared_ptr<EarleyParser<string> > parser);

void Train_test_split(string dataset_path, string filename, double train_ratio=0.7, bool shuffle = true);

void Train_test_exclusive(string train_set_path, string test_set_path);      

double Perplexity_cross_entropy(string learned_tree_path, string test_set_path, double penalty);

double KL_Divergence(string learned_tree_path, string train_set_path);

double Bleu (string reference_file, string hypothesis_file);
void ConvertToSymbols(vector<Symbolic_Rule<string> > rules, ofstream& out);

void GetBestParsingTree(string learned_tree_path, string test_set_path);
/* 
* This function generates a dataset contains only unique entries from the given dataset
* Params:
*   filesc_path: the path to the directory that contains the dataset
*   file_name: the dataset name 
*/
void Generate_unique_data(string file_path, string file_name)
{
    
    vector<SequenceType<string>> dataset = FileParser<string>(file_path + file_name);
    unordered_set<SequenceType<string>> unique_data;
    for(const auto & data : dataset)
        unique_data.insert(data);
    
    ofstream outfile;
    outfile.open(file_path + "../../Unique_Grammar_Example/Output/" + file_name, ofstream :: trunc);
    if(!outfile.is_open())
    {
        std::cerr<<"cannot open file " +file_name + " !\n";
        throw exception();
    }

    for(const auto & data : unique_data)
    {
        for(const auto & state : data)
        {
            outfile << state.GetContent() << " ";
        }
        outfile << "\n";
    }

    outfile.close();
   

}


void Train_test_exclusive(string train_set_path, string test_set_path)
{
    vector<SequenceType<string> > train_set = FileParser<string>(train_set_path);
    vector<SequenceType<string> > test_set = FileParser<string>(test_set_path);

    vector<SequenceType<string> > exclusive_train_set;
    for(const auto & train_data : train_set)
    {
        if(find(test_set.begin(),test_set.end(),train_data) == test_set.end())
        {
            exclusive_train_set.push_back(train_data);
        }
    }

    ofstream train_file;
    boost::filesystem::path p(train_set_path);
    train_file.open(p.parent_path().string()+"/" + p.stem().string() + "_exclusive.txt");
    if(!train_file.is_open())
    {
        cerr <<"cannot open file!\n";
        throw exception();
    }
    printf("haha");
    fflush(stdout);
    for(const auto & data : exclusive_train_set)
    {
        for(const auto & state : data)
        {
            train_file << state.GetContent()<<" ";
        }
        train_file << endl;
    }
    train_file.close();
}                  



/* 
* This function split a dataset into training and testing set with the given train/test ratio
* Params:
*   file_path: the path to the directory that contains the dataset you want to split
*   file_name: the dataset name 
*   train_ratio: the training set ratio with respect to the whole dataset. Default is 0.7 
*   shuffle: shuffle the dataset before splitting. Default is true
*/
void Train_test_split(string dataset_path, string filename, double train_ratio,bool shuffle )
{
    string folder;
    if(filename == "grammar_example5_output.txt")
        folder = "grammar5/";
    else if(filename == "grammar_example6_output.txt")
        folder = "grammar6/";
    else if(filename == "grammar_example7_output.txt")
        folder = "grammar7/";
    else if(filename == "grammar_example8_output.txt")
        folder = "grammar8/";
    else
        folder = filename.substr(0,filename.find('.',0)) + '/';
    
    string split_path = dataset_path + "../Train_test_split/";
    const char* path = split_path.c_str();
    boost::filesystem::path dir(path);
    if(boost::filesystem::create_directory(dir))
    {
        std::cerr<< "Directory Created: "<<split_path<<std::endl;
    
    }
    split_path =  split_path + folder;
    path = split_path.c_str();
    boost::filesystem::path new_dir(path);
    if(boost::filesystem::create_directory(new_dir))
    {
        std::cerr<< "Directory Created: "<<split_path<<std::endl;
    
    }
    vector<SequenceType<string> > dataset = FileParser<string>(dataset_path + filename);
	
    int dataset_size = dataset.size();
    int train_size = floor(dataset_size * train_ratio);
    int test_size = dataset_size - train_size;
    
    //shuffle the dataset
    if(shuffle)
        random_shuffle(dataset.begin(),dataset.end());
    vector<SequenceType<string> > train_set(dataset.begin(),dataset.begin() + train_size);
    vector<SequenceType<string> > test_set(dataset.begin()+train_size,dataset.end());

    //output to file
    ofstream trainfile;
    ofstream testfile;
    trainfile.open(split_path + "train.txt", ofstream :: trunc);
    testfile.open(split_path + "test.txt",ofstream :: trunc);
    if(!trainfile.is_open())
    {
        std::cerr<<"cannot open file !\n";
        throw exception();
    }
    if(!testfile.is_open())
    {
        std::cerr<<"cannot open file !\n";
        throw exception();
    }

    for(const auto & data : train_set)
    {
        for(const auto & state : data)
        {
            trainfile << state.GetContent() << " ";
        }
        trainfile << "\n";
    }

    for(const auto & data : test_set)
    {
        for(const auto & state : data)
        {
            testfile << state.GetContent() << " ";
        }
        testfile << "\n";
    }

    trainfile.close();
    testfile.close();
}

// int Parse(const SequenceType<string> &input_seq,
//                         std::vector<Symbolic_Rule<string> > all_rules, shared_ptr<EarleyParser<string> > parser)
// {

//     //make a copy of the passed in parser
//     shared_ptr<EarleyParser<string> > new_parser = make_shared<EarleyParser<string> >(*parser);

//     // if all rules are eliminated, no rule can parse the input sequence
//     // move to the next input position since there is no parsing for starting at this position
//     if (all_rules.size() == 0)
//     {
//         /** std::cerr << "EarleyParserWrapper:NO RULE TO PARSE THIS INPUT" << std::endl; */
//         return -1;
//     }
//     bool temp = true;
//     int pos = new_parser -> parse(input_seq.begin(), input_seq.end(),std::cout,temp);
   

//     // earley parser parse with rules & states
//     // if return 0, parse with reduced rules, else return position
//     if (pos == 0)
//     {
//         // if cannot parse from start, erase all top level rules
//         SequenceType<string> top_level_rules = new_parser -> get_top_level_rules(all_rules);
//         // if cannot further erase top level rules, no rule can parse the input sequence
//         if (top_level_rules.empty())
//         {
//             std::cerr << "NO WAY TO FURTHER TRUNCATE THE RULES" << std::endl; 
//             return -1;
//         }

//         for (typename std::vector<Symbolic_Rule<string> >::iterator it = all_rules.begin();
//                 it < all_rules.end(); it++)
//         {
//             if (std::find(top_level_rules.begin(), top_level_rules.end(), it->GetSource()) != top_level_rules.end())
//             {        
//                 it = all_rules.erase(it);
//                 it--;               
//             }
//         }
       
//         new_parser ->update(all_rules,input_seq);
//         return Parse(input_seq, all_rules, new_parser);
//     }

//     // else, return the position that has been already parsed
//     return pos;
// }

// double Prob(const SequenceType<string> &input_seq,
//                         std::vector<Symbolic_Rule<string> > all_rules, shared_ptr<EarleyParser<string> > parser, double &prob, std::shared_ptr<T_AOG<std::string> > graph_ptr)
// {

//     //make a copy of the passed in parser
//     shared_ptr<EarleyParser<string> > new_parser = make_shared<EarleyParser<string> >(*parser);

//     // if all rules are eliminated, no rule can parse the input sequence
//     // move to the next input position since there is no parsing for starting at this position
//     if (all_rules.size() == 0)
//     {
//         /** std::cerr << "EarleyParserWrapper:NO RULE TO PARSE THIS INPUT" << std::endl; */
//         return -1;
//     }
//     bool temp = true;
//     int pos = new_parser -> parse2(input_seq.begin(), input_seq.end(),std::cout, temp, prob, graph_ptr);
   

//     // earley parser parse with rules & states
//     // if return 0, parse with reduced rules, else return position
//     if (pos == 0)
//     {
//         // if cannot parse from start, erase all top level rules
//         SequenceType<string> top_level_rules = new_parser -> get_top_level_rules(all_rules);
//         // if cannot further erase top level rules, no rule can parse the input sequence
//         if (top_level_rules.empty())
//         {
//             std::cerr << "NO WAY TO FURTHER TRUNCATE THE RULES" << std::endl; 
//             return -1;
//         }

//         for (typename std::vector<Symbolic_Rule<string> >::iterator it = all_rules.begin();
//                 it < all_rules.end(); it++)
//         {
//             if (std::find(top_level_rules.begin(), top_level_rules.end(), it->GetSource()) != top_level_rules.end())
//             {        
//                 it = all_rules.erase(it);
//                 it--;               
//             }
//         }
       
//         new_parser ->update(all_rules,input_seq);
//         return Prob(input_seq, all_rules, new_parser, prob, graph_ptr);
//     }

//     // else, return the position that has been already parsed
//     return pos;
// }

/*
This function calculates the precision of a learned grammar given ground truth grammar
Param:
    learned_tree_path: the path to learned grammar learned_tree.txt
    ground_truth_path: the path to the ground truth grammar learned_tree.txt
    sample time: the number of time to sample from the learned grammar, the sampled
    data is used to calculate precision of the learned grammar. Default is set to 1000.

*/
double Precision( string learned_tree_path, string ground_truth_path, int sample_time )
{
    shared_ptr<T_AOG<string> > learned_graph_ptr = ReconstructTree(learned_tree_path);
    shared_ptr<T_AOG<string> > gt_graph_ptr = ReconstructTree(ground_truth_path);
    
    vector<Symbolic_Rule<string> > all_rules = gt_graph_ptr->GetRules();
    Symbolic_State<string> root_state = gt_graph_ptr -> GetStateByVertexId(gt_graph_ptr ->GetRoot());
    shared_ptr<grammar<string> > g = make_shared<grammar<string> >(all_rules, vector<Symbolic_State<string> >{root_state});                                        
    
    unordered_map<vector<VertexId>, pair<unordered_set<size_t>, double>, container_hash<vector<VertexId> > > samples_map;
    container_hash<vector<VertexId> > hash_v;

    while(sample_time)
    {
        --sample_time;
        vector<VertexId > seq;
        double prob;
        shared_ptr<vector<VertexId>> parse_tree = learned_graph_ptr->Sample(learned_graph_ptr->GetRoot(), seq, prob);
        size_t parse_tree_hash = hash_v(seq);

        if(samples_map.find(seq) == samples_map.end()){
            unordered_set<size_t> parse_trees_hash ({parse_tree_hash});
            samples_map[seq] = make_pair(parse_trees_hash, prob);
        }
        else{
            auto &value = samples_map[seq];
            unordered_set<size_t> &parse_trees_hash = value.first;
            if(parse_trees_hash.find(parse_tree_hash) == parse_trees_hash.end()){
                parse_trees_hash.insert(parse_tree_hash);
                value.second += prob;
            }
        }
    }

    double weighted_parsed = 0;
    double total_prob = 0;
    unsigned unique_data_count = 0;
    for(const auto &it : samples_map){
        unique_data_count++;
        double prob = it.second.second;
        total_prob += prob;

        vector<Symbolic_State<string> > seq;
        for(const VertexId &iit : it.first){
            seq.push_back(learned_graph_ptr->GetStateByVertexId(iit));
        }
        shared_ptr<EarleyParser<string> > parser = make_shared<EarleyParser<string> >(*g);
        bool parsing_success = false;
        int pos = parser->parse(seq.begin(), seq.end(), std::cout, parsing_success);
        if(pos == seq.size()){
            weighted_parsed += prob;
        }
    }

    cout << "unique data count: " << unique_data_count << endl;

    return weighted_parsed / total_prob;
}

double Precision2( string sample_result_path, string ground_truth_path )
{
    shared_ptr<T_AOG<string> > gt_graph_ptr = ReconstructTree(ground_truth_path);
    vector<Symbolic_Rule<string> > all_rules = gt_graph_ptr->GetRules();
    Symbolic_State<string> root_state = gt_graph_ptr -> GetStateByVertexId(gt_graph_ptr ->GetRoot());
    shared_ptr<grammar<string> > g = make_shared<grammar<string> >(all_rules, vector<Symbolic_State<string> >{root_state});                                        
    shared_ptr<EarleyParser<string> > parser = make_shared<EarleyParser<string> >(*g);

    ifstream f(sample_result_path, std::ifstream::in);

	std::cout << "reading " << sample_result_path << "\n";
    
    unordered_set<vector<Symbolic_State<string>>> unique_sampled_dataset;

    string s = "";

    while (getline(f, s))
    {    
        istringstream iss(s);
        vector<Symbolic_State<string>> sample;
        do
        {
            string subs;
            iss >> subs;
            sample.push_back(Symbolic_State<string>(subs, true));
            
        } while (iss);
        unique_sampled_dataset.insert(sample);
    }


    unsigned total = unique_sampled_dataset.size();
    unsigned parsed = 0;
    for(const auto & data : unique_sampled_dataset)
    {
        shared_ptr<EarleyParser<string> > parser = make_shared<EarleyParser<string> >(*g);
        bool parsing_success = false;
        int pos = parser->parse(data.begin(), data.end(), std::cout, parsing_success);
        // if((!parsing_success && pos == data.size()) || (parsing_success && pos != data.size()))
        // {
        //     std::cerr<<"the parser assumption goes wrong in metrics!\n";
        //     throw std::exception();
        // }
        if(pos == data.size())
            ++parsed;
        // for(const auto & state : data)
        // {
        //     cout <<"("<<state.GetContent()<<", "<<state.GetId()<<")";
        // }
        // cout <<endl;
        // cout <<"pos: "<<pos<<"\n";

        
    }

    cout <<"parsed: "<< parsed << "\n";
    cout << "total unique data size: " << total << endl;
    return ((double) parsed) / total;
}

void Sample_Viterbi ( string learned_tree_path, unsigned sample_time, bool unique )
{
    shared_ptr<T_AOG<string> > learned_graph_ptr = ReconstructTree(learned_tree_path);
    vector<Symbolic_Rule<string> > all_rules = learned_graph_ptr->GetRules();
    Symbolic_State<string> root_state = learned_graph_ptr -> GetStateByVertexId(learned_graph_ptr->GetRoot());
    shared_ptr<grammar<string> > g = make_shared<grammar<string> >(all_rules, vector<Symbolic_State<string> >{root_state});
    
    //sample from the learned AOG and calculate precision
    unordered_map<vector<VertexId>, pair<unordered_set<size_t>, double>, container_hash<vector<VertexId> > > samples_map;
    // set<vector<VertexId> > sampled_dataset;

    container_hash<vector<VertexId> > hash_v;

    if(unique){
        while(sample_time)
        {
            --sample_time;
                vector<VertexId> seq;
            double prob;
            shared_ptr<vector<VertexId>> parse_tree = learned_graph_ptr->Sample(learned_graph_ptr->GetRoot(), seq, prob);
            size_t parse_tree_hash = hash_v(seq);

            if(samples_map.find(seq) == samples_map.end()){

                unordered_set<size_t> parse_trees_hash ({parse_tree_hash});
                samples_map[seq] = make_pair(parse_trees_hash, prob);
                
                //bool parsing_success = false;
                // int pos = parser->parse2(seq.begin(), seq.end(), std::cout, parsing_success, prob, learned_graph_ptr);
                // if (pos != seq.size()){
                //     std::cerr << "sampled data not parsable" << std::endl;
                //     throw exception();
                // }
                // double prob = parser->prob(learned_graph_ptr);
                
            }
            else{
                auto &value = samples_map[seq];
                unordered_set<size_t> &parse_trees_hash = value.first;
                if(parse_trees_hash.find(parse_tree_hash) == parse_trees_hash.end()){
                    parse_trees_hash.insert(parse_tree_hash);
                    value.second += prob;
                }
            }
        }

        double total_prob = 0;
        for(const auto &it : samples_map){
            double prob = it.second.second;
            total_prob += prob;
            cout << prob << ": ";
            for(const VertexId &iit : it.first){
                cout << learned_graph_ptr->GetStateByVertexId(iit).GetContent() << ' ';
            }
            cout << endl;
        }


        cerr << "total prob: " << total_prob << endl;
    }
    else{
        while(sample_time){
            sample_time--;
            vector<VertexId> seq;
            double prob;
            learned_graph_ptr->Sample(learned_graph_ptr->GetRoot(), seq, prob);
            for(const VertexId& id : seq){
                cout << learned_graph_ptr->GetStateByVertexId(id).GetContent() << ' ';
            }
            cout << endl;
        }
    }

}



// void Sample_UniformOrNode(string learned_tree_path, unsigned sample_time){
//     shared_ptr<T_AOG<string> > learned_graph_ptr = ReconstructTree_UniformOrNode(learned_tree_path);
//     vector<Symbolic_Rule<string> > all_rules = learned_graph_ptr->GetRules();
//     Symbolic_State<string> root_state = learned_graph_ptr -> GetStateByVertexId(learned_graph_ptr->GetRoot());
//     shared_ptr<grammar<string> > g = make_shared<grammar<string> >(all_rules, vector<Symbolic_State<string> >{root_state});
    
//     //sample from the learned AOG and calculate precision
//     unordered_set<vector<VertexId>, container_hash<vector<VertexId>> > sampled_dataset;

//     double parsed = 0;
//     double total_sample_time = sample_time;
//     while(sample_time)
//     {
//         vector<Symbolic_State<string> > seq;
//         learned_graph_ptr->Sample(learned_graph_ptr->GetRoot(),seq);
//         vector<VertexId> seq_id;
//         for(const auto& state : seq){
//             seq_id.push_back(learned_graph_ptr->GetVertexIdByState(state));
//         }

//         if(sampled_dataset.find(seq_id) == sampled_dataset.end()){
//             sampled_dataset.insert(seq_id);
//             shared_ptr<EarleyParser<string> > parser = make_shared<EarleyParser<string> >(*g);
//             bool parsing_success = false;
//             int pos = parser->parse(seq.begin(), seq.end(), std::cout, parsing_success);
//             if (pos != seq.size()){
//                 std::cerr << "sampled data not parsable" << std::endl;
//                 throw exception();
//             }
//             double prob = parser->prob(learned_graph_ptr);
//             cout << prob << ": ";
//             for(const auto & state : seq)
//                 cout << state.GetContent() << " " ;
//             cout << endl;
//             --sample_time;
//         }
//     }

// }
Symbolic_State<string> recur_helper(const unordered_map<Symbolic_State<string>, vector<Symbolic_State<string> > >& src_to_tar, Symbolic_State<string> cur_state, bool left,ofstream& out)
{
    if(cur_state.GetIsBasic())
        return cur_state;
    
    auto tar_states = src_to_tar.find(cur_state)->second;
    if(src_to_tar.find(cur_state) == src_to_tar.end())
    {
        cerr << "cur_state not in map!\n";
        throw std::exception();
    }
    Symbolic_State<string> left_child = recur_helper(src_to_tar, tar_states[0], true,out);

    Symbolic_State<string> right_child = recur_helper(src_to_tar, tar_states.back(), false,out);

    //output the current node's set
    out << left_child.GetContent()<<"," << right_child.GetContent()<<",";

    if(left)
        return left_child;
    return right_child;

} 

void ConvertToSymbols(vector<Symbolic_Rule<string> > rules, ofstream& out)
{
    unordered_map<Symbolic_State<string>, vector<Symbolic_State<string> > > src_to_tar; 
    for(auto& rule : rules)
        src_to_tar[rule.GetSource()] = rule.GetResults();

    if(rules.size() == 0)
      {
	cerr << "the rule is empty!\n"<<endl;
	throw exception();
      }
    auto root_state = rules[0].GetSource();
    auto children = src_to_tar.find(root_state)->second;
    auto left_child = recur_helper(src_to_tar, children[0], true, out);
    auto right_child = recur_helper(src_to_tar, children.back(),false,out);
    out << left_child.GetContent()<<"," << right_child.GetContent()<<","<<endl;
}


void GetBestParsingTree(string learned_tree_path, string test_set_path)
{
    ofstream out;
    out.open("learned_parse_tree.txt",ofstream::out | ofstream:: trunc);

     //read in test set
    vector<SequenceType<string>> test_set = FileParser<string>(test_set_path);

    //create parser
    shared_ptr<T_AOG<string> > graph_ptr = ReconstructTree(learned_tree_path);
    vector<Symbolic_Rule<string> > all_rules = graph_ptr->GetRules();
    Symbolic_State<string> root_state = graph_ptr -> GetStateByVertexId(graph_ptr ->GetRoot());
    shared_ptr<grammar<string> > g = make_shared<grammar<string> >(all_rules,
                                                            vector<Symbolic_State<string> >{root_state});

    for(int i = 0; i < test_set.size();i++)
    {
        for(int j = 0; j < test_set[i].size(); j++)
        {
            string word = test_set[i][j].GetContent();
            if(j != test_set[i].size()-1)
                word += " ";
            out << word;

        }
        out << endl;

        shared_ptr<EarleyParser<string> > parser = make_shared<EarleyParser<string> >(*g);
        bool parsing_success = false;
        int pos = parser -> parse(test_set[i].begin(), test_set[i].end(),std::cout,parsing_success, true);
       
        if(pos == test_set[i].size())
        {
	  /* cerr << "pos: "<<pos <<" sent size: "<<test_set[i].size()<<endl; */
            vector<Symbolic_Rule<string> > best_rule_in_one_parse;

            std::vector<std::vector<typename StateList<string>::statelist > > parsed_results = parser->GetPartialParse();
	    /* cerr << "size of parsed results: "<<parsed_results.size()<<endl; */
            for(int idx = 0; idx < parsed_results.size(); ++idx)
            {
                //error check
                if(idx >0)
                {
                    std::cerr<<"Parsing assumption failed!\n";
                    throw std::exception();
                }

                // parse_tree_prob_chart.push_back(std::vector<double>());
                int last_statelist_index = parsed_results[idx].size() - 1;
                double sub_parse_tree_prob = 0;
                std::unordered_map<Symbolic_Rule<string>, double> temp_rule_to_prob;
                double best_parse_tree_prob = INT_MIN;

                for(state<string> st : parsed_results[idx][last_statelist_index])
                {
                    // find top most root
                    if (st.i_ == 0 && st.j_ != 0 && st.rule_->left().GetContent() == "$" && st.completed())
                    {
		      /* cerr<<"find dollar sign!!!!!\n"; */
                        /** std::cerr << "The state is: [" << st.i_ << "," << st.j_ << "], FROM " << st.rule_->left().GetContent() << " TO ";
                        auto results =  st.rule_->right()[0];
                        for (auto result : results)
                        {
                            std::cerr << result.GetContent() << " ";
                        }
                        std::cerr << std::endl; */
                        // @param prob: probability of this particular possible parsing
                        double parse_tree_prob = 1;
                        std::vector<Symbolic_Rule<string> > rules_in_one_parse;
                        // @param q: queue used to backtrack all subparsing in one possible parsing
                        std::queue<state<string> > q;
                        // @param possible_parsing: possible parsing in one chart
                        Symbolic_State<string> possible_parsing_start;
                        q.push(st);
                        while (!q.empty())
                        {
                            st = q.front();
                            q.pop();
                            // this rule must be fully parsed to be considered as a valid parse
                            /** std::cerr << "back pointer size is: " << st.back_pointer_.size() << std::endl; */
                            for (std::pair<int, int> back_pointer : st.back_pointer_)
                            {
                                auto pointed_state = parsed_results[idx][back_pointer.first][back_pointer.second];
                                
                                q.push(parsed_results[idx][back_pointer.first][back_pointer.second]);
                            }

                            // skip the dummy rule from $ -> rest
                            if (st.rule_->left().GetContent() == "$")
                            {
                                assert(st.rule_->right()[st.right_].size() == 1);
                                possible_parsing_start = st.rule_->right()[st.right_][0];
                                continue;
                            }
                            /** std::cerr << "parent state is:" << st.rule_->left().GetContent() << std::endl; */
                            VertexId parentId = graph_ptr->GetVertexIdByState(st.rule_->left());
                            // try
                            // {
                            //     // auto all_states = this->graph_.GetStates();
                            //     /** std::cerr<<"size of states in this->graph_: "<<all_states.size()<<std::endl; */
                            //     parentId = this->graph_.GetVertexIdByState(st.rule_->left());
                                
                            // }
                            // catch(std::exception e)
                            // {
                                
                            //     // for(auto state : all_states)
                            //     //     std::cerr<<state.GetId()<<std::endl;
                            //     std::cerr<<"exception caught in EM!\n";
                            //     throw std::exception();
                            //     // AOFStruct<string> peeker;
                            //     /** std::cerr << "AOF id_: " << peeker.id_; */
                            //     // auto all_states = this->graph_.GetStates();
                            //     /** std::cerr<<"size of states in this->graph_: "<<all_states.size()<<std::endl;
                            //     for(auto state : all_states)
                            //     {
                            //         std::cerr<<"State Content: "<<state.GetContent()<<"_\n"<<"State ID: "<<state.GetId()<<std::endl;
                            //     }

                            //     std::cout << e.what() << std::endl; */
                            // }
                            // if it is Or-node, find the corresponding weight and update likelihood of this possible parsing
                            
                            if (!graph_ptr->GetVertexContent(parentId)->IsAnd())
                            {
                                /** std::cerr << "The rule's source is an Or-node" << std::endl; */
                                // find all out edge weights of the source
                                std::unordered_map<VertexId, double> outEdgeWeights =
                                    graph_ptr->GetOutEdgeWeights(parentId, true);

                                // locate the outedge we are looking for
                                // find all right-hand-side vertexids
                                // @param right_hand_state_ids: the ids of right hand side states (e.g. S->NP VP. ids of NP, VP)
                                std::vector<VertexId> right_hand_state_ids;
                                for (Symbolic_State<string> right_hand_state : st.rule_->right()[st.right_])
                                    right_hand_state_ids.push_back(
                                            graph_ptr->GetVertexIdByState(right_hand_state));
                                // find the children that has all the right-hand-side vertex,
                                // i.e. corresponds to the right-hand-rule
                                // @param dummy_vertices: the ids of the dummy vertices under the Or-node of parent (e.g. S->NP VP. dummy nodes under S)
                                std::vector<VertexId> dummy_vertices =
                                        graph_ptr->ChildrenVertices(parentId);

                                for (VertexId id : dummy_vertices)
                                {
				  /* cerr<<"inside dummy id!\n"; */
                                    // @param found_destination: flag indicating the (e.g. S->NP VP) branch is found
                                    // @param children_vertices: children vertices of the dummy node we are looking at
                                    bool found_destination = true;
                                    std::vector<VertexId> children_vertices =
                                            graph_ptr->ChildrenVertices(id);

                                    // if the children of the dummy node is of different size than our target, then this is not the branch we are looking for
                                    if (children_vertices.size() != right_hand_state_ids.size())
                                        continue ;
                                    // if size is equal, check all children are same
                                    for (int i = 0; i < children_vertices.size(); i++)
                                    {
                                        if (children_vertices[i] != right_hand_state_ids[i])
                                        {
                                            found_destination = false;
                                            break;
                                        }
                                    }
                                    if (found_destination)
                                    {
                                        //record this rule
                                        Symbolic_Rule<string> cur_rule = {st.rule_->left(), st.rule_->right()[st.right_]};                                
                                        rules_in_one_parse.push_back(cur_rule);
					/* cerr<<"rule in one parse pushed in or_node\n"; */
                                        // rules_under_or_node[st.rule_->left()].insert(cur_rule);

                                        //record the mapping from this dummy node to the rule
                                        // dummy_to_rule.insert({id, cur_rule});
                                        double weight = outEdgeWeights[id];
                                        parse_tree_prob *= weight;
                                        // std::cout << "weight is: " << weight << ", probability is: " << prob << std::endl;
                                        break;
                                    }
                                    
                                }
                            
                            }
                            else
			      {
				rules_in_one_parse.push_back({st.rule_->left(), st.rule_->right()[st.right_]});
				/* cerr <<"rule in one parse pushed in and-node\n"; */
			      }


                            
                        }
                        if(parse_tree_prob > best_parse_tree_prob)
                        {
			  /* cerr<<"parse_tree_prob greater!!\n"; */
                            best_parse_tree_prob = parse_tree_prob;
                            best_rule_in_one_parse = rules_in_one_parse;
                        }

                    }
                    


                }


                // //error checking
                // if(sub_parse_tree_prob == 0)
                // {
                //     std::cerr<<"the data is not partial parsable!\n";
                //     throw std::exception();
                // }
                

            }
            
            //already stored all rules in best_rule_in_one_parse
            ConvertToSymbols(best_rule_in_one_parse, out);

        }
        else
            out << "\n";

        
    }


}

/*
This function calculates the recall of the learned grammar given the test set
Param:
    learned_tree_path: the path to learned grammar learned_tree.txt
    test_set: the path to the test dataset
*/
double Recall( string learned_tree_path,  string test_set_path)
{
    //read in test set
    vector<SequenceType<string>> test_set = FileParser<string>(test_set_path);

    //create parser
    shared_ptr<T_AOG<string> > graph_ptr = ReconstructTree(learned_tree_path);
    vector<Symbolic_Rule<string> > all_rules = graph_ptr->GetRules();
    Symbolic_State<string> root_state = graph_ptr -> GetStateByVertexId(graph_ptr ->GetRoot());
    shared_ptr<grammar<string> > g = make_shared<grammar<string> >(all_rules,
                                                            vector<Symbolic_State<string> >{root_state});
    

    double parsed = 0;
    //parse each input 
    for(int i = 0; i < test_set.size(); ++i)
    {
        std::cerr << "[data " << i << "]" << std::endl;
        shared_ptr<EarleyParser<string> > parser = make_shared<EarleyParser<string> >(*g);
        bool parsing_success = false;
        int pos = parser -> parse(test_set[i].begin(), test_set[i].end(),std::cout,parsing_success, false);
        // if((!parsing_success && pos == test_set[i].size()) || (parsing_success && pos != test_set[i].size()))
        // {
        //     std::cerr<<"the parser assumption goes wrong in metrics!\n";
        //     throw std::exception();
        // }
        if(pos == test_set[i].size())
        {
            // std::cerr<<"the parsed seq: \n";
            // for(auto state : test_set[i])
            // {
            //     std::cerr<<"("<<state.GetContent()<<", "<<state.GetId() << ")";
            // }
            // std::cerr<<std::endl;
            ++parsed;

        }
        // else
        // {
        //     std::cerr<<"seq unparsed: \n";
        //     for(auto state : test_set[i])
        //     {
        //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId() << ")";
        //     }
        //     std::cerr<<std::endl;

        //     std::cerr<< "the parsed pos: "<<pos<<std::endl;
        //     int a;
        //     std::cin >>a;
        // }
        // cout << "parsed number: "<<parsed<<endl;
      
    }
    
    cout <<"parsed: "<<parsed<<"\n";
    cout <<"test set size: "<<test_set.size()<<"\n";
    return parsed / test_set.size();

}


shared_ptr<T_AOG<string> > ReconstructTree(string filename)
{
    ifstream f(filename, std::ifstream::in);

    if (f.is_open()) {
	    std::cerr << "reading " << filename << "\n";
    }
    else{
        std::cerr << "file cannot be found" << endl;
        exit(1);
    }
    
    vector<Symbolic_Rule<string> > rules;

    string s = "";
    int line_num = 0;
    unordered_map<Symbolic_Rule<string>, double> m;

    while (getline(f, s))
    {
        line_num++;
        stringstream ss(s);
        string item;
        vector<string> tokens;
        while (getline(ss, item, ',')) 
            tokens.push_back(item);

        if (tokens.size() < 4)
        {
            cerr << "Error: row does not have source id and content\n";
            throw exception();
        }

        

        int num_res = stoi(tokens[0]);
        double weight = stod(tokens[1]);
        int src_id = stoi(tokens[2]);
        string src_ct = tokens[3];
        Symbolic_State<string> src(src_id);

        vector<Symbolic_State<string> > results;
        for (int i = 4; i < tokens.size(); i+=2)
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
        if(weight != 0)
            m[rule] = weight;
    }

    // cout << "All rules: \n";
    // for(auto rule : rules)
    // {
    //     cout <<"("<<rule.GetSource().GetContent()<<", "<<rule.GetSource().GetId()<<") -> ";
    //     for(auto result : rule.GetResults())
    //     {
    //         cout <<"("<<result.GetContent()<<", "<<result.GetId()<<") ";
    //     }
    //     cout <<endl;
    // }

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
        std::cerr<<"top_level_rule: \n";
        for(const auto & state : top_level_rules)
        {
            std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
        }
        std::cerr<<std::endl;
        cerr << "Dangling top level rules beside root\n";
        exit(1);
    }

    shared_ptr<T_AOG<string> > related_graph = make_shared<T_AOG<string> >();
    shared_ptr<AOG_Vertex<string> > source_ptr(new AOG_Vertex<string>(*top_level_rules.begin(), true, true));
	related_graph->AddVertex(source_ptr);

    for (const auto & rule : rules){
        related_graph->AddRule(rule);
    }

    for(const auto& it : m){
        const Symbolic_Rule<string>& rule = it.first;
        double weight = it.second;
        VertexId srcId = related_graph->GetVertexIdByState(rule.GetSource());
        vector<VertexId> resIds;
        for (const Symbolic_State<string>& state : rule.GetResults())
            resIds.push_back(related_graph->GetVertexIdByState(state));
        
        bool found = false;
        VertexId targetDummy;
        for(VertexId dummy : related_graph->ChildrenVertices(srcId)){
            vector<VertexId> children = related_graph->ChildrenVertices(dummy);
            if(children == resIds){
                found = true;
                targetDummy = dummy;
                break;
            }
        }
        if(!found){
            cerr << "rule not found in graph !!" << endl;
            throw exception();
        }
        auto weights = related_graph->GetOutEdgeWeights(srcId, false);
        weights[targetDummy] = weight;
        related_graph->SetOutEdgeWeights(srcId, weights);
    }
    // related_graph->OutputGraph("visualize.txt",PATH,true, false);
	return related_graph;
}

shared_ptr<T_AOG<string> > ReconstructTree_UniformOrNode(string filename)
{
    ifstream f(filename, std::ifstream::in);

    if (f.is_open()) {
	    std::cerr << "reading " << filename << "\n";
    }
    else{
        std::cerr << "file cannot be found" << endl;
        exit(1);
    }
    
    vector<Symbolic_Rule<string> > rules;

    string s = "";
    int line_num = 0;
    // unordered_map<Symbolic_Rule<string>, double> m;

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

        vector<Symbolic_State<string> > results;
        if(s.back() != ',' && tokens.size() == 6){
            Symbolic_State<string> state;
            int state_id = stoi(tokens[3]);
            if (state_id == -1)
            {
                string state_ct = tokens[4];            
                Symbolic_State<string> temp(state_ct, true);
                state = temp;
            }
            else
            {
                Symbolic_State<string> temp(state_id);
                state = temp;
            }
            results.push_back(state);
            Symbolic_Rule<string> rule(src, results);
            // m[rule] = stod(tokens[5]);
            rules.push_back(rule);
        }
        else{
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
    }

    // cout << "All rules: \n";
    // for(auto rule : rules)
    // {
    //     cout <<"("<<rule.GetSource().GetContent()<<", "<<rule.GetSource().GetId()<<") -> ";
    //     for(auto result : rule.GetResults())
    //     {
    //         cout <<"("<<result.GetContent()<<", "<<result.GetId()<<") ";
    //     }
    //     cout <<endl;
    // }

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
        std::cerr<<"top_level_rule: \n";
        for(const auto & state : top_level_rules)
        {
            std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
        }
        std::cerr<<std::endl;
        cerr << "Dangling top level rules beside root\n";
        exit(1);
    }

    shared_ptr<T_AOG<string> > related_graph = make_shared<T_AOG<string> >();
    shared_ptr<AOG_Vertex<string> > source_ptr(new AOG_Vertex<string>(*top_level_rules.begin(), true, true));
	related_graph->AddVertex(source_ptr);

    for (const auto & rule : rules){
        related_graph->AddRule(rule);
    }

    // for (auto rule : rules){
    //     if(m.find(rule) != m.end()){
    //         // cout << "source: " << rule.GetSource().GetId() << endl;
    //         // cout << "result: ";
    //         // for (auto result : rule.GetResults()){
    //         //     cout << "id: " << result.GetId();
    //         // }
    //         // cout << endl;
    //         VertexId srcId = related_graph->GetVertexIdByState(rule.GetSource());
    //         VertexId resId = related_graph->GetVertexIdByState(rule.GetResults()[0]);
    //         VertexId targetDummy;
    //         for(VertexId dummy : related_graph->ChildrenVertices(srcId)){
    //             auto children = related_graph->ChildrenVertices(srcId);
    //             if(related_graph->ChildrenVertices(dummy)[0] == resId){
    //                 targetDummy = dummy;
    //                 break;
    //             }
    //         }
    //         auto weights = related_graph->GetOutEdgeWeights(srcId, false);
    //         weights[targetDummy] = m[rule];
    //         //related_graph->SetOutEdgeWeights(srcId, weights);
    //     }
    // }
    // related_graph->OutputGraph("visualize.txt",PATH,true);
	return related_graph;
}

double Perplexity_cross_entropy(string learned_tree_path, string test_set_path, double penalty){
    
    vector<SequenceType<string>> test_set = FileParser<string>(test_set_path);

    //gernerate learned grammar from learned tree
    shared_ptr<T_AOG<string> > graph_ptr = ReconstructTree(learned_tree_path);
    vector<Symbolic_Rule<string> > all_rules = graph_ptr->GetRules();
    Symbolic_State<string> root_state = graph_ptr -> GetStateByVertexId(graph_ptr ->GetRoot());
    shared_ptr<grammar<string> > g = make_shared<grammar<string> >(all_rules, vector<Symbolic_State<string> >{root_state});                                        



    //stats about test set
    const double type_1_error_penalty = penalty;
    unsigned test_set_size = test_set.size();
    // unsigned unparsed = 0;
    // for(int i = 0; i < test_set_size; ++i){
    //     int pos = Parse(test_set[i], all_rules, parser);
    //     if(pos != test_set[i].size())
    //         ++unparsed;
    //     cout << "unparsed number: "<<unparsed<<endl;  
    // }
    // double smooth_down_ratio = 1 - type_1_error_penalty * unparsed;

    unordered_map<SequenceType<string>, double, container_hash<SequenceType<string>>> test_samples_probabilities;
    for(SequenceType<string> test_sample : test_set){
        if(test_samples_probabilities.find(test_sample) == test_samples_probabilities.end()){
            test_samples_probabilities[test_sample] = 1.0/test_set_size;
        }
        else{
            test_samples_probabilities[test_sample] += 1.0/test_set_size;
        }
    }
    cout << "finished test set stats" << endl;



    // //calculate distribution of sequence that can be generated from the grammar
    // unordered_map<list<VertexId>, double, container_hash<list<VertexId>>> m;
    // queue<list<VertexId>> q;

    // list<VertexId> root_list = list<VertexId>{graph_ptr->GetRoot()};
    // m[root_list] = 1.0;
    // q.push(root_list);

    // unsigned long counter = 0;
    // while(!q.empty()){

    //     cout << "iteration: " << ++counter << endl;
    //     list<VertexId> curr = q.front();
    //     list<VertexId> curr_unchanged = curr;
    //     q.pop();
    //     double currProb = m[curr];

    //     //replace all the and nodes
    //     bool andNodeChanged = false;
    //     for(auto && it = curr.begin(); it != curr.end(); ){
    //         if(!graph_ptr->GetStateByVertexId(*it).GetIsBasic() && graph_ptr->GetVertexContent(*it)->IsAnd()){
    //             vector<VertexId> children = graph_ptr->ChildrenVertices(*it);
    //             it = curr.erase(it);
    //             it = curr.insert(it, children.begin(), children.end());
    //             //advance(it, children.size());
    //             andNodeChanged = true;
    //         }
    //         else{
    //             it++;
    //         }
    //     }

    //     //find or node
    //     bool orNodeChanged = false;
    //     for(auto && it = curr.begin(); it != curr.end(); ){
    //         if(!graph_ptr->GetStateByVertexId(*it).GetIsBasic() && !graph_ptr->GetVertexContent(*it)->IsAnd()){
    //             auto outEdgesProb = graph_ptr->GetOutEdgeWeights(*it, true);
    //             for(VertexId child : graph_ptr->ChildrenVertices(*it)){
    //                 double orProb = outEdgesProb[child];
    //                 it = curr.erase(it);
    //                 it = curr.insert(it, child);
    //                 if(m.find(curr) != m.end()){
    //                     m[curr] += orProb*currProb;
    //                 }
    //                 else{
    //                     m[curr] = orProb*currProb;
    //                 }
    //                 q.push(curr);
    //             }
    //             orNodeChanged = true;
    //             break;
    //         }
    //         else{
    //             it++;
    //         }
    //     }

    //     if(andNodeChanged || orNodeChanged){
    //         m.erase(curr_unchanged);
    //     }
    //     if(andNodeChanged && !orNodeChanged){
    //         if(m.find(curr) != m.end()){
    //             m[curr] += currProb;
    //         }
    //         else{
    //             m[curr] = currProb;
    //         }
    //         q.push(curr);
    //     }

    // }
    // cout << "finished construct distributions of all possible outputs of the model" << endl;



    // //test
    // double total_prob = 0;
    // for(const auto& x : m){
    //     total_prob += x.second;
    // }
    // cout << "total_prob: " << total_prob << endl;



    //calculate perplexity
    unsigned count = 0;
    unsigned count_unparsed = 0;
    double cross_entropy = 0;
    vector<double> test_prob_vector;
    for(const auto& x : test_samples_probabilities){
        cout << "[data " << ++count << "]" << endl;
        
        shared_ptr<EarleyParser<string> > parser = make_shared<EarleyParser<string> >(*g);

        bool parsing_success = false;
        int pos = parser -> parse(x.first.begin(), x.first.end(),std::cout,parsing_success);
        // if((!parsing_success && pos == x.first.size()) || (parsing_success && pos != x.first.size()))
        // {
        //     std::cerr<<"the parser assumption goes wrong in metrics!\n";
        //     throw std::exception();
        // }
        if (pos != x.first.size()){
            cerr << "caution! cannot be parsed!!!!!!!" << endl;
            cross_entropy += -x.second * log2(type_1_error_penalty);
            count_unparsed++;
        }
        else{
            cerr << "parsable !!!!!!!!!" << endl;
            double prob = parser->prob(graph_ptr);
            if (prob == 0){
                std::cerr << "parsable data cannot be of 0 probability" << std::endl;
                throw std::exception();
            }
            cross_entropy += -x.second * log2(prob);
            cerr << "prob: " << prob << endl;
            test_prob_vector.push_back(x.second);
        }
    }
    double smooth_down_ratio = 1 - count_unparsed * type_1_error_penalty;
    for(double test_prob : test_prob_vector){
        cross_entropy += -test_prob * log2(smooth_down_ratio);
    }

    cout << "count_unparsed: " << count_unparsed << endl;
    cout << "smooth_down_ratio: " << smooth_down_ratio << endl;
    cout << "cross_entropy: " << cross_entropy << endl;

    return exp2(cross_entropy);

}

double KL_Divergence_Train_Set(string learned_tree_path, string train_set_path){
    
    vector<SequenceType<string>> train_set = FileParser<string>(train_set_path);

    //gernerate learned grammar from learned tree
    shared_ptr<T_AOG<string> > graph_ptr = ReconstructTree(learned_tree_path);


    //stats about test set
    unsigned train_set_size = train_set.size();
    unordered_map<SequenceType<string>, double, container_hash<SequenceType<string>>> train_samples_probabilities;
    for(const SequenceType<string> &test_sample : train_set){
        if(train_samples_probabilities.find(test_sample) == train_samples_probabilities.end()){
            train_samples_probabilities[test_sample] = 1.0/train_set_size;
        }
        else{
            train_samples_probabilities[test_sample] += 1.0/train_set_size;
        }
    }
    cout << "finished test set stats" << endl;



    //calculate distribution of sequence that can be generated from the grammar
    unordered_map<list<VertexId>, double, container_hash<list<VertexId>>> m;
    queue<list<VertexId>> q;

    list<VertexId> root_list = list<VertexId>{graph_ptr->GetRoot()};
    m[root_list] = 1.0;
    q.push(root_list);

    unsigned long counter = 0;
    while(!q.empty()){

        cout << "iteration: " << ++counter << endl;
        list<VertexId> curr = q.front();
        list<VertexId> curr_unchanged = curr;
        q.pop();
        double currProb = m[curr];

        //replace all the and nodes
        bool andNodeChanged = false;
        for(auto it = curr.begin(); it != curr.end(); ){
            if(!graph_ptr->GetStateByVertexId(*it).GetIsBasic() && graph_ptr->GetVertexContent(*it)->IsAnd()){
                vector<VertexId> children = graph_ptr->ChildrenVertices(*it);
                it = curr.erase(it);
                it = curr.insert(it, children.begin(), children.end());
                //advance(it, children.size());
                andNodeChanged = true;
            }
            else{
                it++;
            }
        }

        //find or node
        bool orNodeChanged = false;
        for(auto it = curr.begin(); it != curr.end(); ){
            if(!graph_ptr->GetStateByVertexId(*it).GetIsBasic() && !graph_ptr->GetVertexContent(*it)->IsAnd()){
                auto outEdgesProb = graph_ptr->GetOutEdgeWeights(*it, true);
                for(VertexId child : graph_ptr->ChildrenVertices(*it)){
                    double orProb = outEdgesProb[child];
                    it = curr.erase(it);
                    it = curr.insert(it, child);
                    if(m.find(curr) != m.end()){
                        m[curr] += orProb*currProb;
                    }
                    else{
                        m[curr] = orProb*currProb;
                    }
                    q.push(curr);
                }
                orNodeChanged = true;
                break;
            }
            else{
                it++;
            }
        }

        if(andNodeChanged || orNodeChanged){
            m.erase(curr_unchanged);
        }
        if(andNodeChanged && !orNodeChanged){
            if(m.find(curr) != m.end()){
                m[curr] += currProb;
            }
            else{
                m[curr] = currProb;
            }
            q.push(curr);
        }

    }
    cout << "finished construct distributions of all possible outputs of the model" << endl;



    //test
    double total_prob = 0;
    for(const auto& x : m){
        total_prob += x.second;
    }
    cout << "total_prob: " << total_prob << endl;



    //calculate perplexity
    unsigned unparsed_count = 0;
    double kl = 0;
    for(const auto& x : train_samples_probabilities){
        list<VertexId> x_vertexids;
        for (const auto& it : x.first){
            x_vertexids.push_back(graph_ptr->GetVertexIdByState(it));
        }
        if(m.find(x_vertexids) == m.end()){
            cerr << "data in training set is not parsable !!!" << endl;
            unparsed_count ++;
        }
        else{
            double prob_train = x.second;
            double prob_learned = m[x_vertexids];
            kl += prob_train * log2(prob_train / prob_learned);
        }
    }

    cout << "unparsed count: " << unparsed_count << endl;

    return kl;

}


double KL_Divergence_Ground_Truth_Tree(string learned_tree_path, string ground_truth_path){

    auto GenerateAllData = [](shared_ptr<T_AOG<string>> graph_ptr, unordered_set<list<VertexId>, container_hash<list<VertexId> > >& s, unordered_map<list<VertexId>, double, container_hash<list<VertexId>>>& m){
        unsigned long counter = 0;
        while(!s.empty()){

            cerr << endl;
            cout << "iteration: " << ++counter << endl;
            cerr << "size of set: " << s.size() << endl;
            cerr << "size of map: " << m.size() << endl << "------" << endl;
            // cerr << "total prob: " << accumulate(m.begin(), m.end(), 0.0, [](const double prev, const pair<list<VertexId>, double> &p){return prev + p.second;}) << endl;
            list<VertexId> curr = *(s.begin());
            list<VertexId> curr_unchanged = curr;
            s.erase(s.begin());

            // if(curr.size() >= 11){
            //     cerr << "longer than 11..." << endl;
            //     for(VertexId x : curr){
            //         cerr << x << " " ;
            //     }
            //     cout << endl;
            //     for(VertexId x: curr){
            //         Symbolic_State<string> state = graph_ptr->GetStateByVertexId(x);
            //         cerr << "( " << state.GetId() << ", " << state.GetContent() << " )" << endl;
            //         if(!state.GetIsBasic()){
            //             for (VertexId child : graph_ptr->ChildrenVertices(x)){
            //                 Symbolic_State<string> child_state = graph_ptr->GetStateByVertexId(child);
            //                 cerr << "\t( " << child_state.GetId() << ", " << child_state.GetContent() << " )" << endl;
            //             }
            //         }
            //     }
            //     throw exception();
            // }

            // cerr << "curr_unchanged: ";
            // for (VertexId x : curr_unchanged){
            //     cerr << x << ' ';
            // }
            // cerr << endl;
            

            double currProb = m[curr];

            // cerr << "currProb: " << currProb << endl;
            if(currProb == 0){
                cerr << "0 currProb is catched!!!" << endl;
            }

            //replace all the and nodes
            bool andNodeChanged = false;
            for(auto it = curr.begin(); it != curr.end(); ){
                if(!graph_ptr->GetStateByVertexId(*it).GetIsBasic() && graph_ptr->GetVertexContent(*it)->IsAnd()){
                    vector<VertexId> children = graph_ptr->ChildrenVertices(*it);
                    it = curr.erase(it);
                    it = curr.insert(it, children.begin(), children.end());
                    //advance(it, children.size());
                    andNodeChanged = true;
                }
                else{
                    it++;
                }
            }

            // if (andNodeChanged){
            //     cerr << "andNodeChanged: ";
            //     for(VertexId x : curr){
            //         cerr << x << ' ';
            //     }
            //     cerr << endl;
            // }

            //find or node
            bool orNodeChanged = false;
            for(auto it = curr.begin(); it != curr.end(); ){
                if(!graph_ptr->GetStateByVertexId(*it).GetIsBasic() && !graph_ptr->GetVertexContent(*it)->IsAnd()){
                    auto outEdgesProb = graph_ptr->GetOutEdgeWeights(*it, true);
                    for(VertexId child : graph_ptr->ChildrenVertices(*it)){
                        double orProb = outEdgesProb[child];
                        it = curr.erase(it);
                        it = curr.insert(it, child);
                        if(m.find(curr) != m.end()){
                            m[curr] += orProb*currProb;
                        }
                        else{
                            m[curr] = orProb*currProb;
                        }
                        s.insert(curr);

                        // cerr << "[ ";
                        // for(VertexId x : curr){
                        //     cerr << x << " ";
                        // }
                        // cerr << " ]" << endl;
                        // cerr << "after or prob: " << orProb*currProb << endl;
                        
                    }
                    orNodeChanged = true;
                    break;
                }
                else{
                    it++;
                }
            }

            if(andNodeChanged || orNodeChanged){
                m.erase(curr_unchanged);
            }
            if(andNodeChanged && !orNodeChanged){
                if(m.find(curr) != m.end()){
                    m[curr] += currProb;
                }
                else{
                    m[curr] = currProb;
                }
                s.insert(curr);
            }

        }
    };


    //gernerate learned grammar from learned tree
    shared_ptr<T_AOG<string> > graph_ptr_learned = ReconstructTree(learned_tree_path);
    graph_ptr_learned->TruncateGraph();
    // graph_ptr_learned->SimplifyGraph();
    shared_ptr<T_AOG<string> > graph_ptr_ground = ReconstructTree(ground_truth_path);


    //calculate distribution of sequence that can be generated from the learned grammar
    unordered_map<list<VertexId>, double, container_hash<list<VertexId>>> m_learned;
    list<VertexId> root_list_learned = list<VertexId>{graph_ptr_learned->GetRoot()};
    unordered_set<list<VertexId>, container_hash<list<VertexId> > > s_learned {root_list_learned};
    m_learned[root_list_learned] = 1.0;
    GenerateAllData(graph_ptr_learned, s_learned, m_learned);


    //calculate distribution of sequence that can be generated from the ground true grammar
    unordered_map<list<VertexId>, double, container_hash<list<VertexId>>> m_ground;
    list<VertexId> root_list_ground = list<VertexId>{graph_ptr_ground->GetRoot()};
    unordered_set<list<VertexId>, container_hash<list<VertexId> > > s_ground {root_list_ground};
    m_ground[root_list_ground] = 1.0;
    GenerateAllData(graph_ptr_ground, s_ground, m_ground);


    //test
    double total_prob = 0;
    for(const auto& x : m_learned){
        total_prob += x.second;
    }
    cout << "total_prob of learned tree: " << total_prob << endl;

    total_prob = 0;
    for(const auto& x : m_ground){
        total_prob += x.second;
    }
    cout << "total_prob of ground truth: " << total_prob << endl;


    //calculate perplexity
    unsigned unparsed_count = 0;
    double kl = 0;
    for(const auto& x : m_learned){
        list<VertexId> l;
        for(VertexId id : x.first){
            l.push_back(graph_ptr_ground->GetVertexIdByState(Symbolic_State<string>(graph_ptr_learned->GetStateByVertexId(id).GetContent(), true)));
        }
        if(m_ground.find(l) == m_ground.end()){
            cerr << "data generated in learned_tree is not parsable in ground truth !!!" << endl;
            unparsed_count ++;
        }
        else{
            double prob_learned = x.second;
            double prob_ground = m_ground[l];
            kl += prob_learned * log2(prob_learned / prob_ground);
        }
    }

    cout << "total size of data generate-able from learned_tree: " << m_learned.size() << endl;

    cout << "unparsed count: " << unparsed_count << endl;

    return kl;

}

double JS_Divergence_By_Sampling(string learned_tree, string ground_truth_tree, unsigned sample_time){


    shared_ptr<T_AOG<string> > learned_graph_ptr = ReconstructTree(learned_tree);
    vector<Symbolic_Rule<string> > all_rules = learned_graph_ptr->GetRules();
    Symbolic_State<string> root_state = learned_graph_ptr -> GetStateByVertexId(learned_graph_ptr->GetRoot());
    shared_ptr<grammar<string> > g = make_shared<grammar<string> >(all_rules, vector<Symbolic_State<string> >{root_state});
    
    //sample from the learned AOG and calculate precision
    unordered_map<vector<VertexId>, double, container_hash<vector<VertexId> > > samples_map;

    double total_prob = 0;

    while(sample_time != 0)
    {
        --sample_time;
        vector<VertexId> seq;
        double prob;
        shared_ptr<vector<VertexId>> parse_tree = learned_graph_ptr->Sample(learned_graph_ptr->GetRoot(), seq, prob);

        if(samples_map.find(seq) == samples_map.end()){

            vector<Symbolic_State<string> > seq_state;
            for(const auto& state : seq){
                seq_state.push_back(learned_graph_ptr->GetStateByVertexId(state));
            }
            shared_ptr<EarleyParser<string> > parser = make_shared<EarleyParser<string> >(*g);
            bool parsing_success = false;
            int pos = parser->parse(seq_state.begin(), seq_state.end(), std::cout, parsing_success);
            if (pos != seq_state.size()){
                std::cerr << "sampled data not parsable" << std::endl;
                throw exception();
            }
            double prob = parser->prob(learned_graph_ptr);
            samples_map[seq] = prob;
            total_prob += prob;
            
        }
    }


    cerr << "total prob: " << total_prob << endl;

    
    unordered_map<vector<VertexId>, double, container_hash<vector<VertexId> > > joint_map;

    for(const auto& pair : samples_map){
        joint_map[pair.first] = pair.second / 2;
    }

    auto GenerateAllData = [](shared_ptr<T_AOG<string>> graph_ptr, unordered_set<list<VertexId>, container_hash<list<VertexId> > >& s, unordered_map<list<VertexId>, double, container_hash<list<VertexId>>>& m){
        unsigned long counter = 0;
        while(!s.empty()){

            // cerr << endl;
            // cout << "iteration: " << ++counter << endl;
            // cerr << "size of set: " << s.size() << endl;
            // cerr << "size of map: " << m.size() << endl << "------" << endl;
            // cerr << "total prob: " << accumulate(m.begin(), m.end(), 0.0, [](const double prev, const pair<list<VertexId>, double> &p){return prev + p.second;}) << endl;
            list<VertexId> curr = *(s.begin());
            list<VertexId> curr_unchanged = curr;
            s.erase(s.begin());

            // if(curr.size() >= 11){
            //     cerr << "longer than 11..." << endl;
            //     for(VertexId x : curr){
            //         cerr << x << " " ;
            //     }
            //     cout << endl;
            //     for(VertexId x: curr){
            //         Symbolic_State<string> state = graph_ptr->GetStateByVertexId(x);
            //         cerr << "( " << state.GetId() << ", " << state.GetContent() << " )" << endl;
            //         if(!state.GetIsBasic()){
            //             for (VertexId child : graph_ptr->ChildrenVertices(x)){
            //                 Symbolic_State<string> child_state = graph_ptr->GetStateByVertexId(child);
            //                 cerr << "\t( " << child_state.GetId() << ", " << child_state.GetContent() << " )" << endl;
            //             }
            //         }
            //     }
            //     throw exception();
            // }

            // cerr << "curr_unchanged: ";
            // for (VertexId x : curr_unchanged){
            //     cerr << x << ' ';
            // }
            // cerr << endl;
            

            double currProb = m[curr];

            // cerr << "currProb: " << currProb << endl;
            if(currProb == 0){
                cerr << "0 currProb is catched!!!" << endl;
            }

            //replace all the and nodes
            bool andNodeChanged = false;
            for(auto it = curr.begin(); it != curr.end(); ){
                if(!graph_ptr->GetStateByVertexId(*it).GetIsBasic() && graph_ptr->GetVertexContent(*it)->IsAnd()){
                    vector<VertexId> children = graph_ptr->ChildrenVertices(*it);
                    it = curr.erase(it);
                    it = curr.insert(it, children.begin(), children.end());
                    //advance(it, children.size());
                    andNodeChanged = true;
                }
                else{
                    it++;
                }
            }

            // if (andNodeChanged){
            //     cerr << "andNodeChanged: ";
            //     for(VertexId x : curr){
            //         cerr << x << ' ';
            //     }
            //     cerr << endl;
            // }

            //find or node
            bool orNodeChanged = false;
            for(auto it = curr.begin(); it != curr.end(); ){
                if(!graph_ptr->GetStateByVertexId(*it).GetIsBasic() && !graph_ptr->GetVertexContent(*it)->IsAnd()){
                    auto outEdgesProb = graph_ptr->GetOutEdgeWeights(*it, true);
                    for(VertexId child : graph_ptr->ChildrenVertices(*it)){
                        double orProb = outEdgesProb[child];
                        it = curr.erase(it);
                        it = curr.insert(it, child);
                        if(m.find(curr) != m.end()){
                            m[curr] += orProb*currProb;
                        }
                        else{
                            m[curr] = orProb*currProb;
                        }
                        s.insert(curr);

                        // cerr << "[ ";
                        // for(VertexId x : curr){
                        //     cerr << x << " ";
                        // }
                        // cerr << " ]" << endl;
                        // cerr << "after or prob: " << orProb*currProb << endl;
                        
                    }
                    orNodeChanged = true;
                    break;
                }
                else{
                    it++;
                }
            }

            if(andNodeChanged || orNodeChanged){
                m.erase(curr_unchanged);
            }
            if(andNodeChanged && !orNodeChanged){
                if(m.find(curr) != m.end()){
                    m[curr] += currProb;
                }
                else{
                    m[curr] = currProb;
                }
                s.insert(curr);
            }

        }
    };

    shared_ptr<T_AOG<string> > graph_ptr_ground = ReconstructTree(ground_truth_tree);
    graph_ptr_ground->TruncateGraph();

    //calculate distribution of sequence that can be generated from the ground true grammar
    unordered_map<list<VertexId>, double, container_hash<list<VertexId> > > m_ground;
    list<VertexId> root_list_ground = list<VertexId>{graph_ptr_ground->GetRoot()};
    unordered_set<list<VertexId>, container_hash<list<VertexId> > > s_ground {root_list_ground};
    m_ground[root_list_ground] = 1.0;
    GenerateAllData(graph_ptr_ground, s_ground, m_ground);

    unordered_map<vector<VertexId>, double, container_hash<vector<VertexId> > > ground_map;

    for(const auto& pair : m_ground){
        vector<VertexId> v(pair.first.begin(), pair.first.end());
        ground_map[v] = pair.second;
        joint_map[v] += pair.second / 2;
    }

    double js_divergence = 0;
    for(const auto& pair : samples_map){
        js_divergence += pair.second * log2(pair.second / joint_map.at(pair.first));
    }
    for(const auto& pair : ground_map){
        js_divergence += pair.second * log2(pair.second / joint_map.at(pair.first));
    }
    js_divergence /= 2;

    return js_divergence;
}

double Bleu (string reference_file, string hypothesis_file){
    istream_iterator<string> eos;
    istreambuf_iterator<char> eoc;
    
    // clock_t start;

    // ifstream reference_stream_tmp(reference_file);
    ifstream reference_stream(reference_file);
    // size_t ref_line_count = count(istreambuf_iterator<char>(reference_stream_tmp), eoc, '\n') + 1;
    
    vector<vector<string>> reference;
    // reference.reserve(ref_line_count);
    unordered_set<string> ref_set;
    // ref_set.reserve(ref_line_count);
    
    if(reference_stream.is_open()){
        string line;
        while(getline(reference_stream, line)){
            if(ref_set.find(line) == ref_set.end()){
                ref_set.insert(line);
                istringstream iss(line);
                reference.emplace_back(istream_iterator<string>(iss), eos);
            }
        }
    }
    else{
        cerr << "data set not found for reference" << endl;
        throw exception();
    }
    // if(reference.empty()){
    //     cerr << "empty reference" << endl;
    //     throw exception();
    // }


    // ifstream hypothesis_stream_tmp(hypothesis_file);
    ifstream hypothesis_stream(hypothesis_file);
    // size_t hyp_line_count = count(istreambuf_iterator<char>(hypothesis_stream_tmp), eoc, '\n') + 1;
    
    vector<vector<string>> hypothesis;
    // hypothesis.reserve(hyp_line_count);
    vector<double> hyp_likelihood;
    // hyp_likelihood.reserve(hyp_line_count);
    
    if(hypothesis_stream.is_open()){
        string line;
        while(getline(hypothesis_stream, line)){
            size_t index;
            double likelihood = stod(line, &index);
            hyp_likelihood.push_back(likelihood);
            istringstream iss(line.substr(index + 2));
            hypothesis.emplace_back(istream_iterator<string>(iss), eos);
        }
    }
    else{
        cerr << "sample data file not found" << endl;
        throw exception();
    }
    // if(hypothesis.empty() || hyp_likelihood.empty()){
    //     cerr << "empty hypothesis or hypothesis likelihood" << endl;
    //     throw exception();
    // }

    // double timenow = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);

    double total_likelihood = 0.0;
    for_each(hyp_likelihood.begin(), hyp_likelihood.end(), [&total_likelihood](double tmp){total_likelihood += tmp;});


    double bleu = 0.0;
    long hyp_index = -1;
    for(const auto& hyp : hypothesis){
        ++hyp_index;

        // for brevity penalty, refering to brevity_penalty in https://www.nltk.org/_modules/nltk/translate/bleu_score.html
        size_t hyp_len = hyp.size();

        pair<unsigned, size_t> absDis_refSize_pair(abs((signed)reference[0].size() - (signed)hyp_len), reference[0].size());
        for(auto& it = ++reference.begin(); it != reference.end(); ++it){
            size_t ref_size = (*it).size();
            unsigned diff = abs((signed)ref_size - (signed)hyp_len);
            if(diff <= absDis_refSize_pair.first && ref_size < absDis_refSize_pair.second){
                absDis_refSize_pair.first = diff;
                absDis_refSize_pair.second = ref_size;
            }
        }

        double brevity_penalty;
        unsigned closed_ref_length = absDis_refSize_pair.second;
        if (hyp_len > closed_ref_length){
            brevity_penalty = 1.0;
        }
        else{
            brevity_penalty = exp(1 -  (double)closed_ref_length / hyp_len);
        }


        size_t len = hyp_len < 4 ? hyp_len : 4;
        
        double hyp_score = 0.0;
        //c++11: If there is no key equivalent to x in the map, inserts value_type(x, T()) into the map.

        for(size_t n_gram = 1; n_gram <= len; n_gram++){
            unordered_map<vector<string>, size_t, container_hash<vector<string> > > counts;
            for(size_t i = 0; i <= hyp_len - n_gram; i++){
                vector<string> word(n_gram);
                for(size_t j = 0; j < n_gram; j++){
                    word[j] = hyp[i + j];
                }
                counts[word]++;
            }
            unordered_map<vector<string> , size_t, container_hash<vector<string> > > clipped_counts;
            for(const auto& ref : reference){
                for(const auto& count : counts){
                    size_t current_clipped_count = 0;
                    auto it = ref.begin();
                    for(int i = 0; i <= (signed)ref.size() - (signed)n_gram; i++){ // must use signed type, ref.size() - n_gram can be -ve
                        if (equal(count.first.begin(), count.first.end(), it))
                            current_clipped_count++;
                        it++;
                    }
                    if(current_clipped_count > clipped_counts[count.first]){
                        clipped_counts[count.first] = current_clipped_count;
                    }
                }
            }

            size_t numerator = 0;
            size_t denominator = 0;
            
            for(const auto& count : counts){
                denominator += count.second;
                unsigned clipped_count = clipped_counts[count.first]; 
                if(clipped_count > count.second)
                    numerator += count.second;
                else
                    numerator += clipped_count;
            }

            if(denominator == 0){
                denominator = 1; // to avoid devision by 0 error
            }

            // if(numerator == 0){
            //     hyp_score = 0.0;
            //     break;
            // }
            // else{
                hyp_score += (1.0/len) * log(((double)numerator / denominator));
            // }
        }
        
        // if(hyp_score != 0.0){
            hyp_score = brevity_penalty * exp(hyp_score);
        // }
        bleu += hyp_likelihood[hyp_index] * hyp_score;
    }

    reference_stream.close();
    hypothesis_stream.close();

    bleu /= total_likelihood;

    return bleu;
}

namespace AOG_LIB_UTIL
{
// sample with total prob given by parse tree
void Sample_Total ( string learned_tree_path, unsigned sample_time )
{
    shared_ptr<T_AOG<string> > learned_graph_ptr = ReconstructTree(learned_tree_path);
    vector<Symbolic_Rule<string> > all_rules = learned_graph_ptr->GetRules();
    Symbolic_State<string> root_state = learned_graph_ptr -> GetStateByVertexId(learned_graph_ptr->GetRoot());
    shared_ptr<grammar<string> > g = make_shared<grammar<string> >(all_rules, vector<Symbolic_State<string> >{root_state});
    
    //sample from the learned AOG and calculate precision
    unordered_set<vector<VertexId>, container_hash<vector<VertexId> > > samples_map;
    // set<vector<VertexId> > sampled_dataset;

    container_hash<vector<VertexId> > hash_v;

    double total_prob = 0;

    while(sample_time != 0)
    {
        --sample_time;
        vector<VertexId> seq;
        double prob;
        shared_ptr<vector<VertexId>> parse_tree = learned_graph_ptr->Sample(learned_graph_ptr->GetRoot(), seq, prob);

        if(samples_map.find(seq) == samples_map.end()){

            samples_map.insert(seq);
            vector<Symbolic_State<string> > seq_state;
            for(const auto& state : seq){
                seq_state.push_back(learned_graph_ptr->GetStateByVertexId(state));
            }
            shared_ptr<EarleyParser<string> > parser = make_shared<EarleyParser<string> >(*g);
            bool parsing_success = false;
            int pos = parser->parse(seq_state.begin(), seq_state.end(), std::cout, parsing_success);
            if (pos != seq_state.size()){
                std::cerr << "sampled data not parsable" << std::endl;
                throw exception();
            }
            double prob = parser->prob(learned_graph_ptr);
            total_prob += prob;

            cout << prob << ": ";
            for(const auto & state : seq_state)
                cout << state.GetContent() << " " ;
            cout << endl;
        
            
        }
    }


    cerr << "total prob: " << total_prob << endl;
}
}



#endif
