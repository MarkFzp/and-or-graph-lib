int FILE_ID = 10;
int AND_STATE_1 = -1;
int AND_STATE_2 = -1;
// TODO: check SN_Ptr point to same address(whether this is indeed address equal)
struct find_checkpoint
{
    SN_Ptr target;

    find_checkpoint(SN_Ptr target) : target(target)
    {}

    bool operator()(const SN_Ptr &candidate)
    {
        return &candidate == &target;
    }
};

int Combination(int n, int k);


template<class StateType>
RdBuffer<StateType> Online_Learner<StateType>::Checkpoint::GenerateCMAndReductions(const ConfigBuffer<StateType> &sorted_configs, 
                                                                                    const ConfigBuffer<StateType>& reduction_space, 
                                                                                    bool generate_cm)
{
    if(generate_cm)
        boost::multi_index::get<3>(this->CM_).clear();
    
    //construct mapping between frequency level and sorted configs
    int config_len = sorted_configs[0].size();
    int num_of_configs = sorted_configs.size();
    RdBuffer<StateType> reductions;

    // std::vector<int> freq_level(num_of_configs);
    // std::vector<int> freq_to_config(num_of_configs);

    // for(int i = 0; i < num_of_configs; ++i)
    // {
    //     freq_level[i] = num_of_configs - i - 1;
    // }
    // std::unordered_map<int, SequenceType<StateType> > idx_to_configs;
    // std::unordered_map<SequenceType<StateType>, int > configs_to_idx;
    // for(int i = num_of_configs ; i >= 1; --i)
    // {
    //     idx_to_configs[i] = sorted_configs[num_of_configs - i];
    //     configs_to_idx[sorted_configs[num_of_configs - i]] = i;
    // }
    
    struct CMP //this is a comparator with an integer map with the same size as the data buffer
	{
		std::vector<std::vector<int> > search_map_;

		CMP(const std::vector<SequenceType<StateType> >& related_data)
		{
			this->search_map_.resize(0);
			for(const auto & data: related_data)
			{
				this->search_map_.push_back(std::vector<int>(data.size(),-1));
			}
		}
		//udpate integer map with the given index to the corresponding config
		void Update(unsigned row, unsigned col, unsigned len, int config_idx)
		{
			for(unsigned i = 0; i < len; ++i)
			this->search_map_[row][col + i] = config_idx;
		}

		//check the boolean map so that lower frequency configs would not overlap with higher frequency configs
		bool Check(unsigned row, unsigned col, unsigned len, int config_idx)
		{
			for(unsigned i = 0; i < len; ++i)
				if(this->search_map_[row][col + i] != -1)
					return false;
			return true;
		}
	};
	CMP cmp(reduction_space);
    
    //iterate through each data to find all possible reductions. Currently implementation not nessarily find the max possible number of reductions
    for (int i = 0; i < reduction_space.size(); ++i) 
    {
        //find all configurations that can be replaced in this data point
        for (int j = 0; j < sorted_configs.size(); ++j)
        {

            auto config_begin = std::search(reduction_space[i].begin(),
                                            reduction_space[i].end(),
                                            sorted_configs[j].begin(),
                                            sorted_configs[j].end());

            //find all occurence of the current configuration in a data point
            while(config_begin != reduction_space[i].end())
            {
                int pos = config_begin - reduction_space[i].begin();
                if(cmp.Check(i,pos,config_len,j))
                {
                    cmp.Update(i,pos,config_len,j);
                }
                
                config_begin = std::search(config_begin + config_len,
                                            reduction_space[i].end(),
                                            sorted_configs[j].begin(),
                                            sorted_configs[j].end());
                
              
            }
        }
      
        //Generate CM and reduction
        if(std::any_of(cmp.search_map_[i].begin(),cmp.search_map_[i].end(),[](int idx){return idx != -1;}))
        {
            //record the indices where search map value changes
            int cur_val = cmp.search_map_[i][0];
            std::vector<int> break_idx(1,0);
            int cur_len = 0;
            for(int k = 0; k < cmp.search_map_[i].size(); ++k)
            {
                if(cmp.search_map_[i][k] != cur_val || (cur_val != -1 && cmp.search_map_[i][k] == cur_val && cur_len == config_len))
                {
                    break_idx.push_back(k);
                    cur_len = 1;
                    cur_val = cmp.search_map_[i][k];
                }
                else
                    ++cur_len;
            }
          
            //if generate both CM and reductions
            if(generate_cm)
            {
                ContextType<StateType> contexts;
                std::vector<SequenceType<StateType> > configurations;
                
                //if there is no context at the beginning, add an empty context            
                if(cmp.search_map_[i][0] != -1)
                    contexts.push_back(SequenceType<StateType>()); 

                //update configurations and contexts found in this data point
                for(int k = 0; k < break_idx.size(); ++ k)
                {
                    //if the current index indicates this fragment is a context
                    if(cmp.search_map_[i][break_idx[k]]== -1)
                    {
                        int context_end;
                        if(k + 1 == break_idx.size())
                            context_end = cmp.search_map_[i].size();                    
                        else    
                            context_end = break_idx[k + 1];
                        
                        contexts.push_back(SequenceType<StateType>(reduction_space[i].begin() + break_idx[k],
                                                                                reduction_space[i].begin() + context_end));
                    }
                    //if the current index indicates this fragment is a configuration
                    else
                    {
                        //if the previous part is also a config, add an empty conetext between them
                        if(break_idx[k] != 0 && cmp.search_map_[i][break_idx[k] -1] != -1 )
                            contexts.push_back(SequenceType<StateType>());
                            
                        
                        int config_idx = cmp.search_map_[i][break_idx[k]];
                        configurations.push_back(sorted_configs[config_idx]);

                        //put the configuration into reductions
                        if(reductions.find(sorted_configs[config_idx]) == reductions.end())
                        {
                            std::vector<std::pair<unsigned, unsigned> > positions(0);
                            reductions[sorted_configs[config_idx]] = positions;
                        }
                        reductions[sorted_configs[config_idx]].push_back(std::pair<unsigned,unsigned>(i,break_idx[k]));
                            
                    }

                }


                //if there is no context at the end, add an empty context
                if(cmp.search_map_[i].back() != -1)
                    contexts.push_back(SequenceType<StateType>());
                
                if(configurations.size() + 1 != contexts.size())
                {
                    std::cout << "configurations size: "<<configurations.size()<<std::endl;
                    std::cout << "contexts size: "<<contexts.size()<<std::endl;
                    std::cout<<"The data: \n";
                    for(const auto & state : reduction_space[i])
                    {
                        std::cout << "("<<state.GetContent()<<", "<<state.GetId()<<")";
                    }
                    std::cout << std::endl;

                    std::cout <<"The configs:\n";
                    for(const auto & config : configurations)
                    {
                        for(const auto & state : config)
                        {
                            std::cout << "("<<state.GetContent()<<", "<<state.GetId()<<")";

                        }
                        std::cout << std :: endl;
                    }
                    std::cout << "The contexts:\n";
                    for(const auto & context : contexts)
                    {
                        for(const auto & state : context)
                        {
                            std::cout << "("<<state.GetContent()<<", "<<state.GetId()<<")";

                        }
                        std::cout << std::endl;
                    }
                }
              
                assert(configurations.size() + 1 == contexts.size());

                //Update CM matrix
                auto &config_idx = boost::multi_index::get<1>(this->CM_);
                auto same_config = config_idx.equal_range(configurations);
                bool entry_exist = false;
                //if same configuration-context setting exists, increases its count
                for (auto& iter = same_config.first; iter != same_config.second; iter++)
                {
                    if ((*iter).context_ == contexts)
                    {
                        config_idx.modify(iter,
                                            [](AOG_LIB::CM_Entry<StateType> &entry)
                                            { ++entry.count_; });
                        entry_exist = true;
                        break;
                    }
                }
                //else create a new CM entry
                if (!entry_exist)
                    this->CM_.insert({contexts, configurations, 1});
            }

            //if only generate reductions
            else
            {
                for(int k = 0; k < break_idx.size(); ++k)
                {
                    int config_idx = cmp.search_map_[i][break_idx[k]];

                    if(config_idx != -1)
                    {
                        if(reductions.find(sorted_configs[config_idx]) == reductions.end())
                        {
                            std::vector<std::pair<unsigned, unsigned> > positions(0);
                            reductions[sorted_configs[config_idx]] = positions;
                        }
                        reductions[sorted_configs[config_idx]].push_back(std::pair<unsigned,unsigned>(i,break_idx[k]));
                    }
                }
            }
            

            
        }    


       
        
        
    }
    return reductions;
    /** std::cerr<<"------------------------Out of GenerateCM()-----------------------------------\n"; */
}

template<class StateType>
SequenceType<StateType> Online_Learner<StateType>::Checkpoint::ParseDelete(const SequenceType<StateType> &reparsed_seq, double &reduced_likelihood)
{
    // std::cerr << "Print the sequence to be reparsed: \n";
    // for(auto && it = reparsed_seq.begin(); it != reparsed_seq.end(); it++)
    // {
    //     std::cerr << "(" << (*it).GetContent() << ", " << (*it).GetId() << ") ";
    // }
    // std::cerr << "\n";
    std::vector<Symbolic_Rule<StateType>> reparse_rules = this->graph_.GetRules();
    for (auto it = reparse_rules.begin(); it != reparse_rules.end(); it++)
    {
        if (it->GetSource() == this->graph_.GetStateByVertexId(this->graph_.GetRoot()))
        {
            it = reparse_rules.erase(it);
            it--;
        }
    }
    std::shared_ptr<grammar<StateType>> t = std::make_shared<grammar<StateType>>(reparse_rules,
                                                                                    this->graph_.GetTopLevelStates(), reparsed_seq);
    // std::cerr << "Print terminals before parseDelete\n";
    // t->print_terminal();
    this->parser_= std::make_shared<EarleyParser<StateType>>(*t);

    // clear position record
    this->parseposition.Clear();

    // std::cerr<<"top level states before parseDelete\n";
    // for(auto state : this->graph_.GetTopLevelStates())
    // {
    //     std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    // }
    // std::cerr<<std::endl;
    // this->graph_.OutputGraph("before_parsedelete", PATH, false);

    this->PartialParseData(reparsed_seq, reparse_rules);
    std::vector<Symbolic_State<StateType>> reparse_replacement = this->GetBestParsing(reduced_likelihood);
    // std::cerr << "Out of GetBestParsing\n";

    // Check if all the nodes in the best parsing exist in the input seq:
    // bool all_in_input = true;
    // for(auto it = reparse_replacement.begin(); it != reparse_replacement.end(); it++)
    // {
    //     if(std::find(reparsed_seq.begin(), reparsed_seq.end(), *it) == reparsed_seq.end())
    //         all_in_input = false;
    // }

    // if(all_in_input)
    // {
    //     reparse_replacement = reparsed_seq;
    // }
    // else

    // std::cerr << "Delete the node if it exists in the input seq:\n";
    // std::cerr << "Sequence before delete: ";
    // for(auto && it=reparse_replacement.begin(); it != reparse_replacement.end(); it++)
    // {
    //     std::cerr <<"(" << (*it).GetContent() << ", " << (*it).GetId() << ") ";
    // }

    // if(reparse_replacement.size() == 2 && reparse_replacement[0].GetId() == 1106 && reparse_replacement[1].GetId() == 1085)
    // {
    //     this->graph_.OutputGraph("debug_graph", PATH, true);
    // }

    // SequenceType<StateType> delete_seq;
    // for(auto it = reparse_replacement.begin(); it != reparse_replacement.end(); it++)
    // {
    //     if(std::find(reparsed_seq.begin(), reparsed_seq.end(), (*it)) != reparsed_seq.end())
    //     {
    //         delete_seq.push_back(*it);
    //     }
    // }

    // for(auto it = delete_seq.begin(); it != delete_seq.end(); it++)
    // {
    //     auto delete_pos = std::find(reparse_replacement.begin(), reparse_replacement.end(), *it);
    //     while(delete_pos != reparse_replacement.end())
    //     {
    //         reparse_replacement.erase(delete_pos);
    //         delete_pos = std::find(reparse_replacement.begin(), reparse_replacement.end(), *it);
    //     }
    // }

    // std::cerr << "Sequence after delete: ";
    // for(auto it=reparse_replacement.begin(); it != reparse_replacement.end(); it++)
    // {
    //     std::cerr <<"(" << (*it).GetContent() << ", " << (*it).GetId() << ") ";
    // }


    // std::cerr << "Print out the corresponding information before insertion: \n";
    // std::cerr << "Size of the unparsed positions: " << this->parseposition.unparsed_position.size() << "\n"
    // << "Print out where to do the insertions: \n";


    // for (int i = 0; i < this->parseposition.unparsed_position.size(); i++)
    // {
    //     // std::cerr << "Insertion round [" << i << "]:\n";
    //     // already parsed + unparsed before itself
    //     int already_parsed_start = this->parseposition.unparsed_position[i].second;
    //     int unparsed_position_in_seq = reparsed_seq.size() - this->parseposition.unparsed_position[i].first;

    //     // std::cerr << "The position to do the insertion: " << already_parsed_start + i << "\n";
    //     // std::cerr << "The position of the node to be inserted in the original sequence: " << unparsed_position_in_seq << "\n";
    // }

    //     for (int i = 0; i < this->parseposition.unparsed_position.size(); i++)
    // {
    //     // already parsed + unparsed before itself
    //     int already_parsed_start = this->parseposition.unparsed_position[i].second;
    //     int unparsed_position_in_seq = input_seq.size() - this->parseposition.unparsed_position[i].first;
    //     //std::cout << "Print out the inserted symbolic state. ID: " << input_seq[unparsed_position_in_seq].GetId() <<
    //     //"  Content: " << input_seq[unparsed_position_in_seq].GetContent() << "\n";
    //     states_under_dummy.insert(states_under_dummy.begin() + already_parsed_start + i, 
    //                                             input_seq[unparsed_position_in_seq]);
        
    // }

    // std::cerr << "Start inserting nodes: \n";
    for (int i = 0; i < this->parseposition.unparsed_position.size(); i++)
    {
    // already parsed + unparsed before itself
    int already_parsed_start = this->parseposition.unparsed_position[i].second;
    int unparsed_position_in_seq = reparsed_seq.size() - this->parseposition.unparsed_position[i].first;
    reparse_replacement.insert(reparse_replacement.begin() + already_parsed_start + i,
                            reparsed_seq[unparsed_position_in_seq]);
    }
    // std::cerr << "Out of ParsingInsertion\n";
    

    // std::cerr << "Print the whole reparsed sequence:\nContent: ";
    // for(auto && it=reparse_replacement.begin(); it != reparse_replacement.end(); it++)
    // {
    //     std::cerr << (*it).GetContent() << " ";
    // }
    // std::cerr << "\nID: ";
    // for(auto && it=reparse_replacement.begin(); it != reparse_replacement.end(); it++)
    // {
    //     std::cerr << (*it).GetId() << " ";
    // }

    return reparse_replacement;
}


template<class StateType>
bool Online_Learner<StateType>::Checkpoint:: Delete()
{
    std::cerr << "-------- Inside Delete --------\n";
    assert(this->AOF_cache_map_.empty());
    // this->graph_.OutputGraph("before_delete", PATH, false);
   
    std::unordered_map<std::vector<Symbolic_State<StateType> >, double> third_level = this->GetThirdLevel();
    // std::cerr << "Print data buffer: \n";
    // for(auto it = this->data_buffer_.begin(); it != this->data_buffer_.end(); it++)
    // {
    //     for(auto iter = (*it).begin(); iter != (*it).end(); iter++)
    //     {
    //         std::cerr << "(" << (*iter).GetContent() << ", " << (*iter).GetId() << ") ";
    //     }
    //     std::cerr << "\n";
    // }

    // this->PrintThirdLevel();

    //choose unique AOF roots that only appear in the third level
    std::unordered_set<Symbolic_State<StateType> > unique_aof_root;
    std::vector<Symbolic_State<StateType> > all_aof_root;
    for( auto && iter = third_level.begin(); iter != third_level.end(); ++iter)
    {
        for(Symbolic_State<StateType> state : iter->first)
        {
            //if the state is not a leaf state, treat it as an aof root
            if(state.GetId() != -1)
            {
                //if this aof root has not appeared before
                if(unique_aof_root.find(state) == unique_aof_root.end())
                {
                    unique_aof_root.insert(state);
                    
                    //if the aof root found was just added by its parent, skip it
                    if(state == this->AOFragment_.and_node_)
                        continue;

                    //if the state does not only exist in third level, skip it
                    //check all parents of the and node to see if they are either dummy or root
                    VertexId and_node_id = this->graph_.GetVertexIdByState(state);
                    std::vector<VertexId> and_node_parents_id = this->graph_.ParentsVertices(and_node_id);
                    
                    //throw exception if aof root has no parent
                    if(!and_node_parents_id.size())
                    {
                        std::cerr<<"the aof root does not have any parents!\n";
                        throw std::exception();
                    }
                    
                    //get the root's state id
                    VertexId root_vt_id = this->graph_.GetRoot();
                    int root_state_id = this->graph_.GetStateByVertexId(root_vt_id).GetId();

                    bool true_third_level = true;
                    //check if all parents has state id the same as root id
                    for(VertexId parent_vt_id: and_node_parents_id)
                    {
                        if(this->graph_.GetStateByVertexId(parent_vt_id).GetId() != root_state_id)
                        {
                            true_third_level = false;
                            break;
                        }
                    }

                    //if this aof only exists in third level, add it as a candidate
                    if(true_third_level)
                    {
                        all_aof_root.push_back(state);
                    }
                    
                }
                
            }
        }
    }

    int aof_root_size = all_aof_root.size();
    // std::cerr << "The size of the aof_roots: " << aof_root_size << "\n";
    if(!aof_root_size)
    {
        std::cerr<<"no aof root can be selected!\n";
        return false;
    }

    // std::cerr<<"all aof roots size: "<<aof_root_size<<std::endl;
    // for(auto root : all_aof_root)
    // {
    //     std::cerr<<"("<<root.GetContent()<<","<<root.GetId()<<")\n";
    // }
    //choose a aof to delete
    std::vector<int> uniform_weight(aof_root_size,1);
    std::discrete_distribution<int> uniform_dis(uniform_weight.begin(),uniform_weight.end());
    // Symbolic_State<std::string> AndUnderRoot(1, "Body");
    // Symbolic_State<StateType> chosen_aof_root = AndUnderRoot;
    Symbolic_State<StateType> chosen_aof_root = all_aof_root[uniform_dis(gen)];

    // Find the sequences with this symbolic state in the third level: (assumption: there will not be same sequences in the third level)
    std::vector<SequenceType<StateType>> sequences_with_target;
    for(auto && it = third_level.begin(); it != third_level.end(); it++)
    {
        SequenceType<StateType> third_level_sequence;
        bool target_found = false;
        for(auto && iter = it->first.begin(); iter != it->first.end(); iter++)
        {   
            third_level_sequence.push_back(*iter);
            if((*iter) == chosen_aof_root)
                target_found = true;
        }

        if(target_found)
            sequences_with_target.push_back(third_level_sequence);
    }

    // Symbolic_State<std::string> chosen_aof_root(2, "dummy");
    AOFStruct<StateType> chosen_aof = FindAOFFromGraph(chosen_aof_root);

    // for(auto or_node : chosen_aof.or_children_)
    // {
    //     for(auto state : or_node)
    //     {
    //         if(state.GetId() == 4561)
    //         {
    //             std::cerr<<"find 4561\n";
    //             int a;
    //             std::cin >>a;
    //         }
    //     }
    // }
    // if(chosen_aof_root.GetId() == 4561)
    // {
    //     std::cerr << "Debug node 4561\n";
    //     for(auto or_node : this->AOFragment_.or_children_)
    //     {
    //         for(auto state : or_node)
    //             std::cerr<< "("<<state.GetContent()<<","<< state.GetId()<<")";
    //         std::cerr<<std::endl;
    //     }
    //     std::cerr<<"find 4561 as chosen AOF root to be deleted\n";
    //     int a ;
    //     std::cin >> a;
    // }
    // if(chosen_aof.and_node_.GetId() == 4561)
    // {
    //     std::cerr << "Debug node 4561\n";
    //     for(auto or_node : this->AOFragment_.or_children_)
    //     {
    //         for(auto state : or_node)
    //             std::cerr<< "("<<state.GetContent()<<","<< state.GetId()<<")";
    //         std::cerr<<std::endl;
    //     }

    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<< "("<<state.GetContent()<<","<< state.GetId()<<")";

    //         }
    //         std::cerr<<std::endl;
    //     }
        
    //     this->graph_.OutputGraph("debug_delete", PATH, true);
    //     std::cerr<<"find 4561 as and-node\n";
    //     int a;
    //     std::cin >> a;
    // }

    std::cerr<<"the chosen AOF to be deleted:\n";
    std::cerr<<"and_node: "<<"("<<chosen_aof.and_node_.GetContent()<<", "<<chosen_aof.and_node_.GetId()<<")\n" << "VertexID: " << this->graph_.GetVertexIdByState(chosen_aof_root) << "\n";
    std::cerr<<"Is the And-node a leaf? :"<<chosen_aof.and_node_.GetIsBasic()<<std::endl;

    for(const auto & or_node : chosen_aof.or_children_)
    {
        for(const auto & state : or_node)
        {
            std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
        }
        std::cerr<<std::endl;
    }
    // this->graph_.OutputGraph("before_delete",PATH,true,SIMPLIFY);
    // std::cerr<<"the corresponding rules:\n (size: " << chosen_aof.all_rules_.size() << ")\n";
    // for(auto rule : chosen_aof.all_rules_)
    // {
    //     AOG_LIB::Symbolic_State<StateType> source = rule.GetSource();
    //     std::vector<AOG_LIB::Symbolic_State<StateType> > results = rule.GetResults();
    //     for (AOG_LIB::Symbolic_State<StateType> result : results)
    //     {
    //         std::cerr << "\t" << "(" << source.GetContent() << "," << source.GetId() << ")" 
    //             << " -> " << "(" << result.GetContent() << "," << result.GetId() << ")" << "\n";
    //     }
    // }
   
    //find rrmaps used to update graph and buffer
    this->AOFragment_ = chosen_aof;
    // this-> graph_.OutputGraph("visuaize_meaningless",PATH,true);
    // this-> graph_.OutputGraph("visuaize_meaningless_no_trunc",PATH,false);
    // std::vector<VertexId> chosen_aof_parents_id = this->graph_.ParentsVertices(this->graph_.GetVertexIdByState(chosen_aof_root));
    // std::cerr << "The parents of the chosen AOF: \n";
    // for(int i=0; i < chosen_aof_parents_id.size(); i++)
    // {
    //     std::cerr << "state ID: " << this->graph_.GetStateByVertexId(chosen_aof_parents_id[i]).GetId() << " VertexID: " << chosen_aof_parents_id[i] << "\n";
    // }
    // std::cerr << "\n";

    // this->graph_.OutputGraph("before_delete",PATH,true,SIMPLIFY);
    RRMap<StateType> rrmap = FindRRMap(chosen_aof);
    Replacement<StateType> delete_replacement; 
    delete_replacement.rrmap_ = rrmap;
    
    // if(chosen_aof.and_node_.GetId() == 44007)
    // {
    //     for(auto iter = rrmap.map_.begin(); iter != rrmap.map_.end();++iter)
    //     {
    //         for(auto state : iter -> first)
    //         {
    //                 std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<< "->\n";

    //         for(auto seq : iter -> second)
    //         {
    //             for(auto state : seq)
    //             {
    //                 std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //             }
    //             std::cerr<<std::endl;
    //         }
    //     }
    //     std::cerr<<std::endl;
        
    //     int a;
    //     std :: cin >> a;
    // }

    // for(auto iter = rrmap.map_.begin(); iter != rrmap.map_.end();++iter)
    // {
    //     for(auto state : iter -> first)
    //     {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //     }
    //     std::cerr<< "->\n";

    //     for(auto seq : iter -> second)
    //     {
    //         for(auto state : seq)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    // }
    //update graph and buffer
    UpdateGraphAndBuffer(delete_replacement);
    
    // if(chosen_aof.and_node_.GetId() == 59107)
    // {
    //     for(auto rule : this->graph_.GetRules())
    //     {
    //         std::cerr<<"("<<rule.GetSource().GetContent()<<", "<<rule.GetSource().GetId()<<")->\n";
    //         for(auto state : rule.GetResults())
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }

    //     // auto parents = this->graph_.ParentsVertices(406);
    //     auto children = this->graph_.ChildrenVertices(685);
    //     // std::cerr<<"406's parent vertices: \n";
    //     // for(auto parent : parents)
    //     //     std::cerr << parent<<" ";
    //     // std::cerr <<std::endl;

    //     // std::cerr<<"685's children vertices: \n";                
    //     // for(auto child : children)
    //     //     std::cerr<<child<<" ";
    //     // std::cerr<<std::endl;

    //     // std::cerr << "find you in delete()"<<std::endl;
    //     // this->graph_.OutputGraph("delete_44007_after",PATH,true,SIMPLIFY);

    //     // int a;
    //     // std:: cin >> a;
    // }


    //handle possible nodes that are not sampled in AOG
    for(const auto & or_node : this->AOFragment_.or_children_)
    {
        if(or_node.size() > 1)
        {
            for(const auto & state : or_node)
                this->graph_.DeleteNoParentRules(state);
        }
    }

    // if(chosen_aof.and_node_.GetId() == 937)
    // {
    //     std::cerr<<"the data buffer: \n";
    //     for(auto data : this-> data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     std::cerr << "find you in delete()"<<std::endl;
    //     this->graph_.OutputGraph("middle_delete_937",PATH,true);
    //     // int a;
    //     // std:: cin >> a;
    // }

    // Now suppose we already found the possible replacement for all the sequences in the third level:


    
    // if(chosen_aof.and_node_.GetId() == 24)
    // {
    //     std::cerr<<"graph after deleting 24\n";
    //     this->graph_.OutputGraph("visualize_delete_24",PATH,true);
    //     // int a;
    //     // std::cin >>a;
    // }

    ConfigsMap<StateType> config_map = GenerateKConfigs(chosen_aof.or_children_.size());
    ConfigBuffer<StateType> all_configs = chosen_aof.GenerateAllConfigsAndWeights(config_map);
    RdBuffer<StateType> reductions = GenerateCMAndReductions(all_configs,this->data_buffer_);
    
    // if(reductions.empty())
    // {
        // std::cerr<<"the chosen AOF:\n";
        // for(const auto & or_node : chosen_aof.or_children_)
        // {
        //     for(const auto & state : or_node)
        //     {
        //         std::cerr<< "("<<state.GetContent()<<","<<state.GetId()<<")";
        //     }
        //     std::cerr<<std::endl;
        // }

        // std::cerr<<"the data buffer:\n";
        // for(const auto & data : this->data_buffer_)
        // {
        //     for(const auto & state : data)
        //     {
        //         std::cerr<< "("<<state.GetContent()<<","<<state.GetId()<<")";
        //     }
        //     std::cerr<<std::endl;
        // }

        // std::cerr<<"the third level:\n";
        // std::unordered_map<std::vector<Symbolic_State<StateType> >, double> third_lvl = this->GetThirdLevel();
        // for(auto && iter = third_lvl.begin();iter != third_lvl.end();++iter)
        // {
        //     for(const auto & state : iter->first)
        //     {
        //         std::cerr<< "("<<state.GetContent()<<","<<state.GetId()<<")";

        //     }
        // }

        // int a;
        // std::cin >> a;
        // int b;
        // std::cin >> b;
    // }

    double post_gain = CalculatePosterior(chosen_aof,reductions,false,true);
    this->posterior_ -= post_gain;  
    std::cerr << "posterior after delete AOF: " << this->posterior_<<std::endl;

    boost::multi_index::get<3>(this->CM_).clear();

    // if(this->AOFragment_.and_node_.GetId() == 18370)
    // {
    //     std::cerr<<"all rules about to use to parse 18370:\n";
    //     auto all_rules = this->graph_.GetRules();
    //     for(auto rule : all_rules)
    //     {
    //         std::cerr<<"("<<rule.GetSource().GetContent()<<", "<< rule.GetSource().GetId()<<")->";
    //         for(auto target : rule.GetResults())
    //         {
    //             std::cerr<<"("<<target.GetContent()<<", "<<target.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     std::cerr<<"before reparse in delete 18370!\n";
    //     this->graph_.OutputGraph("before_reparse",PATH,true,SIMPLIFY);
    //     int a;
    //     std::cin >>a;
    // }


    // Find possible replacements for the sequences with the target node:
    std::vector<SequenceType<StateType>> possible_third_level_sequences;
    for(auto && it = rrmap.map_.begin(); it != rrmap.map_.end(); it++)
    {
        for(auto && iter = it->second.begin(); iter != it->second.end(); iter++)
        {
            //if(std::find(possible_third_level_sequences.begin(), possible_third_level_sequences.end(), *iter) == possible_third_level_sequences.end())
            possible_third_level_sequences.push_back(*iter);
        }
    }

    //Get new third level and find all the sequences to be reparsed:
    std::unordered_map<std::vector<Symbolic_State<StateType> >, double> new_third_level = this->GetThirdLevel();
    std::vector<SequenceType<StateType>> sequences_to_be_reparsed;
    for(auto && it = new_third_level.begin(); it != new_third_level.end(); it++)
    {
        if(std::find(possible_third_level_sequences.begin(), possible_third_level_sequences.end(), it->first) != possible_third_level_sequences.end())
            sequences_to_be_reparsed.push_back(it->first);
    }

    for(auto && it = sequences_to_be_reparsed.begin(); it != sequences_to_be_reparsed.end(); it++)
    {
        // std::cerr << "\n--------- Print third level before updating databuffer and graph---------\n";
        // this->PrintThirdLevel();
        // std::cerr << "\n-------------------------------------------------------------------\n";

        double old_prior = -ALPHA_PRIOR * this->graph_.GrammarSize(1);
        double old_posterior = this->posterior_;

        double reduced_likelihood = 1;
        //TODO: check whether reduced_likelihood would be returned as -1
        SequenceType<StateType> reparsed_sequence = ParseDelete(*it, reduced_likelihood);

        // If the reparsed sequence is the same as before, do nothing, o.w do the following:
        if(reparsed_sequence != *it)
        {
            double increased_likelihood = 1;
            int count_data_buffer_old = 0;
            int count_data_buffer_new = 0;
            for(auto && iter = this->data_buffer_.begin(); iter != this->data_buffer_.end(); iter++)
            {
                if(*iter == reparsed_sequence)
                    count_data_buffer_old++;
                
                if(*iter == *it)
                    count_data_buffer_new++;
            }

            if(count_data_buffer_new != 0 && count_data_buffer_old != 0)
                increased_likelihood = (count_data_buffer_new + count_data_buffer_old) / count_data_buffer_old;

            double likelihood_gain = ALPHA_LIKELIHOOD * (log(increased_likelihood) + log(reduced_likelihood));
            std::cerr << "The reduced likelihood is: " << reduced_likelihood << "\n";
            std::cerr << "The increased likelihood is: " << increased_likelihood << "\n";
            std::cerr << "The likelihood gain is: " << likelihood_gain << "\n";

            //Delete in the graph
            // if there exist the sequence in the graph
            std::vector<VertexId> dummy_children = this->graph_.ChildrenVertices(this->graph_.GetRoot());
            std::unordered_map<VertexId, double> weights = this->graph_.GetOutEdgeWeights(this->graph_.GetRoot(), false);

            bool reparsed_sequence_found = false;
            int reparsed_sequence_weight = 0;
            int sequence_to_be_reparsed_weight = 0;
            VertexId reparsed_sequence_vertex;
            VertexId sequence_to_be_reparsed_vertex;

            for (int i = 0; i < dummy_children.size(); i++)
            {   
                // record third level sequence == reparsed_sequence, also third level sequence == sequence_to_be_reparsed
                std::vector<VertexId> third_level_children = this->graph_.ChildrenVertices(dummy_children[i]);
                SequenceType<StateType> third_level_child;
                for (VertexId child : third_level_children)
                    third_level_child.push_back(this->graph_.GetStateByVertexId(child));

                if(third_level_child == reparsed_sequence)
                {
                    reparsed_sequence_found = true;
                    reparsed_sequence_weight = weights[dummy_children[i]];
                    reparsed_sequence_vertex = i;
                }

                if(third_level_child == *it)
                {
                    sequence_to_be_reparsed_weight = weights[dummy_children[i]];
                    sequence_to_be_reparsed_vertex = i;
                }
            }


            //if found, add weight to it and delete the sequence_to_be_reparsed
            if(reparsed_sequence_found)
            {  
                // std::cerr << "Sequence found in third level!\n";
                // Delete sequence_to_be_reparsed
                std::vector<Symbolic_State<std::string>> under_root = *it;
                Symbolic_Rule<std::string> rule_to_data(this->graph_.GetRoot(), under_root);
                // if(rule_to_data.GetResults()[0].GetContent() == "a21")
                // {
                //     std::cerr<<"the rule: \n"<< "("<<rule_to_data.GetSource().GetContent()<<", "<<rule_to_data.GetSource().GetId()<<") ->";
                //     for(auto target : rule_to_data.GetResults())
                //     {
                //         std::cerr<<"("<<target.GetContent()<<", "<<target.GetId()<<")";
                //     }
                //     std::cerr<<std::endl;
                //     std::cerr<<"find you in reparse found!\n";
                //     std::cerr<<"the aof to be deleted: \n";
                //     for(auto or_node : this->AOFragment_.or_children_)
                //     {
                //         for(auto state : or_node)
                //         {
                //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
                //         }
                //         std::cerr<<std::endl;
                //     }
                //     int a;
                //     std::cin >>a;
                // }
               
                this->graph_.DeleteRule(rule_to_data);
                
                weights = this->graph_.GetOutEdgeWeights(this->graph_.GetRoot(), false);
                weights[dummy_children[reparsed_sequence_vertex]] = reparsed_sequence_weight + sequence_to_be_reparsed_weight;
                // std::cerr<<"the updated weight in delete: "<< weights[dummy_children[reparsed_sequence_vertex]]<<std::endl;
                // std::cerr<<"all the weights to be updated: \n";
                // for(auto iter = weights.begin(); iter != weights.end(); ++iter)
                // {
                //     std::cerr<<iter->second<<" ";
                //     if(iter->second == 0)
                //     {
                //         std::cerr<<"update weight in delete if branch error!\n";
                //         throw std::exception();
                //     }
                // }
                // std::cerr<<std::endl;
                this->graph_.SetOutEdgeWeights(this->graph_.GetRoot(), weights);
            }
            else  // add a new rule and set its weight
            {
                //Delete the old rule
                std::vector<Symbolic_State<std::string>> under_root = *it;
                Symbolic_Rule<std::string> rule_to_data(this->graph_.GetRoot(), under_root);
                this->graph_.DeleteRule(rule_to_data);
                //  if(rule_to_data.GetResults()[0].GetContent() == "a21")
                // {
                //     std::cerr<<"the rule: \n"<< "("<<rule_to_data.GetSource().GetContent()<<", "<<rule_to_data.GetSource().GetId()<<") ->";
                //     for(auto target : rule_to_data.GetResults())
                //     {
                //         std::cerr<<"("<<target.GetContent()<<", "<<target.GetId()<<")";
                //     }
                //     std::cerr<<std::endl;
                //     std::cerr<<"find you in reparse not found!\n";
                //     std::cerr<<"the aof to be deleted: \n";
                //     for(auto or_node : this->AOFragment_.or_children_)
                //     {
                //         for(auto state : or_node)
                //         {
                //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
                //         }
                //         std::cerr<<std::endl;
                //     }
                //     int a;
                //     std::cin >>a;
                // }

                // std::cerr << "Sequence not found in third level, so add the reparsed sequence to the graph\n";
                Symbolic_Rule<StateType> root_to_new_parsing(this->graph_.GetRoot(), reparsed_sequence);
                this->graph_.AddRule(root_to_new_parsing);

                // set its weight
                std::vector<VertexId> dummy_children = this->graph_.ChildrenVertices(this->graph_.GetRoot());
                std::unordered_map<VertexId, double> weights = this->graph_.GetOutEdgeWeights(this->graph_.GetRoot(), false);

                for (int i = 0; i < dummy_children.size(); i++)
                {
                    std::vector<VertexId> third_level_children = this->graph_.ChildrenVertices(dummy_children[i]);
                    SequenceType<StateType> third_level_child;
                    for (VertexId child : third_level_children)
                        third_level_child.push_back(this->graph_.GetStateByVertexId(child));

                    if(third_level_child == reparsed_sequence)
                    {
                        weights[dummy_children[i]] = sequence_to_be_reparsed_weight;

                        // std::cerr<<"the updated weight in delete: "<< weights[dummy_children[i]]<<std::endl;
                        // std::cerr<<"all the weights to be updated: \n";
                        // for(auto iter = weights.begin(); iter != weights.end(); ++iter)
                        // {
                        //     std::cerr<<iter->second<<" ";
                        //     if(iter->second == 0)
                        //     {
                        //         std::cerr<<"update weight in delete else branch error!\n";
                        //         throw std::exception();
                        //     }
                        // }
                        // std::cerr<<std::endl;
                        this->graph_.SetOutEdgeWeights(this->graph_.GetRoot(), weights);
                    }
                }
                // this->graph_.OutputGraph("after_reparse",PATH,true);
            }

            // std::cerr << "Print the data buffer: \n";
            // for(auto it = this->data_buffer_.begin(); it != this->data_buffer_.end(); it++)
            // {
            //     for(auto iter = (*it).begin(); iter != (*it).end(); iter++)
            //     {
            //         std::cerr << "(" << (*iter).GetContent() << ", " << (*iter).GetId() << ") ";
            //     }
            //     std::cerr << "\n";
            // }


            // update databuffer:
            for(auto iter = this->data_buffer_.begin(); iter != this->data_buffer_.end(); iter++)
            {
                if((*iter) == (*it))
                    (*iter) = reparsed_sequence;
            }

            // std::cerr << "Print the whole reparsed sequence:\nContent: ";
            // for(auto iter=reparsed_sequence.begin(); iter != reparsed_sequence.end(); iter++)
            // {
            //     std::cerr << (*iter).GetContent() << " ";
            // }
            // std::cerr << "\nID: ";
            // for(auto iter=reparsed_sequence.begin(); iter != reparsed_sequence.end(); iter++)
            // {
            //     std::cerr << (*iter).GetId() << " ";
            // }

            // std::cerr << "\nPrint the sequence to be replaced: \n";
            // for(auto iter=(*it).begin(); iter != (*it).end(); iter++)
            // {
            //     std::cerr << (*iter).GetContent() << " ";
            // }
            // std::cerr << "\nID: ";
            // for(auto iter=(*it).begin(); iter != (*it).end(); iter++)
            // {
            //     std::cerr << (*iter).GetId() << " ";
            // }


            // Delete in the data buffer
            // auto match_pos = std::find(this->data_buffer_.begin(), this->data_buffer_.end(), *it);
            // while(match_pos != this->data_buffer_.end())
            // {
            //     std::cerr << "sequence found in data buffer\n";
            //     this->data_buffer_[match_pos - this->data_buffer_.begin()] = reparsed_sequence;
            //     auto match_pos = std::find(this->data_buffer_.begin(), this->data_buffer_.end(), *it);
            // }

            // std::cerr << "Landmark\n";



            // update the posterior after reparsing:
            double new_prior = -ALPHA_PRIOR * this->graph_.GrammarSize(1);
            double old_likelihood = old_posterior - old_prior;
            double new_likelihood = old_likelihood + likelihood_gain;
            std::cerr << "Old posterior: " << this->posterior_ << "\n";
            this->posterior_ = new_likelihood + new_prior; //new posterior
            std::cerr << "New posterior: " << this->posterior_ << "\n";

            // std::cerr << "\n--------- Print third level after updating databuffer and graph ---------\n";
            // this->PrintThirdLevel();
            // std::cerr << "\n-------------------------------------------------------------------\n";

        }

    }
    std::cerr << "posterior after delete: " << this->posterior_<<std::endl;
    // std::cerr << "-------- Out of Delete --------\n";
    // if(chosen_aof.and_node_.GetId() == 937)
    // {
    //     std::cerr << "find you in delete()"<<std::endl;
    //     this->graph_.OutputGraph("after_delete_937",PATH,true);
    //     // int a;
    //     // std:: cin >> a;
    // }
    // Symbolic_State<std::string> a31("a31",true);
    // bool print_buffer = false;
    // for(auto data : this->data_buffer_)
    // {
    //     for(int i = 0; i < data.size(); ++i)
    //     {
    //         if(data[i] == a31 && i < data.size()-1)
    //         {
    //             std::cerr<<"find the reverse order !\n";
    //             print_buffer = true;
    //             break;
    //         }
    //     }
    // }
    // if(print_buffer)
    // {
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     int a;
    //     std::cin >>a;
    // }
    // if(chosen_aof.and_node_.GetId() == 73981)
    // {
    //     std::cerr<<"the data buffer: \n";
    //     for(auto data : this-> data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }

    //     for(auto rule : this->graph_.GetRules())
    //     {
    //         std::cerr<<"("<<rule.GetSource().GetContent()<<", "<<rule.GetSource().GetId()<<")->\n";
    //         for(auto state : rule.GetResults())
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     std::cerr << "find you in delete()"<<std::endl;
    //     this->graph_.OutputGraph("after_delete_73981",PATH,true);
    //     int a;
    //     std:: cin >> a;
    // }
    // bool wrong_buffer = false;
    // for(auto data : this->data_buffer_)
    // {
    //     if(data[data.size()-1].GetId() == 261)
    //     {
    //         std::cerr<<"wrong data buffer in delete\n";
    //         wrong_buffer = true;
    //         break;
    //     }
    // }
    // if(wrong_buffer)
    // {
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     this->graph_.OutputGraph("wrong_buffer",PATH,true,false);
    //     int a;
    //     std::cin >>a;
    // }
    if(NEG_GAIN_COUNT != -1)
    {
        if(this->posterior_ > this->best_posterior_)
        {
            this->best_posterior_ = this->posterior_;
            this->best_aog_ = T_AOG<StateType>(this->graph_);
            if(!OFFLINE)
                this -> best_buffer_ = this -> data_buffer_;
        }
    }
    
    std::cerr<<"-----------------------Out of Delete()-----------------------\n";
    return true;


}

template<class StateType>
bool Online_Learner<StateType>::Checkpoint::Generate(int iterations, bool use_stochastic)
{
    //now assume use same mechanism to find AOF in exoand() and simulate()
    // std::cerr<<"----------------Inside Generate()-------------------------\n"<<std::endl;
    // int c;
    // std::cin >> c;
    //get all configurations and configs map of size 2 from data buffer
    // std::random_device rd;
    // std::mt19937 gen(417);
    // NUM_OF_VARIATION_EACH_ROUND = 0;
    auto ori_CM = this->CM_;
    ConfigsMap <StateType> one_configs_map = GenerateKConfigs(1);
    ConfigsMap<StateType> two_configs_map = GenerateKConfigs(2);

    // used to generate statistics about num of bigram in each call of Generate()
    int num_of_two_configs= two_configs_map.size();


    std::vector<SequenceType<StateType> > one_frags;
    std::vector<SequenceType<StateType> > two_frags;
    

    for(const auto & config: one_configs_map)
        one_frags.push_back(config.first);
        
    /** std::cerr<<"Just inside Generate()\n";
    std::cerr<<"size of two_configs_map: "<<two_configs_map.size()<<std::endl; */
    for (const auto & config : two_configs_map)
        two_frags.push_back(config.first);
    

    //generate frags and config maps of size 3 to be used in AddOr()
    std::vector<SequenceType<StateType> > three_frags;
    ConfigsMap<StateType> three_configs_map =  GenerateKConfigs(3);
    
    
    for (const auto & config : three_configs_map)
        three_frags.push_back(config.first);


    //store configurations and configs_map of known size k for search
    std::vector<ConfigsMap<StateType> > all_k_maps = {one_configs_map,two_configs_map,three_configs_map};
    std::vector<ConfigBuffer<StateType> >  all_k_frags = {one_frags,two_frags, three_frags};

    //std::cerr<<"size of counts: "<<counts.size()<<std::endl;

    std::vector<AOFStruct<StateType> > best_AOF_of_all_rds;
    std::vector<double> best_lg_of_all_rds;
    std::vector<RdBuffer<StateType> > best_reductions_of_all_rds;
    int num_of_iter = iterations;
    int test = 0;
    while(num_of_iter > 0 )
    {    
        std::vector<SequenceType<StateType> > curr_frag = all_k_frags[1];
        std::vector<SequenceType<StateType> > next_frag = all_k_frags[2];
        ConfigsMap<StateType> prev_map = all_k_maps[0];
        ConfigsMap<StateType> curr_map = all_k_maps[1];
        ConfigsMap<StateType> next_map = all_k_maps[2];
        
        RdBuffer<StateType> best_reductions;
        //choose the initial bigram, weights are uninitialized
        
        AOFStruct<StateType> AOF = ChooseBigram();
        // std::cerr<<"current round: "<<num_of_iter<<std::endl;
        assert(AOF.all_rules_.size() == 0);
        // std::cerr<<"the chosen initial bigram: \n";
        // bool changed = 0;

        // for(auto or_node : AOF.or_children_)
        // {
        //     for(auto state : or_node)
        //     {
        //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
        //         // if(state.GetId() == 16)
        //         // {
        //         //     changed = 1;
        //         // }
        //     }
        //     std::cerr<<std::endl;
        // }
        // if(changed == 1)
        // {
        //     std::vector<Symbolic_State<std::string>> test_or_1 = {Symbolic_State<std::string>("Round", true)};
        //     std::vector<Symbolic_State<std::string>> test_or_2 = {Symbolic_State<std::string>("Square", true), Symbolic_State<std::string>( "Round",true)};
        //     AOFStruct<std::string> test_bigram;
        //     test_bigram.or_children_.push_back(test_or_1);
        //     test_bigram.or_children_.push_back(test_or_2);
        //     AOF.or_children_ = test_bigram.or_children_;
        //     AOF.weights_[0] = std::vector<double>(AOF.or_children_[0].size(), 0);
        //     AOF.weights_[1] = std::vector<double>(AOF.or_children_[1].size(), 0);

        // }
       
        // if(AOF == test_bigram)
        // {
        //     std::cerr<<"about to change AOF\n";
        //     int a;
        //     std::cin >>a;
        //     changed = 1;
        //     AOF.or_children_[0].push_back(Symbolic_State<std::string>("Square", true));
        // }
        // std::cerr<<std::endl;
        // //test AOF
        // AOFStruct<StateType> test_aof;
        // std::vector<Symbolic_State<std::string> > test_or_child_1;        
        // std::vector<Symbolic_State<std::string> > test_or_child_2;
        // // std::vector<Symbolic_State<std::string> > test_or_child_4;
        
        // test_or_child_1.push_back(Symbolic_State<std::string> ("b1B2",true));       
        // test_or_child_2.push_back(Symbolic_State<std::string> ("b2A1",true));

        // // test_or_child_2.push_back(Symbolic_State<std::string>(child1));
        // // test_or_child_2.push_back(Symbolic_State<std::string>(child2));

        // // test_or_child_4.push_back(Symbolic_State<std::string>("a31",true));
        // // test_or_child_4.push_back(Symbolic_State<std::string>("a32",true));
        
        // test_aof.or_children_.push_back(test_or_child_1);       
        // test_aof.or_children_.push_back(test_or_child_2);
        // // test_aof.or_children_.push_back(test_or_child_3);
        // // test_aof.or_children_.push_back(test_or_child_4);
        
        // test_aof.weights_.push_back(std::vector<double>(1,0)) ;
        // test_aof.weights_.push_back(std::vector<double>(1,0)) ;
        // // test_aof.weights_.push_back(std::vector<int>(1,0)) ;        
        // // test_aof.weights_.push_back(std::vector<int>(2,0)) ;
        // // test_aof.weights_.push_back(std::vector<int>(2,0)) ;
        

        // // //ground truth AOF
        // // std::vector<Symbolic_State<std::string> > truth_or_child_1;
        // // std::vector<Symbolic_State<std::string> > truth_or_child_2;
        // // std::vector<Symbolic_State<std::string> > truth_or_child_3;
        // // std::vector<Symbolic_State<std::string> > truth_or_child_4;

        // // truth_or_child_1.push_back(Symbolic_State<std::string> ("a21",true));
        // // truth_or_child_1.push_back(Symbolic_State<std::string> ("a22",true));
        // // truth_or_child_2.push_back(Symbolic_State<std::string> ("a31",true));
        // // truth_or_child_2.push_back(Symbolic_State<std::string> ("a32",true));

        // // // truth_or_child_3.push_back(Symbolic_State<std::string> ("a1A1",true));
        // // // truth_or_child_4.push_back(Symbolic_State<std::string> ("a1A2",true));
        
        
        // // AOFStruct<std::string> truth_aof;
        // // truth_aof.or_children_.push_back(truth_or_child_1);
        // // truth_aof.or_children_.push_back(truth_or_child_2);
        // // truth_aof.or_children_.push_back(truth_or_child_3);
        // // truth_aof.or_children_.push_back(truth_or_child_4);
        
        // // truth_aof.and_node_ = Symbolic_State<std::string>(43);
        // // truth_aof.weights_.push_back(std::vector<int>(2,0));
        // // truth_aof.weights_.push_back(std::vector<int>(2,0));
        // // truth_aof.weights_.push_back(std::vector<int>(1,0));
        // // truth_aof.weights_.push_back(std::vector<int>(1,0));
        
        // // ConfigsMap<StateType> four_configs_map = GenerateKConfigs(4);
        
        // ConfigBuffer<std::string> test_configs = test_aof.GenerateAllConfigsAndWeights(two_configs_map);
        // // ConfigBuffer<std::string> truth_configs = truth_aof.GenerateAllConfigsAndWeights(four_configs_map);

        // GenerateCM(test_configs);

        // RdBuffer<std::string> test_reductions = FindReduction(test_configs);
        // double test_lg = CalculatePosterior(test_aof,test_reductions);
        // std::cerr << "test_aof posterior: " << test_lg << std::endl;
        // std::cerr << "------------------------------" << std::endl;

        // GenerateCM(truth_configs);
        // RdBuffer<std::string> truth_reductions = FindReduction(truth_configs);

        // std::cerr<<"truth configs: \n";
        // for(auto config : truth_configs)
        // {
        //     for(auto state : config)
        //         std::cerr<<state.GetContent()<<"_";
        //     std::cerr<<std::endl;
        // }

        // std::cerr<<"truth reductions: \n";
        // if(truth_reductions.empty())
        //     std::cerr<<"reduction is empty!\n";
        // for(auto rd : truth_reductions)
        // {
        //     for(auto state : rd.first)
        //     {
        //         std::cerr<<state.GetContent()<<"_";
                
        //     }
        //     std::cerr<<std::endl;
        // }
        // double truth_lg = CalculatePosterior(truth_aof,truth_reductions);

        


        //if cannot generate a bigram from data anymore, directly return
        if(AOF.or_children_.empty())
        {
            //print out stats
            // NUM_OF_BIGRAM.push_back(0);
            // NUM_OF_VARIATION.push_back(0);
            // if(this->posterior_ > BEST_POSTERIOR.back())
            // {
            //     BEST_POSTERIOR.push_back(this->posterior_);
            //     KL_DIVERGENCE.push_back(true);
            // }            
            // else
            // {
            //     BEST_POSTERIOR.push_back(BEST_POSTERIOR.back());
            //     KL_DIVERGENCE.push_back(false);

            // }
            // OutputToFile("../AOG_Sample_Results/sample_"+std::to_string(BEST_POSTERIOR.size()-1),this->graph_,this->graph_.GetRoot(),SAMPLING);
            std::cerr<<"return from empty bigram!\n";
            this->AOF_cache_map_.clear();
            return false;
        }
            
        // std::cerr<<std::endl;
        // std::cerr<<"Bigram summary:\n";
        
        /** std::cerr<<"weight in bigram: \n"; */
        
        // bool fstExt = false;
        // bool sndExt = false;
        // for(SequenceType<StateType> it : AOF.or_children_){
        //     for(Symbolic_State<StateType> iit : it){
        //         if(iit.GetContent() == "1"){
        //             fstExt = true;
        //         }
        //         else if (iit.GetContent() == "2"){
        //             sndExt = true;
        //         }
        //     }
        // }
        // bool exist = fstExt && sndExt;
        // if(exist){
        //     std::cerr<<"attention:!!!!!!!!!found!!!!!!!!!!! \n";
        //     for(int j = 0; j <  AOF.or_children_.size();++j)
        //     {
        //         std::cerr<<"or_node "<<j<<": ";
        //         for(auto child : AOF.or_children_[j])
        //         {
                    
        //             std::cerr<<child.GetContent()<<"_";
                    
        //         }
                    
        //         std::cerr<<std::endl;
        //     }
        // }
        
        //find all configurations represented by the bigram and update bigram's weights
        ConfigBuffer<StateType> configs = AOF.GenerateAllConfigsAndWeights(two_configs_map);
        RdBuffer<StateType> reductions;

        if(CACHE)
        {
             auto iter = this->AOF_cache_map_.find(AOF.or_children_);
   
            //if current AOF's corresponding CM and reductions have been cached before, use them directly
            if(iter != this->AOF_cache_map_.end())
            {
                // std::cerr << "cache find!\n";
                // int a;
                // std::cin >> a;
                Cache<StateType> cache = this -> AOF_cache_map_[AOF.or_children_];
                reductions = cache.rd_cache_;
                this->CM_ = cache.cm_cache_;
            }
            else
            {
                reductions =  GenerateCMAndReductions(configs,this->data_buffer_);
                
                //cache the generated reductions and CM
                Cache<StateType> cache = {this->CM_, reductions};
                this->AOF_cache_map_[AOF.or_children_] = cache;
            }
        }
        else
            reductions =  GenerateCMAndReductions(configs,this->data_buffer_);

       
        
        //calculate the likelihood of the bigram as the baseline
        AOFStruct<StateType> best_AOF = AOF;
        // std::cerr<<"bigram likelihood:\n";
        double best_lg = CalculatePosterior(best_AOF,reductions);
        double prev_best_lg = std::numeric_limits<double>::lowest();

        // std::cout<<"begin four operation=============\n";
        int num_of_or_node;
        double best_del_leaf_lg, best_add_leaf_lg, best_del_or_lg, best_add_or_lg;
        best_reductions = reductions;
        // RdBuffer<StateType> temp_rd = best_reductions;
        RdBuffer<StateType> best_del_leaf_rd, best_add_leaf_rd, best_del_or_rd,best_add_or_rd;
        //while likelihood keeps increasing
        while (best_lg > prev_best_lg)
        {
            // std::cerr<<"current best_lg: "<<best_lg<<std::endl;
            // std::cerr << "0:" << std::endl;
            // std::cerr<<"Start 4 operations\n";
            // std::cerr<<"Current AOF Related info:\n";
            // // std::cerr<<"num of or nodes: "<<AOF.or_children_.size()<<std::endl;
            // // for(int i = 0; i < AOF.or_children_.size();++i)
            // // {
            // //     std::cerr<<"or_node "<<i<<" size: "<<AOF.or_children_[i].size()<<std::endl;
                
            // // }
            // std::cerr<<"Each or_node's children: \n";
            // for(int j = 0; j <  AOF.or_children_.size();++j)
            // {
            //     for(const auto & child : AOF.or_children_[j])
            //     {
            //         std::cerr<<"("<<child.GetContent()<<", "<<child.GetId()<<")";
                    
            //     }
                    
            //     std::cerr<<std::endl;
            // } 
            
            /** std::cerr<<"Configurations generated from the AOF:\n";
            for(auto config : configs)
            {
                for(auto state : config)
                   std::cerr<<state.GetContent()<<"_";
                std::cerr<<std::endl;
            } */
                
            //print configurations
            /** std::cerr<<"size of configs: "<<configs.size()<<std::endl; */
            prev_best_lg = best_lg;
            //set 4 opeartions likelihood
            best_del_leaf_lg = best_add_leaf_lg = best_del_or_lg = best_add_or_lg = best_lg;
            best_del_leaf_rd = best_add_leaf_rd = best_del_or_rd = best_add_or_rd = best_reductions;
            // std::cout << test<<std::endl;
            // ++test;
            //calculate best AOF obtained by delete leaf node
            // bool target_found = false;
            // bool badaof_found = false;
            AOFStruct<StateType> best_aof_del_leaf = DeleteLeaf(AOF, best_del_leaf_rd,curr_map, best_del_leaf_lg);
            // std::cerr << "[------DeleteLeaf------]" << "\n";
            // if(best_AOF.or_children_[0].size() == 1 && best_AOF.or_children_[1].size() == 1 && best_AOF.or_children_[0][0].GetContent() == "b2B1" && best_AOF.or_children_[1][0].GetContent() == "b2B2")
            //     target_found = true;
            // if(best_aof_del_leaf.or_children_.size() == 3 && best_aof_del_leaf.or_children_[1][0].GetContent() == "b2B1" && best_aof_del_leaf.or_children_[2][0].GetContent() == "b2B2")
            //     badaof_found = true;

            // if(target_found && badaof_found)
            // {
            //     ConfigBuffer<StateType> all_configs_previous = best_AOF.GenerateAllConfigsWrapper();
            //     std::cerr << "\n";
            //     std::cerr << "Previous AOF size:" << best_AOF.or_children_.size() << "\n";
            //     std::cerr << "Previous Children: " << "\n";
            //     for (int i = 0; i < best_AOF.or_children_.size(); i++)
            //     {
            //         for (int j = 0; j < best_AOF.or_children_[i].size(); j++)
            //             std::cerr << "(" << best_AOF.or_children_[i][j].GetContent() << "," << best_AOF.or_children_[i][j].GetId() << ")";
            //             std::cerr << "\n";
            //     }
            //     std::cerr << "Previous configs: " << "\n";
            //     for(auto it = all_configs_previous.begin(); it != all_configs_previous.end(); it++)
            //     {
            //         for(auto iter = (*it).begin(); iter != (*it).end(); iter++)
            //         {
            //             std::cerr << iter->GetContent() << "_";
            //         }
            //         std::cerr << "\n";
            //     }

            //     // print its posterior: 
            //     std::cerr << "Previous AOF's Posterior Gain: " << std::endl;
            //     GenerateCM(all_configs_previous);
            //     RdBuffer<StateType> reductions_previous = FindReduction(all_configs_previous);
            //     std::cerr << CalculatePosterior(best_AOF, reductions_previous, false) << "\n"; 


            //     //print DeleteLeaf AOF:
            //     std::cerr << "\n";
            //     std::cerr << "DeleteLeaf AOF size:" << best_aof_del_leaf.or_children_.size() << "\n";
            //     std::cerr << "DeleteLeaf Children: " << "\n";
            //     for (int i = 0; i < best_aof_del_leaf.or_children_.size(); i++)
            //     {
            //         for (int j = 0; j < best_aof_del_leaf.or_children_[i].size(); j++)
            //         {
            //             std::cerr << "(" << best_aof_del_leaf.or_children_[i][j].GetContent() << "," << best_aof_del_leaf.or_children_[i][j].GetId() << ")";
            //             std::cerr << "\n";
            //         }
            //     }

            //     // print its posterior: 
            //     ConfigBuffer<StateType> all_configs = best_aof_del_leaf.GenerateAllConfigsWrapper();
            //     std::cerr << "DeleteLeaf AOF's Posterior Gain: " << std::endl;
            //     GenerateCM(all_configs);
            //     RdBuffer<StateType> reductions = FindReduction(all_configs);
            //     std::cerr << CalculatePosterior(best_aof_del_leaf, reductions, false); 
            //     std::cerr << "\n";
            // }
            if (best_del_leaf_lg > best_lg)
            {
                
                best_lg = best_del_leaf_lg;
                best_AOF = best_aof_del_leaf;
                best_reductions = best_del_leaf_rd;
                // std::cerr<<"best AOF or children size: " << best_AOF.or_children_.size();
                // std::cerr << "AOF or children:\n";
                // for(auto or_node : best_AOF.or_children_)
                // {
                //     for(auto state : or_node)
                //     {
                //         std::cerr << "(" << state.GetContent() << ", " << state.GetId() << ")";
                //     }
                //     std::cerr << std::endl;
                // }
                AOFStruct<StateType> test_aof;
                test_aof.or_children_.push_back({Symbolic_State<StateType>("B_t2",true)});
                test_aof.or_children_.push_back({Symbolic_State<StateType>("L_t1",true)});
                // if(test_aof == best_AOF)
                // {
                //     std::cerr<<"data buffer: \n";
                //     for(auto data : this-> data_buffer_)
                //     {
                //         for(auto state : data)
                //         {
                //             std::cerr <<"(" << state.GetContent() << ", "<<state.GetId()<<")";
                //         }
                //         std::cerr << std::endl;
                //     }
                // }
                // std::cerr << "best reduction size: " << best_reductions.begin() -> first.size();
               

                // std::cerr<<"first reduction:\n";
                // for(auto state : best_reductions.begin()->first)
                // {
                //     std::cerr<<"(" << state.GetContent() <<", "<<state.GetId() << ")";
                // }
                // std::cerr << std::endl;
                assert(best_AOF.or_children_.size() == best_reductions.begin()->first.size() );
                
            }
            else if(best_del_leaf_lg == best_lg)
            {

                if(gen() % 2)
                {
                    best_lg = best_del_leaf_lg;
                    best_AOF = best_aof_del_leaf;
                    best_reductions = best_del_leaf_rd;
                    assert(best_AOF.or_children_.size() == best_reductions.begin()->first.size() );
                    
                }
                // std::cerr << "\nWhen equal posterior\n";
                // std::cerr << "Chosen AOF size:" << best_AOF.or_children_.size() << "\n";
                // std::cerr << "Chosen AOF Children: " << "\n";
                // for (int i = 0; i < best_AOF.or_children_.size(); i++)
                // {
                //     for (int j = 0; j < best_AOF.or_children_[i].size(); j++)
                //     {
                //         std::cerr << "(" << best_AOF.or_children_[i][j].GetContent() << "," << best_AOF.or_children_[i][j].GetId() << ")";
                //         std::cerr << "\n";
                //     }
                // }
            }
            // std::cerr << "\n" << "###############################" << "\n" << "###############################" << "\n";

            
            //calculate best AOF obtained by delete an or-node
            // target_found = false;
            // badaof_found = false;
            AOFStruct<StateType> best_aof_del_or = DeleteOr(AOF, best_del_or_rd, prev_map,best_del_or_lg);
            // std::cerr << "[------DeleteOr------]" << "\n";
            // target_found = false;
            // if(best_AOF.or_children_[0].size() == 1 && best_AOF.or_children_[1].size() == 1 && best_AOF.or_children_[0][0].GetContent() == "b2B1" && best_AOF.or_children_[1][0].GetContent() == "b2B2")
            //     target_found = true;
            // if(best_aof_del_or.or_children_.size() == 3 && best_aof_del_or.or_children_[1][0].GetContent() == "b2B1" && best_aof_del_or.or_children_[2][0].GetContent() == "b2B2")
            //     badaof_found = true;
                

            // if(target_found && badaof_found)
            // {
            //     ConfigBuffer<StateType> all_configs_previous = best_AOF.GenerateAllConfigsWrapper();
            //     std::cerr << "\n";
            //     std::cerr << "Previous AOF size:" << best_AOF.or_children_.size() << "\n";
            //     std::cerr << "Previous Children: " << "\n";
            //     for (int i = 0; i < best_AOF.or_children_.size(); i++)
            //     {
            //         for (int j = 0; j < best_AOF.or_children_[i].size(); j++)
            //             std::cerr << "(" << best_AOF.or_children_[i][j].GetContent() << "," << best_AOF.or_children_[i][j].GetId() << ")";
            //             std::cerr << "\n";
            //     }
            //     std::cerr << "Previous configs: " << "\n";
            //     for(auto it = all_configs_previous.begin(); it != all_configs_previous.end(); it++)
            //     {
            //         for(auto iter = (*it).begin(); iter != (*it).end(); iter++)
            //         {
            //             std::cerr << iter->GetContent() << "_";
            //         }
            //         std::cerr << "\n";
            //     }

            //     // print its posterior: 
            //     std::cerr << "Previous AOF's Posterior Gain: " << std::endl;
            //     GenerateCM(all_configs_previous);
            //     RdBuffer<StateType> reductions_previous = FindReduction(all_configs_previous);
            //     std::cerr << CalculatePosterior(best_AOF, reductions_previous, false) << "\n"; 


            //     //print DeleteLeaf AOF:
            //     std::cerr << "\n";
            //     std::cerr << "DeleteOr AOF size:" << best_aof_del_or.or_children_.size() << "\n";
            //     std::cerr << "DeleteOr Children: " << "\n";
            //     for (int i = 0; i < best_aof_del_or.or_children_.size(); i++)
            //     {
            //         for (int j = 0; j < best_aof_del_or.or_children_[i].size(); j++)
            //         {
            //             std::cerr << "(" << best_aof_del_or.or_children_[i][j].GetContent() << "," << best_aof_del_or.or_children_[i][j].GetId() << ")";
            //             std::cerr << "\n";
            //         }

            //     }

            //     // print its posterior: 
            //     ConfigBuffer<StateType> all_configs = best_aof_del_or.GenerateAllConfigsWrapper();
            //     std::cerr << "DeleteOr AOF's Posterior Gain: " << std::endl;
            //     GenerateCM(all_configs);
            //     RdBuffer<StateType> reductions = FindReduction(all_configs);
            //     std::cerr << CalculatePosterior(best_aof_del_or, reductions, false); 
            //     std::cerr << "\n";
            // }
            if (best_del_or_lg > best_lg)
            {
                best_lg = best_del_or_lg;
                best_AOF = best_aof_del_or;
                best_reductions = best_del_or_rd;
                assert(best_AOF.or_children_.size() == best_reductions.begin()->first.size() );    
            }
            else if(best_del_or_lg == best_lg)
            {
                if(gen()%2)
                {
                    best_lg = best_del_or_lg;
                    best_AOF = best_aof_del_or;
                    best_reductions = best_del_or_rd;
                    assert(best_AOF.or_children_.size() == best_reductions.begin()->first.size() );
                    
                }
                //                 std::cerr << "\nWhen equal posterior\n";
                // std::cerr << "Chosen AOF size:" << best_AOF.or_children_.size() << "\n";
                // std::cerr << "Chosen AOF Children: " << "\n";
                // for (int i = 0; i < best_AOF.or_children_.size(); i++)
                // {
                //     for (int j = 0; j < best_AOF.or_children_[i].size(); j++)
                //     {
                //         std::cerr << "(" << best_AOF.or_children_[i][j].GetContent() << "," << best_AOF.or_children_[i][j].GetId() << ")";
                //         std::cerr << "\n";
                //     }
                // }
            }
            // std::cerr << "\n" << "###############################" << "\n" << "###############################" << "\n";

            // target_found = false;
            // badaof_found = false;
            //calculate best AOF obtained by add an or-node
            AOFStruct<StateType> best_aof_add_or = AddOr(AOF, next_frag, best_add_or_rd, next_map,two_configs_map ,best_add_or_lg);
            // std::cerr << "[------AddOr------]" << "\n";
            // target_found = false;
            // if(best_AOF.or_children_[0].size() == 1 && best_AOF.or_children_[1].size() == 1 && best_AOF.or_children_[0][0].GetContent() == "b2B1" && best_AOF.or_children_[1][0].GetContent() == "b2B2")
            //     target_found = true;
            // if(best_aof_add_or.or_children_.size() == 3 && best_aof_add_or.or_children_[1][0].GetContent() == "b2B1" && best_aof_add_or.or_children_[2][0].GetContent() == "b2B2")
            //     badaof_found = true;

            // if(target_found && badaof_found)
            // {
            //     ConfigBuffer<StateType> all_configs_previous = best_AOF.GenerateAllConfigsWrapper();
            //     std::cerr << "\n";
            //     std::cerr << "Previous AOF size:" << best_AOF.or_children_.size() << "\n";
            //     std::cerr << "Previous Children: " << "\n";
            //     for (int i = 0; i < best_AOF.or_children_.size(); i++)
            //     {
            //         for (int j = 0; j < best_AOF.or_children_[i].size(); j++)
            //             std::cerr << "(" << best_AOF.or_children_[i][j].GetContent() << "," << best_AOF.or_children_[i][j].GetId() << ")";
            //             std::cerr << "\n";
            //     }
            //     std::cerr << "Previous configs: " << "\n";
            //     for(auto it = all_configs_previous.begin(); it != all_configs_previous.end(); it++)
            //     {
            //         for(auto iter = (*it).begin(); iter != (*it).end(); iter++)
            //         {
            //             std::cerr << iter->GetContent() << "_";
            //         }
            //         std::cerr << "\n";
            //     }

            //     // print its posterior: 
            //     std::cerr << "Previous AOF's Posterior Gain: " << std::endl;
            //     GenerateCM(all_configs_previous);
            //     RdBuffer<StateType> reductions_previous = FindReduction(all_configs_previous);
            //     std::cerr << CalculatePosterior(best_AOF, reductions_previous, false) << "\n"; 


            //     //print AddOr AOF:
            //     std::cerr << "\n";
            //     std::cerr << "AddOr AOF size:" << best_aof_add_or.or_children_.size() << "\n";
            //     std::cerr << "AddOr Children: " << "\n";
            //     for (int i = 0; i < best_aof_add_or.or_children_.size(); i++)
            //     {
            //         for (int j = 0; j < best_aof_add_or.or_children_[i].size(); j++)
            //         {
            //             std::cerr << "(" << best_aof_add_or.or_children_[i][j].GetContent() << "," << best_aof_add_or.or_children_[i][j].GetId() << ")";
            //             std::cerr << "\n";
            //         }
            //     }

            //     // print its posterior: 
            //     ConfigBuffer<StateType> all_configs = best_aof_add_or.GenerateAllConfigsWrapper();
            //     std::cerr << "AddOr AOF's Posterior Gain: " << std::endl;
            //     GenerateCM(all_configs);
            //     RdBuffer<StateType> reductions = FindReduction(all_configs);
            //     std::cerr << CalculatePosterior(best_aof_add_or, reductions, false); 
            //     std::cerr << "\n";
            // }
            if(best_add_or_lg > best_lg)
            {       
                best_AOF = best_aof_add_or;
                best_lg = best_add_or_lg;
                best_reductions = best_add_or_rd;
                assert(best_AOF.or_children_.size() == best_reductions.begin()->first.size() );

            }
            else if(best_add_or_lg == best_lg)
            {
                if(gen() % 2)
                {
                    best_AOF = best_aof_add_or;
                    best_lg = best_add_or_lg;
                    best_reductions = best_add_or_rd;
                    assert(best_AOF.or_children_.size() == best_reductions.begin()->first.size() );
                    
                }
                // std::cerr << "\nWhen equal posterior\n";
                // std::cerr << "Chosen AOF size:" << best_AOF.or_children_.size() << "\n";
                // std::cerr << "Chosen AOF Children: " << "\n";
                // for (int i = 0; i < best_AOF.or_children_.size(); i++)
                // {
                //     for (int j = 0; j < best_AOF.or_children_[i].size(); j++)
                //     {
                //         std::cerr << "(" << best_AOF.or_children_[i][j].GetContent() << "," << best_AOF.or_children_[i][j].GetId() << ")";
                //         std::cerr << "\n";
                //     }
                // }
            }
            // std::cerr << "\n" << "###############################" << "\n" << "###############################" << "\n";


            // std::cout<<"Best AOF before ADDLEAF : \n";
            
            // for(auto or_node : best_AOF.or_children_)
            // {
            //     for(auto state : or_node)
            //     {
            //         auto content = state.GetContent();    
            //         std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
            //     }
            //     std::cout<<std::endl;
            // }

            // std::cout<<std::endl;
            
            // std::cout<<"Best AOF before addleaf's corresponding rd:\n";
            // for(auto rd : best_reductions)
            // {
                
            //     for(auto state : rd.first)
            //     {
            //         std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
            //     }
            //     std::cout << std:: endl; 
                        
            // }
            // std::cout << std::endl;
            //calculate best AOF obtained by add a leaf node
            // target_found = false;
            // badaof_found = false;
            AOFStruct<StateType> best_aof_add_leaf = AddLeaf(AOF, curr_frag, best_add_leaf_rd,curr_map,best_add_leaf_lg);
            // std::cerr << "[------AddLeaf------]" << "\n";
            // if(best_AOF.or_children_[0].size() == 1 && best_AOF.or_children_[1].size() == 1 && best_AOF.or_children_[0][0].GetContent() == "b2B1" && best_AOF.or_children_[1][0].GetContent() == "b2B2")
            //     target_found = true;
            // if(best_aof_add_leaf.or_children_.size() == 3 && best_aof_add_leaf.or_children_[1][0].GetContent() == "b2B1" && best_aof_add_leaf.or_children_[2][0].GetContent() == "b2B2")
            //     badaof_found = true;

            // if(target_found && badaof_found)
            // {
            //     ConfigBuffer<StateType> all_configs_previous = best_AOF.GenerateAllConfigsWrapper();
            //     std::cerr << "\n";
            //     std::cerr << "Previous AOF size:" << best_AOF.or_children_.size() << "\n";
            //     std::cerr << "Previous Children: " << "\n";
            //     for (int i = 0; i < best_AOF.or_children_.size(); i++)
            //     {
            //         for (int j = 0; j < best_AOF.or_children_[i].size(); j++)
            //             std::cerr << "(" << best_AOF.or_children_[i][j].GetContent() << "," << best_AOF.or_children_[i][j].GetId() << ")";
            //             std::cerr << "\n";
            //     }
            //     std::cerr << "Previous configs: " << "\n";
            //     for(auto it = all_configs_previous.begin(); it != all_configs_previous.end(); it++)
            //     {
            //         for(auto iter = (*it).begin(); iter != (*it).end(); iter++)
            //         {
            //             std::cerr << iter->GetContent() << "_";
            //         }
            //         std::cerr << "\n";
            //     }

            //     // print its posterior: 
            //     std::cerr << "Previous AOF's Posterior Gain: " << std::endl;
            //     GenerateCM(all_configs_previous);
            //     RdBuffer<StateType> reductions_previous = FindReduction(all_configs_previous);
            //     std::cerr << CalculatePosterior(best_AOF, reductions_previous, false) << "\n"; 


            //     //print AddLeaf AOF:
            //     std::cerr << "\n";
            //     std::cerr << "AddLeaf AOF size:" << best_aof_add_leaf.or_children_.size() << "\n";
            //     std::cerr << "AddLeaf Children: " << "\n";
            //     for (int i = 0; i < best_aof_add_leaf.or_children_.size(); i++)
            //     {
            //         for (int j = 0; j < best_aof_add_leaf.or_children_[i].size(); j++)
            //         {
            //             std::cerr << "(" << best_aof_add_leaf.or_children_[i][j].GetContent() << "," << best_aof_add_leaf.or_children_[i][j].GetId() << ")";
            //             std::cerr << "\n";
            //         }
            //     }

            //     // print its posterior: 
            //     ConfigBuffer<StateType> all_configs = best_aof_add_leaf.GenerateAllConfigsWrapper();
            //     std::cerr << "AddLeaf AOF's Posterior Gain: " << std::endl;
            //     GenerateCM(all_configs);
            //     RdBuffer<StateType> reductions = FindReduction(all_configs);
            //     std::cerr << CalculatePosterior(best_aof_add_leaf, reductions, false); 
            //     std::cerr << "\n";
            // }
            if (best_add_leaf_lg > best_lg)
            {
                best_lg = best_add_leaf_lg;
                best_AOF = best_aof_add_leaf;
                best_reductions = best_add_leaf_rd;
                assert(best_AOF.or_children_.size() == best_reductions.begin()->first.size() );
                
            }
            else if(best_add_leaf_lg == best_lg)
            {
                if(gen() % 2)
                {
                    best_lg = best_add_leaf_lg;
                    best_AOF = best_aof_add_leaf;
                    best_reductions = best_add_leaf_rd;

                    // std::cout<<"Best AOF after ADDLEAF : \n";
                    // for(auto or_node : best_AOF.or_children_)
                    // {
                    //     for(auto state : or_node)
                    //     {
                    //         auto content = state.GetContent();    
                    //         std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
                    //     }
                    //     std::cout<<std::endl;
                    // }
                    // std::cout<<std::endl;
                    
                    // std::cout<<"Best AOF after addleaf's corresponding rd:\n";
                    // for(auto rd : best_reductions)
                    // {
                        
                    //     for(auto state : rd.first)
                    //     {
                    //         std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
                    //     }
                    //     std::cout << std:: endl; 
                                
                    // }
                    // std::cout << std::endl;

                    assert(best_AOF.or_children_.size() == best_reductions.begin()->first.size());
                
                }
                // std::cerr << "\nWhen equal posterior\n";
                // std::cerr << "Chosen AOF size:" << best_AOF.or_children_.size() << "\n";
                // std::cerr << "Chosen AOF Children: " << "\n";
                // for (int i = 0; i < best_AOF.or_children_.size(); i++)
                // {
                //     for (int j = 0; j < best_AOF.or_children_[i].size(); j++)
                //     {
                //         std::cerr << "(" << best_AOF.or_children_[i][j].GetContent() << "," << best_AOF.or_children_[i][j].GetId() << ")";
                //         std::cerr << "\n";
                //     }
                // }
            }
            // std::cerr << "\n" << "###############################" << "\n" << "###############################" << "\n";

            // if(exist){
            //     std::cerr<<"attention:!!!!!!!!!after 4 operations!!!!!!!!!!! \n";
            //     for(int j = 0; j <  best_AOF.or_children_.size();++j)
            //     {
            //         std::cerr<<"or_node "<<j<<": ";
            //         for(auto child : best_AOF.or_children_[j])
            //         {
                        
            //             std::cerr<<child.GetContent()<<"_";
                        
            //         }
                        
            //         std::cerr<<std::endl;
            //     }
            // }

            //let the best AOF and corresponding reductinos of this round to be the baseline of next round
            AOF = best_AOF;


            // if(changed)              
            // {
            //     std::cout<<"current best AOF or_children: \n";
            //     for(const auto& or_node : AOF.or_children_)
            //     {
            //         for(const auto& state : or_node)
            //         {
            //             std::cout<<"("<< state.GetContent()<<", "<<state.GetId()<<")";                
            //         }
            //         std::cout<<std::endl;
            //     }
            //     std::cout<<"current best_lg: "<<best_lg<<std::endl
            //         <<"prev best lg: "<<prev_best_lg<<std::endl; 
            //     int b;
            //     std::cin >>b;
            
            // }
            
            
            //generate related k-frags and configs map for next round
            num_of_or_node = AOF.or_children_.size();
            if (num_of_or_node >= all_k_frags.size())
            {
                std::vector<SequenceType<StateType>> temp;
                ConfigsMap<StateType> k_map = GenerateKConfigs(num_of_or_node + 1);
                for (const auto & frag : k_map)
                    temp.push_back(frag.first);
                all_k_frags.push_back(temp);
                all_k_maps.push_back(k_map);

            }

            curr_frag = all_k_frags[num_of_or_node - 1];
            next_frag = all_k_frags[num_of_or_node];
            prev_map = (num_of_or_node > 1) ? all_k_maps[num_of_or_node - 2]: all_k_maps[0];
            curr_map = all_k_maps[num_of_or_node - 1];
            next_map = all_k_maps[num_of_or_node];
        }
        
       
        //store each round's best result and sample from them
        best_lg_of_all_rds.push_back(best_lg);
        best_AOF_of_all_rds.push_back(best_AOF);
        best_reductions_of_all_rds.push_back(best_reductions);
        assert(best_AOF.or_children_.size() == best_reductions.begin()->first.size() );
        
        --num_of_iter;
    }

    // NUM_OF_BIGRAM.push_back(Combination(num_of_two_configs,2));
    // NUM_OF_VARIATION.push_back((double)NUM_OF_VARIATION_EACH_ROUND/iterations);
    // TOTAL_VARIATIONS += NUM_OF_VARIATION_EACH_ROUND;
    // std::cerr<<"Full AOF of all rounds: \n";
    // for(auto aof : best_AOF_of_all_rds )
    // {
    //     std::cerr<<std::endl;
    //     for(auto or_node : aof.or_children_)
    //     {
    //         for(auto state : or_node)
    //         {
    //             std::cerr<<state.GetContent()<<"_";
    //         }
    //         std::cerr<<std::endl;
    //     }
    // }

    // std::cout<<"Best AOF of all rounds: \n";
    // for(auto aof : best_AOF_of_all_rds )
    // {
    //     for(auto or_node : aof.or_children_)
    //     {
    //         for(auto state : or_node)
    //         {
    //             auto content = state.GetContent();    
    //             std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
    //         }
    //         std::cout<<std::endl;
    //     }

    //     std::cout<<std::endl;
    // }
    // std::cout<<"Best AOF's corresponding rd:\n";
    // for(auto rd : best_reductions_of_all_rds)
    // {
    //     for(auto iter : rd )
    //     {
    //         for(auto state : iter.first)
    //         {
    //             std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
    //         }
    //         std::cout << std:: endl; 
    //     }
    //     std::cout<<std::endl;
       
    // }

    //remove repeated AOF
    std::unordered_set<AOFStruct<StateType> > unique_aofs;
    std::vector<bool> is_unique(best_AOF_of_all_rds.size(),true);
    for(int i = 0; i < best_AOF_of_all_rds.size();++i)
    {
        if(unique_aofs.find(best_AOF_of_all_rds[i]) == unique_aofs.end())
            unique_aofs.insert(best_AOF_of_all_rds[i]);
        else
            is_unique[i] = false;
    }

    std::vector<AOFStruct<StateType> > unique_AOF_of_all_rds;
    std::vector<double> unique_lg_of_all_rds;
    std::vector<RdBuffer<StateType> > unique_reduction_of_all_rds;
    for(int i = 0; i < is_unique.size();++i)
    {
        if(is_unique[i])
        {
            unique_AOF_of_all_rds.push_back(best_AOF_of_all_rds[i]);
            unique_lg_of_all_rds.push_back(best_lg_of_all_rds[i]);
            unique_reduction_of_all_rds.push_back(best_reductions_of_all_rds[i]);
        }

    }

    

    //  std::cerr<<"Unique AOF of all rounds: \n";
    // for(auto aof : unique_AOF_of_all_rds )
    // {
    //     std::cerr<<std::endl;
    //     for(auto or_node : aof.or_children_)
    //     {
    //         for(auto state : or_node)
    //         {
    //             std::cerr<<state.GetContent()<<"_";
    //         }
    //         std::cerr<<std::endl;
    //     }
    // }

    //remove AOF that are in the current checkpoint's children
    std::vector<int> remove_helper(unique_AOF_of_all_rds.size(),0);
    SN_Ptr cpt_parent = this->GetParent();
    if(cpt_parent != 0)
    {
        std::vector<SN_Ptr> cpt_siblings = cpt_parent->GetChildren();
        if(!cpt_siblings.empty())
        {
            std::vector<AOFStruct<StateType> > siblings_AOF;
            for(const auto & sibling : cpt_siblings)
                siblings_AOF.push_back(static_cast<Checkpoint *>(sibling.get())->AOFragment_);

            for(const auto & sibling : siblings_AOF)
            {
                for(int i = 0; i < unique_AOF_of_all_rds.size(); ++i)
                {
                    if(unique_AOF_of_all_rds[i] == sibling)
                        remove_helper[i] = -1;
                    
                }
            
            }

            std::vector<AOFStruct<StateType> > best_AOF_for_sample;
            std::vector<double> best_lg_for_sample;
            std::vector<RdBuffer<StateType> > best_rd_for_sample;

            for(int i = 0; i < remove_helper.size();++i)
            {
                if(remove_helper[i] != -1)
                {
                    best_AOF_for_sample.push_back(unique_AOF_of_all_rds[i]);
                    best_lg_for_sample.push_back(unique_lg_of_all_rds[i]);
                    best_rd_for_sample.push_back(unique_reduction_of_all_rds[i]);
                }
            }

            unique_AOF_of_all_rds = best_AOF_for_sample;
            unique_lg_of_all_rds = best_lg_for_sample;
            unique_reduction_of_all_rds = best_rd_for_sample;
        }
    
    }

    //if the likelihood gain of all best AOF found is less than 1, return empty reductions
    bool find_good_aof = false;
    // std::cerr<<"all found aof's post gain:\n";
    for(const auto & postgain : unique_lg_of_all_rds)
    {
        // std::cerr<< postgain;
        if(postgain >= THRESHOLD)
        {
            find_good_aof = true; 
            break;
        }
    }
    // std::cerr<<std::endl;

    if(!find_good_aof)
    {
        this->CM_ = ori_CM;
        this->AOF_cache_map_.clear();

        // if(this->posterior_ > BEST_POSTERIOR.back())
        // {
        //     BEST_POSTERIOR.push_back(this->posterior_);
        //     KL_DIVERGENCE.push_back(true);
        // }
        // else
        // {
        //     BEST_POSTERIOR.push_back(BEST_POSTERIOR.back());
        //     KL_DIVERGENCE.push_back(false);
        // }
        std::cerr<<"cannot find good aof\n";
    
        // OutputToFile("../AOG_Sample_Results/sample_"+std::to_string(BEST_POSTERIOR.size()-1),this->graph_,this->graph_.GetRoot(),SAMPLING);                    
        return false;
    }
    
    int best_idx;
    if (use_stochastic)
    {
        std::vector<double> exp_lg_for_sample;
        auto iter = std::min_element(unique_lg_of_all_rds.begin(), unique_lg_of_all_rds.end());
        double min = unique_lg_of_all_rds[iter - unique_lg_of_all_rds.begin()];

        for(const auto & lg : unique_lg_of_all_rds)
        {
            exp_lg_for_sample.push_back(exp(lg-min));
            // std::cerr<<"log posterior: "<<lg<<" ";            
            // std::cerr<<"exponential posterior: "<<exp(lg - min)<<" ";
            // exp_lg_for_sample.push_back(lg);

            // std::cerr<<"normalized exponential posterior: "<<exp(lg - min)<<" ";
            // std::cerr<<"\n";

        }
        // int g;
        // std::cin >>g;


        //else make a distribution of all possible best AOF
        std::discrete_distribution<int> lg_dis(exp_lg_for_sample.begin(),exp_lg_for_sample.end());

        //select best AOF
        best_idx = lg_dis(gen);
    }
    else
    {
        auto it = std::max_element(unique_lg_of_all_rds.begin(), unique_lg_of_all_rds.end());
        best_idx = it-unique_lg_of_all_rds.begin();
    }

    //record the best posterior gain found as statistics
    double max_post_gain = *std::max_element(unique_lg_of_all_rds.begin(),unique_lg_of_all_rds.end());
    // compare to the last recorded best posterior
    // if(this->posterior_ + max_post_gain > BEST_POSTERIOR.back())
    // {
    //     // if(max_post_gain + this->posterior_ > 0)
    //     // {
    //     //     std::cerr<<"current posterior larger than previous!\n";
    //     //     std::cerr<<"current posterior: "<<this->posterior_<<std::endl;
    //     //     std::cerr<<"max post gain: "<<max_post_gain<<std::endl; 
    //     //     std::cerr<<"previous posterior"<< BEST_POSTERIOR.back()<<std::endl;   
    //     //     std::cerr<< this->posterior_ + max_post_gain<<std::endl;

    //     //     std::cerr<<"unique post gain of all rounds:\n";
    //     //     for(auto lg : unique_lg_of_all_rds)
    //     //     {
    //     //         std::cerr<<lg<<" ";
    //     //     }
    //     //     std::cerr<<std::endl;
    //     //     int a;
    //     //     std::cin >> a;
    //     // }
        
    //     BEST_POSTERIOR.push_back(this->posterior_+max_post_gain);
    //     KL_DIVERGENCE.push_back(true);
    //     // int a;
    //     // std::cin >>a;
    // }
    // else
    // {
    //     // std::cerr<<"previous larger!\n";
    //     // std::cerr<< BEST_POSTERIOR.back();
    //     BEST_POSTERIOR.push_back(BEST_POSTERIOR.back());
    //     KL_DIVERGENCE.push_back(false);
    //     // int a;
    //     // std::cin >> a;

    // }

   
    this->AOFragment_ = unique_AOF_of_all_rds[best_idx];
    //remove weight 0 state
    for(int i = 0; i < this->AOFragment_.weights_.size(); ++i)
    {
        for(auto iter = this->AOFragment_.weights_[i].begin(); iter != this->AOFragment_.weights_[i].end();)
        {
            if((*iter) == 0)
            {
                if(this->AOFragment_.weights_[i].size() == 1)
                {
                    std::cerr<<"find the only weight 0 AOF!\n";
                    throw std::exception();
                }
               iter = this->AOFragment_.weights_[i].erase(iter);
               this->AOFragment_.or_children_[i].erase( iter - this->AOFragment_.weights_[i].begin() + this->AOFragment_.or_children_[i].begin());
            }
            else
                ++iter;
        }
    }

    for(auto weights : this->AOFragment_.weights_)
    {
    for(auto weight : weights)
        {
        if(weight == 0)
            {
            std::cerr<<"find you after removing weight 0 state in aof\n";
            throw std::exception();
            }

        }

    }
    
    // std::vector<Symbolic_State<std::string> > test_or_child_1 = {Symbolic_State<std::string> ("1",true), Symbolic_State<std::string>("2",true)};
    // std::vector<Symbolic_State<std::string> > test_or_child_2 = {Symbolic_State<std::string>("6",true)};
    // std::vector<Symbolic_State<std::string> > test_or_child_3 = {Symbolic_State<std::string>("7",true)};
    // std::vector<Symbolic_State<std::string> > test_or_child_4 = {Symbolic_State<std::string>("8",true)};
    // AOFStruct<StateType> test_aof;
    // test_aof.or_children_ = {test_or_child_1, test_or_child_2, test_or_child_3, test_or_child_4};
    // if(test_aof == this->AOFragment_)
    // {
    //     unsigned best_or_size = this->AOFragment_.or_children_.size();
    //     ConfigBuffer<StateType> best_configs = this->AOFragment_.GenerateAllConfigsAndWeights(all_k_maps[best_or_size-1]);
    //     RdBuffer<StateType> best_rd1 = GenerateCMAndReductions(best_configs,this->data_buffer_);

    //     std::cout << "catch !!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    //     //std::cout << "factor 1: " << this->CalculateNewFactor1(this->AOFragment_, best_rd1) << std::endl;

    //     AOFStruct<StateType> test2_aof;
    //     test2_aof.and_node_ = Symbolic_State<StateType>(test2_aof.id_);
	//     //++(test2_aof.id_);
    //     std::vector<Symbolic_State<std::string> > test2_or_child_1 = {Symbolic_State<std::string> ("1",true), Symbolic_State<std::string>("2",true)};
    //     test2_aof.or_children_ = {test2_or_child_1};
    //     test2_aof.weights_ = {std::vector<double>(test2_or_child_1.size(), 0)};
    //     ConfigBuffer<StateType> test2_configs = test2_aof.GenerateAllConfigsAndWeights(one_configs_map);
    //     RdBuffer<StateType> test2_reductions = this->GenerateCMAndReductions(test2_configs,this->data_buffer_);
    //     std::cout << "test aof:\n";
    //     std::cout << "Posterior gain with beta: " << this->CalculatePosterior(test2_aof, test2_reductions, true) << std::endl;
    //     CalculatePosterior(test2_aof, test2_reductions, false);
        

    // }


    std::cerr<<"The AOF to be added:\n";
    std::cerr<<"and node: ("<<this->AOFragment_.and_node_.GetContent()<<", "<<this->AOFragment_.and_node_.GetId()<<")\n";
    for(const auto & or_node : this->AOFragment_.or_children_)
    {
        for(const auto & state : or_node)
            std::cerr<<"( "<<state.GetContent()<<", "<<state.GetId()<<")";
        std::cerr<<std::endl;
    }
    std::cerr<<std::endl;

    unsigned best_or_size = this->AOFragment_.or_children_.size();
    // std::cerr<<"it's post gain: "<<unique_lg_of_all_rds[best_idx]<<std::endl;

   
    
    // std::cerr<<"the number of neg post gain: "<<this->num_of_neg_post_gain_<<std::endl;
    // if(this->AOFragment_.and_node_.GetId() == 5714)
    // {
    //     std::cerr <<"find 5714\n";
    //     for(auto weights : this->AOFragment_.weights_)
    //     {
    //         for(auto weight : weights)
    //         {
    //             std::cerr <<weight <<" ";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     std::cerr << "the current data buffer:\n";
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }

    //      std::cerr<<"all best aof: \n";
    //     for(auto unique_aof : unique_AOF_of_all_rds)
    //     {
    //         for(auto or_node : unique_aof.or_children_)
    //         {
    //             for(auto state : or_node)
    //             {
    //                 std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //             }
    //             std::cerr<<std::endl;
    //         }
    //     }
    //     std::vector<Symbolic_State<std::string> > test_or_child_1 = {Symbolic_State<std::string> ("J_con1",true), Symbolic_State<std::string>("J_con2",true)};
    //     std::vector<Symbolic_State<std::string> > test_or_child_2 = {Symbolic_State<std::string>("J_move1",true), Symbolic_State<std::string>("J_move2",true)};
    // // std::vector<Symbolic_State<std::string> > test_or_child_3 = {Symbolic_State<std::string>("7",true)};
    // // std::vector<Symbolic_State<std::string> > test_or_child_4 = {Symbolic_State<std::string>("8",true)};
    //     AOFStruct<StateType> test_aof;
    //     test_aof.or_children_ = {test_or_child_1, test_or_child_2};
    //     test_aof.weights_ = {std::vector<double>(test_or_child_1.size(), 0), std::vector<double>(test_or_child_1.size(), 0)};
    //     ConfigBuffer<StateType> test_configs = test_aof.GenerateAllConfigsAndWeights(two_configs_map);
    //     RdBuffer<StateType> test_reductions = this->GenerateCMAndReductions(test_configs,this->data_buffer_);

    //     std::cout << "test aof:\n";
    //     for(const auto & or_node : test_aof.or_children_)
    //     {
    //         for(const auto & state : or_node)
    //             std::cerr<<"( "<<state.GetContent()<<", "<<state.GetId()<<")";
    //         std::cerr<<std::endl;
    //     }

    //     std::cerr<<"test aof's weight: \n";
    //     for(const auto & weights : test_aof.weights_)
    //     {
    //         for(const auto & weight : weights)
    //             std::cerr<< weight<<" ";
    //         std::cerr<<std::endl;
    //     }

    //     std::cout << "Posterior gain with beta: " << this->CalculatePosterior(test_aof, test_reductions, true) << std::endl;

    //     this->graph_.OutputGraph("peeker",PATH,true,SIMPLIFY);

    //     this->AOFragment_.or_children_[1] = { Symbolic_State<std::string>("J_move2",true)};
    //     this->AOFragment_.weights_[1] = std::vector<double>(this->AOFragment_.or_children_[1].size(), 0);
    //     std::cout << "new test aof:\n";
    //     for(const auto & or_node : this->AOFragment_.or_children_)
    //     {
    //         for(const auto & state : or_node)
    //             std::cerr<<"( "<<state.GetContent()<<", "<<state.GetId()<<")";
    //         std::cerr<<std::endl;
    //     }
    //     ConfigBuffer<StateType> changed_configs = this->AOFragment_.GenerateAllConfigsAndWeights(three_configs_map);
    //     RdBuffer<StateType> changed_reductions = this->GenerateCMAndReductions(changed_configs,this->data_buffer_);
    //     std::cout << "Posterior gain with beta: " << this->CalculatePosterior(this->AOFragment_, changed_reductions, true) << std::endl;


    //     int a;
    //     std::cin >>a;
    // }
    // if (this->AOFragment_.or_children_.size() == 2
    // && this->AOFragment_.or_children_[0].size() == 1
    // && this->AOFragment_.or_children_[1].size() == 1
    // && this->AOFragment_.or_children_[0][0].GetId() == 9489
    // && this->AOFragment_.or_children_[1][0].GetId() == 29198
    // )
    // {

    //     std::cerr << "9489 29198 found!\n";
    //     this->graph_.OutputGraph("before_bug",PATH,false);
    // }

    //  if (this->AOFragment_.or_children_.size() == 2
    // && this->AOFragment_.or_children_[0].size() == 1
    // && this->AOFragment_.or_children_[1].size() == 1
    // && this->AOFragment_.or_children_[0][0].GetId() == 29198
    // && this->AOFragment_.or_children_[1][0].GetId() == 29849
    // )
    // {
    //     std::cerr << "29198 29848 found!\n";
    //     this->graph_.OutputGraph("before_bug",PATH,false);
        
    // }
    
    RdBuffer<StateType> best_rd;
    ConfigBuffer<StateType> best_configs = this->AOFragment_.GenerateAllConfigsAndWeights(all_k_maps[best_or_size-1]);
    
    if(CACHE)
    {
        auto iter = this->AOF_cache_map_.find(this->AOFragment_.or_children_);
        //if current AOF's corresponding CM and reductions have been cached before, use them directly
        if(iter != this->AOF_cache_map_.end())
        {
            // std::cerr << "cache find!\n";
            // int a;
            // std::cin >> a;
            Cache<StateType> cache = this->AOF_cache_map_[this->AOFragment_.or_children_];
            best_rd = cache.rd_cache_;
            this->CM_ = cache.cm_cache_;
        }
        else
            best_rd =  GenerateCMAndReductions(best_configs,this->data_buffer_);
            

        
    }
    else
        best_rd =  GenerateCMAndReductions(best_configs,this->data_buffer_);


    
    double best_post_gain = CalculatePosterior(this->AOFragment_, best_rd, false);

    // std::cerr<<"actual post gain of the added aof: "<<best_post_gain<<std::endl;
    this->posterior_ +=  best_post_gain;
    std::cerr << "posterior after generate: " << this->posterior_<<std::endl;

   
    
    // for(auto or_node : this->AOFragment_.or_children_)
    // {
    //     for(auto state : or_node)
    //     {
    //         if(state.GetId() == 73861)
    //         {
    //             std::cerr<<"find you in or_children!\n";
    //             this->graph_.OutputGraph("add_73861_before",PATH,true,SIMPLIFY);
    //             int a;
    //             std::cin >>a;
    //         }
    //     }
    // }
    // int q;
    // std :: cin >> q;

    // if(this->AOFragment_.and_node_.GetId() == 73841)
    // {
    //     std::cerr<<"find 73841";
    //     this->graph_.OutputGraph("73841_before",PATH,true,true);
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     int a;
    //     std::cin >> a;
    // }
    //test AOF
    // if (this->AOFragment_.and_node_.GetId() == 72)
    // {
    //     std::cerr<<"Best AOF chosen from Generate():\n";
    //     for(auto or_node : this->AOFragment_.or_children_)
    //     {
    //         for(auto state : or_node)
    //             std::cerr<< state.GetContent()<<"_";
    //         std::cerr<<std::endl;
    //     }
    //     std::cerr<<"the current posterior gain:"<<best_post_gain<<std::endl;

    //     AOFStruct<StateType> test_aof;
    //     // test_aof.or_children_[0].push_back(Symbolic_State<std::string>("Square",true));
        
    //     std::vector<Symbolic_State<std::string> > test_or_child_1;        
    //     std::vector<Symbolic_State<std::string> > test_or_child_2;
    //     std::vector<Symbolic_State<std::string> > test_or_child_3;

    //     // test_aof.or_children_.push_back(this->AOFragment_.or_children_[1]);

    //     // std::vector<Symbolic_State<std::string> > test_or_child_3;
        
    //     test_or_child_1 = {Symbolic_State<std::string> ("J_m1",true)};  
    //     test_or_child_2 = {Symbolic_State<std::string> ("J_m2",true)};
    //     test_or_child_3 = {Symbolic_State<std::string> ("J_m3",true)};

    //     // test_or_child_2.push_back(Symbolic_State<std::string> ("Round",true));
    //     // test_or_child_2.push_back(Symbolic_State<std::string> ("Square",true));
    //     test_aof.or_children_.push_back(test_or_child_1);
    //     test_aof.or_children_.push_back(test_or_child_2);
    //     test_aof.or_children_.push_back(test_or_child_3);
    //     // test_aof.or_children_ = {this->AOFragment_.or_children_[0]};
    //     // test_aof.or_children_[1].push_back(Symbolic_State<std::string> ("Blue",true));
    //     test_aof.and_node_ = Symbolic_State<std::string>(1000);
    //     // test_or_child_1.push_back(Symbolic_State<std::string> ("c1A3",true));
    //     test_aof.weights_.push_back(std::vector<double>(1,0)) ;
    //     test_aof.weights_.push_back(std::vector<double>(1,0)) ;
    //     test_aof.weights_.push_back(std::vector<double>(1,0)) ;


    //     std::cerr<<"Desired AOF:\n";
    //     for(auto or_node : test_aof.or_children_)
    //     {
    //         for(auto state : or_node)
    //         {
    //             std::cerr<<"(" <<state.GetContent()<<", "<<state.GetId()<<")";

    //         }
    //         std::cerr<<std::endl;
    //     }

    //     ConfigBuffer<std::string> test_configs = test_aof.GenerateAllConfigsAndWeights(three_configs_map);
    //     std::cerr<<"the test aof's weight is: \n";
    //     for(auto weights : test_aof.weights_)
    //     {
    //         for(auto weight : weights)
    //             std::cerr<<weight<<" ";
    //         std::cerr<<std::endl;
    //     }
    //     RdBuffer<std::string> test_reductions = GenerateCMAndReductions(test_configs,this->data_buffer_);
    //     double test_lg = CalculatePosterior(test_aof,test_reductions);
    //     std::cerr<<"the test aof post gain: "<<test_lg<<std::endl;
    //     int a;
    //     std:: cin >> a;
    // }
    // {
   
    //     // try
    //     // {
    //     //     this->graph_.GetVertexIdByState(Symbolic_State<std::string> (225));
    //     // } catch (...)
    //     // {
    //     //     std::cerr << "or node hasn't been formed yet\n";
    //     // }

    //     // test_or_child_1.push_back(Symbolic_State<std::string> (83));
    //     // test_or_child_2.push_back(Symbolic_State<std::string> (225));        

    //     test_aof.or_children_.push_back(test_or_child_1);       
    //     // test_aof.or_children_.push_back(test_or_child_2);
    //     // test_aof.or_children_.push_back(test_or_child_3);
    
    //     test_aof.weights_.push_back(std::vector<double>(2,1)) ;
   

    // }



    /** std::cerr<<"before generateCM()\n"; */
    // std::cerr<<"Size of chosen AOF's Or children: "<<best_or_size<<std::endl;
    // std::cerr << "Chosen AOF's or children:" << this->AOFragment_.or_children_.size() << "\n";
    // for (int i = 0; i < this->AOFragment_.or_children_.size(); i++)
    // {
    //     for (int j = 0; j < this->AOFragment_.or_children_[i].size(); j++)
    //         std::cerr << "(" << this->AOFragment_.or_children_[i][j].GetContent() << "," << this->AOFragment_.or_children_[i][j].GetId() << ")";
    //         std::cerr << "\n";
    // }
    
    RRMap<StateType> rrmap = FindRRMap(this->AOFragment_);
    // for(auto entry = rrmap.map_.begin(); entry != rrmap.map_.end(); ++entry)
	// {
	// 	for(auto state : entry->first)
	// 	{
	// 		std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
	// 	}
	// 	std::cerr<<"->\n";
	// 	for(auto data : entry->second)
	// 	{
	// 		for(auto state : data)
	// 		{
	// 			std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
	// 		}
	// 		std::cerr<<std::endl;
	// 	}
	// }
    Replacement<StateType> best_replacement;
    best_replacement.rrmap_ = rrmap;
    best_replacement.reductions_ = unique_reduction_of_all_rds[best_idx];
        


    
    
    // test_aof.weights_.push_back(std::vector<int>(2,0)  ;
        // test_aof.weights_.push_back(std::vector<int>(1,0) );
    //     test_aof.and_node_ = Symbolic_State<std::string>(2000000);

    //     ConfigBuffer<std::string> test_configs = test_aof.GenerateAllConfigsAndWeights(one_configs_map);
    //     GenerateCM(test_configs);
    //     RdBuffer<std::string> test_reductions = FindReduction(test_configs);
    //     double test_lg = CalculatePosterior(test_aof,test_reductions);


    // std::cout << "REDUCTION\n";
    // for (auto it = red.begin(); it != red.end(); it++)
    // {
    //     auto sequence_to_be_replaced = it->first;
    //     for (int i = 0; i < sequence_to_be_replaced.size(); i++)
    //             std::cout << "(" << sequence_to_be_replaced[i].GetContent() << "," << sequence_to_be_replaced[i].GetId() << ")";
    //     std::cout << "\n";
    // }
    // std::cout << "AOF size:" << this->AOFragment_.or_children_.size() << "\n";
    // for (int i = 0; i < this->AOFragment_.or_children_.size(); i++)
    // {
    //     for (int j = 0; j < this->AOFragment_.or_children_[i].size(); j++)
    //         std::cout << "(" << this->AOFragment_.or_children_[i][j].GetContent() << "," << this->AOFragment_.or_children_[i][j].GetId() << ")";
    //     std::cout << "\n";
    // }
    // std::cout << "----------";
    //if mistake_allowed reaches 0, return false
    /** std::cerr<<"the size of reductions returned from Generate(): "<<best_reductions.size()<<std::endl; */
    // cout <<"posterior gain: "<<unique_lg_of_all_rds[best_idx];
    UpdateGraphAndBuffer(best_replacement);

    if(NEG_GAIN_COUNT != -1)
    {
        if(best_post_gain < 0)
            ++num_of_neg_post_gain_;
        else
            num_of_neg_post_gain_ = 0;
        std::cerr<<"current neg gain count: "<<num_of_neg_post_gain_<<std::endl;
        if(this->posterior_ > this->best_posterior_)
        {
            this->best_posterior_ = this->posterior_;
            this->best_aog_ = T_AOG<StateType>(this->graph_);
            if(!OFFLINE)
                this->best_buffer_ = this->data_buffer_;
            // std::cerr<<"the best aog's root: "<< this->best_aog_.GetRoot()<<std::endl;
        }
    }
    
    // for(auto or_node : this->AOFragment_.or_children_)
    // {
    //     for(auto state : or_node)
    //     {
    //         if(state.GetId() == 73861)
    //         {
    //             std::cerr<<"find you in or_children!\n";
    //             this->graph_.OutputGraph("add_73861_after",PATH,true,SIMPLIFY);
    //             int a;
    //             std::cin >>a;
    //         }
    //     }
    // }
    // if(this->AOFragment_.and_node_.GetId() == 3536)
    // {
    //     std::cerr<<"find 3536";
    //     this->graph_.OutputGraph("3536_after",PATH,true,true);
    //     int a;
    //     std::cin >> a;
    // }
    // OutputToFile("../AOG_Sample_Results/sample_"+std::to_string(BEST_POSTERIOR.size()-1),this->graph_,this->graph_.GetRoot(),SAMPLING);
    
    // std::cerr<<"---------------------Outside Generate()-----------------------------\n";
    //  bool wrong_buffer = false;
    // for(auto data : this->data_buffer_)
    // {
    //     if(data[data.size()-1].GetId() == 261)
    //     {
    //         std::cerr<<"wrong data buffer in generate\n";
    //         wrong_buffer = true;
    //         break;
    //     }
    // }
    // if(this->AOFragment_.and_node_.GetId() == 73841)
    // {
    //     std::cerr<<"find 73841";
    //     this->graph_.OutputGraph("73841_after",PATH,true,SIMPLIFY);
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }

    //     auto weights_map = this->graph_.GetOutEdgeWeights(0,false);
    //     bool find_weight = false;
    //     for(auto iter : weights_map)
    //     {
    //         auto dummy_children = this->graph_.ChildrenVertices(iter.first);
    //         for(auto child : dummy_children)
    //         {
    //             if(this->graph_.GetStateByVertexId(child).GetId() == 73861)
    //             {
    //                 find_weight = true;
    //                 break;
    //             }
    //         }
    //         if(find_weight)
    //         {
    //             std::cerr<<"73861 weight: "<<iter.second<<std::endl;
    //             break;
    //         }
    //     }

    //     int a;
    //     std::cin >> a;
    // }
    // if(wrong_buffer)
    // {
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     this->graph_.OutputGraph("wrong_buffer",PATH,true,false);
    //     int a;
    //     std::cin >>a;
    // }

    this->AOF_cache_map_.clear();
    return true;

}


template<class StateType>
void Online_Learner<StateType>::Checkpoint::UpdateGraphAndBuffer(const Replacement<StateType> &replacement)
{
    // std::cerr<<"inside update graph and buffer\n";
    // if(this->AOFragment_.and_node_.GetId() == 18760)
    // {
    //     std::cerr<<"before update 18760\n";
    //     this->graph_.OutputGraph("before_18760",PATH,true,SIMPLIFY);
    //     int a;
    //     std::cin >>a;
    // }
    // this->graph_.OutputGraph("before_update",PATH,false);
    // auto third_level_map = this->GetThirdLevel();
    // std::cerr<<"Third level in updategraph and buffer:\n";
    // for(auto iter = third_level_map.begin(); iter != third_level_map.end();++iter)
    // {
    //     for(auto state : iter->first)
    //     {
    //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //     }
    //     std::cerr<<std::endl;
    // }

    // if(this->AOFragment_.and_node_.GetId() == 3207)
    // {
    //     std::cout << "*****DATA BUFFER IS*****\n";
    //     for (int i = 0; i < this->data_buffer_.size(); i++)
    //     {
    //         std::cout << "[" << i << "]:";
    //         for (int j = 0; j < this->data_buffer_[i].size(); j++)
    //             std::cout << "(" <<  this->data_buffer_[i][j].GetContent() << "," << this->data_buffer_[i][j].GetId() << ")";
    //         std::cout << "\n";
    //     }
    //     auto rrmap = replacement.rrmap_;
    //     std::cout << "*****RRMap IS*****\n";
    //     for (auto it = rrmap.map_.begin(); it != rrmap.map_.end(); it++)
    //     {
    //         for (int i = 0; i < it->first.size(); i++)
    //         {
    //             std::cout << "(" << it->first[i].GetContent() << "," << it->first[i].GetId() << ")";
    //         }
    //         std::cout << "->\n";
    //         for (int i = 0; i < it->second.size(); i++)
    //         {
    //             for (int j = 0; j < it->second[i].size(); j++)
    //                 std::cout << "\t(" << it->second[i][j].GetContent() << "," << it->second[i][j].GetId() << ")";
    //             std::cout << "\n";
    //         }
    //     }

    //     std::cerr<<"find you as and_node!\n";
    //     std::cerr<<"is_delete: "<<replacement.rrmap_.is_delete_<<std::endl;
    //     this->graph_.OutputGraph("add3207_before",PATH,false);
    //     std::cerr<<"and node: ("<<this->AOFragment_.and_node_.GetContent()<<", "<<this->AOFragment_.and_node_.GetId()<<")"<<std::endl;
    //     for(auto nodes: this->AOFragment_.or_children_)
    //     {
    //         for(auto node_state : nodes)
    //         {
    //             std::cerr<<"("<<node_state.GetContent()<<", "<<node_state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }

    //     int a;
    //     std::cin >> a;
    // }
    // for(auto or_node : this->AOFragment_.or_children_)
    // {
    //     for(auto state : or_node)
    //     {
    //         if(state.GetId() == 3207)
    //         {
    //             std::cerr<<"find you from or children!\n";
    //             std::cerr<<"is_delete: "<<replacement.rrmap_.is_delete_<<std::endl;                
    //             this->graph_.OutputGraph("add3207_before",PATH,false);
    //             std::cerr<<"and node: ("<<this->AOFragment_.and_node_.GetContent()<<", "<<this->AOFragment_.and_node_.GetId()<<")"<<std::endl;
    //             for(auto nodes: this->AOFragment_.or_children_)
    //             {
    //                 for(auto node_state : nodes)
    //                 {
    //                     std::cerr<<"("<<node_state.GetContent()<<", "<<node_state.GetId()<<")";
    //                 }
    //                 std::cerr<<std::endl;
    //             }
    //             int a;
    //             std:: cin >> a;                
    //         }
    //     }
    // }

    // std::cout << "-----Entering update buffer-----\n";
    std::vector<Symbolic_Rule<StateType> > all_rules = this->graph_.GetRules();
    // std::cout << "All rules in graph:\n";    
    // for (auto rule : all_rules)
    // {
    //     Symbolic_State<StateType> source = rule.GetSource();
        
    //     std::cout << "(" << source.GetContent() << "," << source.GetId() << ")" << "->";        
    //     std::vector<Symbolic_State<StateType> > results = rule.GetResults();
    //     for (int i = 0; i < results.size(); i++)
    //     {
    //         std::cout << "(" << results[i].GetContent() << "," << results[i].GetId() << ")";
    //     }
    //     std::cout << "\n";
        
    // }

    RdBuffer<StateType> reduction = replacement.reductions_;
    RRMap<StateType> rrmap = replacement.rrmap_;
    
    // std::cerr<<"aof added to the graph:\n";
    // for(auto or_node : this->AOFragment_.or_children_)
    // {
    //     for(auto state : or_node)
    //     {
    //         std::cerr<<"("<<state.GetContent()<<","<<state.GetId()<<")";
    //     }
    //     std::cerr<<std::endl;
    // }
    // AOFStruct<StateType> test_aof;
    // std::vector<Symbolic_State<std::string> > test_or_child_1;        
    
    // test_or_child_1.push_back(Symbolic_State<std::string> ("Hawk",true));       
    // test_or_child_1.push_back(Symbolic_State<std::string> ("Roman",true));
    // test_or_child_1.push_back(Symbolic_State<std::string> ("Bumpy",true));

    
    // test_aof.or_children_.push_back(test_or_child_1);       
    
    
    // test_aof.weights_.push_back(std::vector<double>(3,0)) ;
        

    // if(this->AOFragment_ == test_aof)
    // {
    //     std::cerr<<"before 693\n";
    //     std::cerr<<"data buffer: \n";
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     this->graph_.OutputGraph("before_add_nose",PATH,false);

        // int a;
        // std::cin >>a;
    // }

    // if(this -> AOFragment_.and_node_.GetId() == 24)
    // {
    //     std::cerr<<"before 24\n";
    //     std::cerr<<"size of or_children: "<<this->AOFragment_.or_children_.size()<<  std::endl;
    //     this->graph_.OutputGraph("visualize_24_before",PATH,true);
    //     // int a;
    //     // std::cin >>a;
    // }

    // std::cout<<"rrmap: \n";
    // for(auto iter = rrmap.map_.begin(); iter != rrmap.map_.end(); ++ iter)
    // {
    //     for(auto ori : iter->first)
    //     {
    //         std::cout<<"("<<ori.GetContent()<<","<<ori.GetId()<<")";
    //     }
    //     std::cout<<" -->\n";
    //     for(auto replace : iter->second)
    //     {
    //         for(auto data : replace)
    //         {
               
    //             std::cout<<"("<<data.GetContent()<<","<<data.GetId()<<")";
    //         }
    //         std::cout<<std::endl;

    //     }
    // }

    // std::cout << "Before update buffer\n";
    // for (int i = 0; i < this->data_buffer_.size(); i++)
    // {
    //     std::cout << i << " ";
    //     for (int j = 0; j < this->data_buffer_[i].size(); j++)
    //         std::cout << "(" << this->data_buffer_[i][j].GetContent() << "," << this->data_buffer_[i][j].GetId() << ") ";
    //     std::cout << "\n";
    // }

    // @param pos: a vector indicating the position need to move to find the corresponding reduction (caused by reducing previous reduction)
	Symbolic_State<StateType> replacing_state;

	std::vector<std::vector<int> > adjusted_pos;    
    for (int i = 0; i < this->data_buffer_.size(); i++)
        adjusted_pos.push_back(std::vector<int>(this->data_buffer_[i].size(), 0));

    std::vector<Symbolic_Rule<StateType> > rules_in_AOF = this->AOFragment_.GetAllRules();
  
    // for(auto rule : rules_in_AOF)
    // {
    //     if(rule.GetSource().GetId() == 261)
    //     {
    //         std::cerr<<"find 261 in updategraph before upate!\n";
    //         this->graph_.OutputGraph("before_261",PATH,true,true);
    //         for(auto data : this->data_buffer_)
    //         {
    //             for(auto state : data)
    //             {
    //                 std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //             }
    //             std::cerr<<"\n";
    //         }
    //         int a;
    //         std::cin >>a;
    //     }
    // }
    // Add all rules in AOF
    for (Symbolic_Rule<StateType> rule : rules_in_AOF)
        this->graph_.AddRule(rule);

    // if(this->AOFragment_.and_node_.GetId() == 192)
    // {
    //     this->graph_.OutputGraph("just_aof",PATH,false);
    //     std::cerr << "find you!\n";
    //     auto children = this->graph_.ChildrenVertices(453);
    //     for(auto child : children)
    //     {
    //         auto state = this->graph_.GetStateByVertexId(child);
    //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //     }
    //     std::cerr<<"\n";
    //     int a;
    //     std :: cin >> a;
    // }
    // std::cout << "Rules in AOF:\n";
    // for (Symbolic_Rule<StateType> rule : rules_in_AOF)
    // {
    //     // std::cout << "Rule is:\n";
    //     AOG_LIB::Symbolic_State<StateType> source = rule.GetSource();
    //     std::vector<AOG_LIB::Symbolic_State<StateType> > results = rule.GetResults();
    //     // for (AOG_LIB::Symbolic_State<StateType> result : results)
    //     // {
    //     //     std::cout << "\t" << "(" << source.GetContent() << "," << source.GetId() << ")" 
    //     //         << " -> " << "(" << result.GetContent() << "," << result.GetId() << ")" << "\n";
    //     // }
    // }

    // if (this->AOFragment_.or_children_.size() == 2 && this->AOFragment_.or_children_[0].size() == 1 && this->AOFragment_.or_children_[1].size() == 1)
    // {
    //     std::vector<VertexId> children1 = this->graph_.ChildrenVertices(this->graph_.GetVertexIdByState(this->AOFragment_.or_children_[0][0]));
    //     std::vector<VertexId> children2 = this->graph_.ChildrenVertices(this->graph_.GetVertexIdByState(this->AOFragment_.or_children_[1][0]));
    //     std::unordered_set<VertexId> debug(children1.begin(), children1.end());
    //     for (int i = 0; i < children2.size(); i++)
    //     {
    //         if (debug.find(children2[i]) != debug.end())
    //             std::cout << "bug\n";
    //         debug.insert(children2[i]);
    //     }    
    // }
    


    if (rrmap.is_delete_)
    {
        for (Symbolic_Rule<StateType> rule : rules_in_AOF)
        {
            this->graph_.DeleteRule(rule);
            // std::cout << "is delete AOF:\n";
            // AOG_LIB::Symbolic_State<StateType> source = rule.GetSource();
            // std::vector<AOG_LIB::Symbolic_State<StateType> > results = rule.GetResults();
            // for (AOG_LIB::Symbolic_State<StateType> result : results)
            // {
            //     std::cout << "\t" << "(" << source.GetContent() << "," << source.GetId() << ")" 
            //         << " -> " << "(" << result.GetContent() << "," << result.GetId() << ")" << "\n";
            // }
        }
    }
    // this->graph_.OutputGraph("after_delete_aof",PATH,false);
    // Update AOF weight
    std::vector<Symbolic_Rule<StateType> > aof_rules = this->AOFragment_.GetAllRules();
    std::unordered_set<Symbolic_State<StateType> > unique_sources;

    // the first rule is from And to all nodes    
    for (int i = 0; i < aof_rules.size(); i++)
        unique_sources.insert(aof_rules[i].GetSource());

    for (Symbolic_State<StateType> source_state : unique_sources)
    {
        // if this is an Or-node, we need to set weight 
        VertexId source_id = this->graph_.GetVertexIdByState(source_state);
        int pos = -1;
        // record the position of each or children's weight
        std::vector<int> weight_pos;
        if (!this->graph_.GetVertexContent(source_id)->IsAnd())
        {
            std::unordered_map<VertexId, double> weights;
            // all dummy children under the Or-node
            std::vector<VertexId> children_id = this->graph_.ChildrenVertices(source_id);

            // find corresponding Or-children in the AOF
            for (int i = 0; i < this->AOFragment_.or_children_.size(); i++)
            {
                weight_pos.clear();
                bool is_match = true;

                for (VertexId child_id : children_id)
                {
                    // search matching state in AOF or-children
                    int j;
                    for (j = 0; j < this->AOFragment_.or_children_[i].size(); j++)
                    {
                        assert(this->graph_.ChildrenVertices(child_id).size() == 1);
                        if (this->graph_.GetStateByVertexId(this->graph_.ChildrenVertices(child_id)[0]) == this->AOFragment_.or_children_[i][j])
                        {
                            weight_pos.push_back(j);
                            break;
                        }
                    }
                    // unmatched state (child_id)
                    if (j == this->AOFragment_.or_children_[i].size())
                    {
                        is_match = false;
                        break;
                    }
                }
                // find matching or children
                if (is_match)
                {
                    pos = i;
                    break;
                }
            }
            assert(pos != -1);
            for (int i = 0; i < children_id.size(); i++)
                weights.insert(std::make_pair(children_id[i], AOFragment_.weights_[pos][weight_pos[i]]));
            
            // std::cerr<<"all the weights to be updated in first update graph and buffer: \n";
            // for(auto iter = weights.begin(); iter != weights.end(); ++iter)
            // {
            //     std::cerr<<iter->second<<" ";
            //     if(iter->second == 0)
            //     {
            //         std::cerr<<"update weight in update graph and buffer first error!\n";
            //         throw std::exception();
            //     }
            // }
            // std::cerr<<std::endl;
            this->graph_.SetOutEdgeWeights(source_id, weights);
        }
    }

    // this->PrintThirdLevel();


    // Update weights in third level rules
    std::unordered_map<SequenceType<StateType>, int> rule_weight;
    Symbolic_State<StateType> root_state = this->graph_.GetStateByVertexId(this->graph_.GetRoot());
    
    // std::cout << "Rules added in RRMap:\n";
    
    for (auto && it = rrmap.map_.begin(); it != rrmap.map_.end(); it++)
    {
        // std::cout << "From rrmap\n";
        // std::cout << "deleting rules:\n";
        AOG_LIB::Symbolic_State<StateType> source = root_state;
        std::vector<AOG_LIB::Symbolic_State<StateType> > results = it->first;
        // for (AOG_LIB::Symbolic_State<StateType> result : results)
        // {
        //     std::cout << "\t" << "(" << source.GetContent() << "," << source.GetId() << ")" 
        //         << " -> " << "(" << result.GetContent() << "," << result.GetId() << ")" << "\n";
        // }

        // delete root to AOF head
        this->graph_.DeleteRule(Symbolic_Rule<StateType>(root_state, it->first));
    }
        // this->graph_.OutputGraph("after_delete_root_to_aof",PATH,false);
        //delete AOF head to its or_children
        // if (it->first.size() == 1)
        //     for (Symbolic_Rule<StateType> rule : all_rules)
        //         if (rule.GetSource() == it->first[0] && rule.GetResults().size() != 1)
        //         {
        //             this->graph_.DeleteRule(rule);
        //             std::cout << "is delete AOF:\n";
        //             AOG_LIB::Symbolic_State<StateType> source = rule.GetSource();
        //             std::vector<AOG_LIB::Symbolic_State<StateType> > results = rule.GetResults();
        //             for (AOG_LIB::Symbolic_State<StateType> result : results)
        //             {
        //                 std::cout << "\t" << "(" << source.GetContent() << "," << source.GetId() << ")" 
        //                     << " -> " << "(" << result.GetContent() << "," << result.GetId() << ")" << "\n";
        //             }
        //         }

    // get rules that remain in the graph
    // std::cout << "origin_weights.........." << std::endl;
    std::unordered_map<VertexId, double> origin_weights = this->graph_.GetOutEdgeWeights(this->graph_.GetRoot(), false);

    // add root to or_children
    for (auto && it = rrmap.map_.begin(); it != rrmap.map_.end(); it++)
    {
        ConfigBuffer<StateType> new_rules = it->second;
        for (int i = 0; i < new_rules.size(); i++)
        {
            auto rule_weight_it = rule_weight.find(new_rules[i]);
            if (rule_weight_it == rule_weight.end())
            {
                rule_weight.insert(std::make_pair(new_rules[i], 1));
                Symbolic_Rule<StateType> rule(root_state, new_rules[i]);
                this->graph_.AddRule(rule);

                // std::cout << "adding rules:\n";
                // AOG_LIB::Symbolic_State<StateType> source = rule.GetSource();
                // std::vector<AOG_LIB::Symbolic_State<StateType> > results = rule.GetResults();
                // for (AOG_LIB::Symbolic_State<StateType> result : results)
                // {
                //     std::cout << "\t" << "(" << source.GetContent() << "," << source.GetId() << ")" 
                //         << " -> " << "(" << result.GetContent() << "," << result.GetId() << ")" << "\n";
                // }
            }
            else
                rule_weight_it->second++;
        }
        
    }
    // this->graph_.OutputGraph("after_add",PATH,false);

     
    std::vector<VertexId> dummy_children = this->graph_.ChildrenVertices(this->graph_.GetRoot());
    std::unordered_map<VertexId, double> weights = this->graph_.GetOutEdgeWeights(this->graph_.GetRoot(), false);

        // std::cout << "inside updateweights............." << std::endl;
        // for (auto i : weights){
        //     std::cout << "id :" << i.first << " weight: " << i.second << std::endl;
        // }

    for (int i = 0; i < dummy_children.size(); i++)
    {
        std::vector<VertexId> third_level_children = this->graph_.ChildrenVertices(dummy_children[i]);
        SequenceType<StateType> third_level_child;
        for (VertexId child : third_level_children)
            third_level_child.push_back(this->graph_.GetStateByVertexId(child));
        

        auto it = rule_weight.find(third_level_child);
        if (it != rule_weight.end())
        {
            // if this third level sequence already exists in the original sequence, weight is old + new
            if (origin_weights.find(dummy_children[i]) != origin_weights.end())
                weights[dummy_children[i]] = it->second + origin_weights[dummy_children[i]];
            // if it is a new sequence, the weight is just the number of occurence in RRMap
            else
                weights[dummy_children[i]] = it->second;
        }
    }
    
    // std::cerr<<"all the weights to be updated in second update graph and buffer: \n";
    // for(auto iter = weights.begin(); iter != weights.end(); ++iter)
    // {
    //     std::cerr<<iter->second<<" ";
    //     if(iter->second == 0)
    //     {
    //         std::cerr<<"update weight in update graph and buffer second error!\n";
    //         throw std::exception();
    //     }
    // }
    // std::cerr<<std::endl;
    this->graph_.SetOutEdgeWeights(this->graph_.GetRoot(), weights);

    // std::cout << "*****DATA BUFFER IS*****\n";
    // for (int i = 0; i < this->data_buffer_.size(); i++)
    // {
    //     std::cout << "[" << i << "]:";
    //     for (int j = 0; j < this->data_buffer_[i].size(); j++)
    //         std::cout << "(" <<  this->data_buffer_[i][j].GetContent() << "," << this->data_buffer_[i][j].GetId() << ")";
    //     std::cout << "\n";
    // }

    // std::cout << "*****RRMap IS*****\n";
    // for (auto it = rrmap.map_.begin(); it != rrmap.map_.end(); it++)
    // {
    //     for (int i = 0; i < it->first.size(); i++)
    //     {
    //         std::cout << "(" << it->first[i].GetContent() << "," << it->first[i].GetId() << ")";
    //     }
    //     std::cout << "->\n";
    //     for (int i = 0; i < it->second.size(); i++)
    //     {
    //         for (int j = 0; j < it->second[i].size(); j++)
    //             std::cout << "(" << it->second[i][j].GetContent() << "," << it->second[i][j].GetId() << ")";
    //         std::cout << "\n";
    //     }
    // }

    // this->PrintThirdLevel();
    // std::cerr << "Landmark\n";
    
    // Update DataBuffer
    int count = 0;
    for (auto it = rrmap.map_.begin(); it != rrmap.map_.end(); it++)
    {
        // std::cerr << "Inside debug loop\n";
        int substitute_pos = 0;
        int real_pos = 0;
        int total_substitute = it->second.size();
        auto match_it = std::find(this->data_buffer_.begin(), this->data_buffer_.end(), it->first);
        while (match_it != this->data_buffer_.end())
        {
            // std::cerr << "The pos of the sequence in databuffer: " << match_it - this->data_buffer_.begin() << "\n";
            // std::cerr << "The sequence to be replaced in data buffer: \n";
            // for(auto iter = this->data_buffer_[match_it - this->data_buffer_.begin()].begin(); iter != this->data_buffer_[match_it - this->data_buffer_.begin()].end(); iter++)
            // {
            //     std::cerr << "(" << (*iter).GetContent() << ", " << (*iter).GetId() << ") ";
            // }
            // std::cerr << "The size of the substitution container: " << it->second.size() << "\n";
            // std::cerr << "The pos of the substitution: " << substitute_pos << "\n";
            // std::cerr << "The sequence is replaced by: \n";

            // for(auto iter = it->second[substitute_pos].begin(); iter != it->second[substitute_pos].end(); iter++)
            // { 
            //     std::cerr << "(" << (*iter).GetContent() << ", " << (*iter).GetId() << ") ";
            // }
            real_pos = substitute_pos;
            if(substitute_pos >= it->second.size())
            {
                int idx = gen() % it->second.size();
                real_pos = idx;
            }
            this->data_buffer_[match_it - this->data_buffer_.begin()] = it->second[real_pos];
            substitute_pos++;
            match_it = std::find(this->data_buffer_.begin(), this->data_buffer_.end(), it->first);
        }
        count++;
    }
    // std::cout << "After update buffer\n";
    // for (int i = 0; i < this->data_buffer_.size(); i++)
    // {
    //     std::cout << i << " ";
    //     for (int j = 0; j < this->data_buffer_[i].size(); j++)
    //         std::cout << "(" << this->data_buffer_[i][j].GetContent() << "," << this->data_buffer_[i][j].GetId() << ") ";
    //     std::cout << "\n";
    // }

    // this->graph_.OutputGraph("after_bug",PATH,false);
    
    // // Add weights in AOF
    // std::vector<VertexId> or_children;
    // // if source is or node
    // if (!this->graph_.GetVertexContent(this->graph_.GetVertexIdByState(replacing_state))->IsAnd())
    // {
    //     or_children.push_back(this->graph_.GetVertexIdByState(replacing_state));
    //     // std::cerr<<"Root is an or-node\n";
    // }
    // else
    // {
    //     std::vector<VertexId> children = this->graph_.ChildrenVertices(this->graph_.GetVertexIdByState(replacing_state));
    //     // if children is or node
    //     for (int i = 0; i < children.size(); i++)
    //         if (!this->graph_.GetVertexContent(children[i])->IsAnd())
    //         {
    //             or_children.push_back(children[i]);
    //             // std::cerr<<"children is or node!\n";
    //         }
    // }

    // for (int i = 0; i < or_children.size(); i++)
    // {
    //     std::unordered_map<VertexId, double> weights;
    //     std::vector<VertexId> and_under_or = this->graph_.ChildrenVertices(or_children[i]);
        
    //     // if size is 1, then this is an and node, no processing
    //     if (and_under_or.size() > 1 && !this->graph_.GetVertexContent(or_children[i])->IsAnd())
    //     {

    //         // find corresponding aof or children position in AOF, because it can mismatch with position in TAOG
    //         int aof_pos = 0;
    //         for (; aof_pos < this->AOFragment_.or_children_.size(); aof_pos++)
    //         {
    //             auto it = and_under_or.begin();
    //             for (; it != and_under_or.end(); it++)
    //             {
    //                 assert(this->graph_.ChildrenVertices(*it).size() == 1);
    //                 Symbolic_State<StateType> state = this->graph_.GetStateByVertexId(this->graph_.ChildrenVertices(*it)[0]);
    //                 // if cannot find one of the state, means this is not the AOF.or_children
    //                 if (std::find(this->AOFragment_.or_children_[aof_pos].begin(), 
    //                     this->AOFragment_.or_children_[aof_pos].end(), state) == this->AOFragment_.or_children_[aof_pos].end())
    //                     break;
    //             }
    //             // found AOF.or_children that matches all
    //             if (it == and_under_or.end())
    //                 break;
    //         }

    //         // case when no matching AOF.or_children is found, which is impossible
    //         if (aof_pos == this->AOFragment_.or_children_.size())
    //             continue;   

    //         for (int j = 0; j < and_under_or.size(); j++)
    //         {
    //             Symbolic_State<StateType> leaf_state = this->graph_.
    //                         GetStateByVertexId(this->graph_.ChildrenVertices(and_under_or[j])[0]);

    //             for (int pos = 0; pos < this->AOFragment_.or_children_[aof_pos].size(); pos++)
    //             {
    //                 if (this->AOFragment_.or_children_[aof_pos][pos] == leaf_state)
    //                 {
    //                     weights.insert(std::make_pair(and_under_or[j], this->AOFragment_.weights_[aof_pos][pos]));
    //                     break;
    //                 }
    //             }
    //         }
    //         this->graph_.SetOutEdgeWeights(or_children[i], weights);
    //     }
    // }

    /** std::cerr <<"after change rules:\n";
    std::vector<Symbolic_Rule<StateType> > all_rules = this->graph_.GetRules();
    for (Symbolic_Rule<StateType> rule : all_rules)
    {
        std::cerr << "From " << rule.GetSource().GetContent() << " To ";
        for (Symbolic_State<StateType> result : rule.GetResults())
        {
            std::cerr << result.GetContent() << " ";
        }
        std::cerr << "\n";
    } 

    std::cerr << "finish updategraph\n"; */
    // unsigned root_id = this->graph_.GetRoot();
    // auto root_state = this->graph_.GetStateByVertexId(root_id);
    // assert(root_state.GetId() == 0 );
    /** std::cerr<<"root id: "<<root_state.GetId()<<std::endl; */

    // std::cerr<<"Outside UpdateGraph()\n";
    // if(this -> AOFragment_.and_node_.GetId() == 5413)
    // {
    //     std::cerr<<"after 5413\n";
    //     std::cerr<<"size of or_children: "<<this->AOFragment_.or_children_.size()<<  std::endl;
    //     this->graph_.OutputGraph("visualize_5413_after",PATH,true);
    //     // int a;
    //     // std::cin >>a;

    // }
    // if(this -> AOFragment_.and_node_.GetId() == 24)
    // {
    //     std::cerr<<"after 24\n";
    //     std::cerr<<"size of or_children: " << this->AOFragment_.or_children_.size()<<  std::endl;
    //     this->graph_.OutputGraph("visualize_24_after",PATH,true);
    //     // int a;
    //     // std::cin >>a;
    // // }
    // if(this->AOFragment_.and_node_.GetId() == 3207)
    // {
    //     std::cout << "*****DATA BUFFER IS*****\n";
    //     for (int i = 0; i < this->data_buffer_.size(); i++)
    //     {
    //         std::cout << "[" << i << "]:";
    //         for (int j = 0; j < this->data_buffer_[i].size(); j++)
    //             std::cout << "(" <<  this->data_buffer_[i][j].GetContent() << "," << this->data_buffer_[i][j].GetId() << ")";
    //         std::cout << "\n";
    //     }
    //     std::cerr<<"find you as and node after update databuffer!\n";
    //     std::cerr<<"is_delete: "<<replacement.rrmap_.is_delete_<<std::endl;                   
    //     this->graph_.OutputGraph("add3207_after",PATH,false);
    //     std::cerr<<"and node: ("<<this->AOFragment_.and_node_.GetContent()<<", "<<this->AOFragment_.and_node_.GetId()<<")"<<std::endl;                
    //     for(auto nodes: this->AOFragment_.or_children_)
    //     {
    //         for(auto node_state : nodes)
    //         {
    //             std::cerr<<"("<<node_state.GetContent()<<", "<<node_state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
                
    //     int a;
    //     std::cin >> a;
    // }
    // for(auto or_node : this->AOFragment_.or_children_)
    // {
    //     for(auto state : or_node)
    //     {
    //         if(state.GetId() == 3207)
    //         {
    //             std::cerr<<"find you from or children!\n"; 
    //             std::cerr<<"is_delete: "<<replacement.rrmap_.is_delete_<<std::endl;                                              
    //             this->graph_.OutputGraph("add3207_after",PATH,false);
    //             std::cerr<<"and node: ("<<this->AOFragment_.and_node_.GetContent()<<", "<<this->AOFragment_.and_node_.GetId()<<")"<<std::endl;                
    //             for(auto nodes: this->AOFragment_.or_children_)
    //             {
    //                 for(auto node_state : nodes)
    //                 {
    //                     std::cerr<<"("<<node_state.GetContent()<<", "<<node_state.GetId()<<")";
    //                 }
    //                 std::cerr<<std::endl;
    //             }
                
    //             int a;
    //             std :: cin >> a;

    //         }
    //     }
    // }

    // if(this->graph_.IsValidVertex(115) && this->graph_.ChildrenVertices(115).size() == 1 && this->graph_.ChildrenVertices(115)[0] == 91)
    // {
    //     std::cerr<<"find the bad bug!!\n";
    //     auto node_115 = this->graph_.GetStateByVertexId(115);
    //     int node_id = node_115.GetId();
    //     std::cerr<<"vertice 115's state id: "<<node_id<<std::endl;
    //     auto children_vertices = this->graph_.ChildrenVertices(115);
    //     auto child_state = this->graph_.GetStateByVertexId(children_vertices[0]);
    //     // auto test = children_vertices[1];
    //     // std::cerr<<"the unknown one: "<<test<<std::endl;
    //     std::cerr<<"the only child: "<<children_vertices[0]<<std::endl;
    //     std::cerr<<"the only child's state id: "<<child_state.GetId()<<std::endl;
    //     this->graph_.OutputGraph("peek",PATH,false);
    //     int a;
    //     std::cin >> a;
    // }

    // //check if data buffer and aog are consistent
    // auto third_level = this -> GetThirdLevel();
    // for(auto data : this->data_buffer_)
    // {
    //     if (third_level.find(data) == third_level.end())
    //     {
    //         std::cerr<<"find you!\n";
    //         this->graph_.OutputGraph("buggy_graph",PATH,false);
    //         int a;
    //         std::cin >> a;

    //     }
    // }
    // if(this->AOFragment_ == test_aof)    
    // {
    //     std::cerr<<"after 693\n";
    //     std::cerr<<"data buffer: \n";
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     this->graph_.OutputGraph("after_add_nose",PATH,true);
        
    // }
    // else
    // {
    //     std::cerr << "this round's data buffer:\n";
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    // }
    // if(this->graph_.IsValidVertex(470))
    // {
    //     std::cerr <<"find you1\n";
    //     this->graph_.OutputGraph("find470_false",PATH,false);
    //     this->graph_.OutputGraph("find470_true",PATH,true);
    //      std::cerr << "this round's data buffer:\n";
    //     for(auto data : this->data_buffer_)
    //     {
    //         for(auto state : data)
    //         {
    //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     int a;
    //     std::cin >> a;
    // }
    // if(this->graph_.IsValidVertex(470))
    // {
        
    //     auto children_vtxs = this->graph_.ChildrenVertices(470);
    //     if(children_vtxs.size() == 3)
    //     {
    //         std::cerr<< "find 3 children\n";
            
    //         for(auto child : children_vtxs)
    //         {
    //             auto state = this->graph_.GetStateByVertexId(child);
    //             std :: cerr << "("<<state.GetContent()<<", "<<state.GetId()<<")";
    //             std::cerr << std::endl;
    //         }
            
    //         int b;
    //         std::cin >> b;
    //         for(auto child : children_vtxs)
    //         {
    //             auto state = this->graph_.GetStateByVertexId(child);
    //             if(state.GetContent() == "Bumpy")
    //             {
    //                 this->graph_.OutputGraph("find_bumpy",PATH,false);
    //                 int a;
    //                 std::cin >>a;

    //             }
    //         }
    //     }

    // }
    
    // SequenceType<std::string> test_seq;
    // test_seq.push_back(Symbolic_State<std::string>(652));
    // // test_seq.push_back(Symbolic_State<std::string>("Bumpy",true));
    // test_seq.push_back(Symbolic_State<std::string>(693));


    // for(auto data : this->data_buffer_)
    // {
    //     if(data == test_seq)
    //     {
    //         std::cerr<<"current AOF:\n";
    //         std::cerr<<"and node: ("<<this->AOFragment_.and_node_.GetContent()<<", "<<this->AOFragment_.and_node_.GetId()<<")\n";
    //         for(auto or_node : this->AOFragment_.or_children_)
    //         {
    //             for(auto state : or_node)
    //             {
    //                 std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //             }
    //             std::cerr<<std::endl;
    //         }

    //         this->graph_.OutputGraph("problem",PATH,false);
    //         int a;
    //         std::cin >> a;
    //     }
    // }
    // if(this->graph_.IsValidVertex(453) && this->graph_.ChildrenVertices(453).size() == 2)
    // {
    //     std::cerr << "find you!\n";
    //     auto children = this->graph_.ChildrenVertices(453);
    //     for(auto child : children)
    //     {
    //         auto state = this->graph_.GetStateByVertexId(child);
    //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //     }
    //     std::cerr<<"\n";
    //     std::cerr <<"all rules: \n";
    //     for(auto rule : aof_rules)
    //     {
    //         auto source = rule.GetSource();
    //         auto targets = rule.GetResults();
    //         std::cerr<<"("<<source.GetContent()<<", "<<source.GetId()<<") -> ";
    //         for(auto target : targets)
    //         {
    //             std::cerr<<"("<<target.GetContent()<<", "<<target.GetId()<<")  ";
                
    //         }
    //         std::cerr<<std::endl;
    //     }
    //     this->graph_.OutputGraph("wrong_order",PATH,false);
    //     int a;
    //     std :: cin >>a;
    // }
    // std::cerr<<"arrive here!\n";
    // for(auto rule : rules_in_AOF)
    // {
        // if(rule.GetSource().GetId() == 261)
        // {
        //     std::cerr<<"find 261 in updategraph after update!\n";
        //     this->graph_.OutputGraph("after_261",PATH,true,false);
        //     for(auto data : this->data_buffer_)
        //     {
        //         for(auto state : data)
        //         {
        //             std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
        //         }
        //         std::cerr<<"\n";
        //     }
        //     int a;
        //     std::cin >>a;
        // }
    // }
    // if(this->graph_.IsValidVertex(452) && this->graph_.GetStateByVertexId(452).GetId() == -1)
    // {
    //     std::cerr<<"find you!\n";
    //     this->graph_.OutputGraph("452_-1",PATH,true,SIMPLIFY);
    //     int a;
    //     std::cin >>a;
    // }
    // if(this->graph_.IsValidVertex(453))
    // {
    //     std::cerr<<"find 453\n";
    //     this->graph_.OutputGraph("453",PATH,true,SIMPLIFY);
    //     if(this->graph_.IsValidVertex(452) )
    //     {
    //         std::cerr<<"452 valid inside 453!\n";
    //         this->graph_.OutputGraph("452_inside_453",PATH,true,SIMPLIFY);
    //         if(this->graph_.GetStateByVertexId(452).GetId() == -1)
    //         {
    //             std::cerr<<"find 452 to be -1!\n";
    //             this->graph_.OutputGraph("452_-1_453",PATH,true,SIMPLIFY);
    //         }
    //     }
    //     int a;
    //     std::cin >>a;
    // }
    // if(this->AOFragment_.and_node_.GetId() == 18760)
    // {
    //     std::cerr<<"after update 18760\n";
    //     this->graph_.OutputGraph("after_18760",PATH,true,SIMPLIFY);
    //     if(this->graph_.IsValidVertex(452) && this->graph_.GetStateByVertexId(452).GetId() == -1)
    //     {
    //         std::cerr<<"hahaaaaaaaaaaaaaaaaaaa!\n";
    //     }
    //     int a;
    //     std::cin >>a;
    // }
    // this ->graph_.OutputGraph("after_update_graph",PATH,true,SIMPLIFY);
    // std::cerr<<"out of update graph and buffer!\n";
    return ;
}

template <class StateType>
int Online_Learner<StateType>::Checkpoint::CountGraphReplacement(const RdBuffer<StateType> &reductions)
{
    std::unordered_set<std::pair<SequenceType<StateType>, int > > unique_rd;
    for(const auto & rd : reductions)
        for(const auto & pos : rd.second)
            unique_rd.insert(std::pair<SequenceType<StateType>, int>(this->data_buffer_[pos.first], pos.second));

    return unique_rd.size();
        
}

//Note this function calculates the log posterior gain rather than the actual posterior gain
template <class StateType>
double Online_Learner<StateType>::Checkpoint::CalculatePosterior(const AOFStruct<StateType> &AOF,
                                                                      const RdBuffer<StateType> &reduction_full, bool use_beta, bool is_delete)
{
    //Need: # of Or-Nodes in the AOF and the # of child nodes under each of the or nodes
    // std::cerr << "Inside Posterior: " << std::endl;
    double factor1_log = 0, factor2_log = 0;
    int sizeOfReductions = 0;
    for (auto && it = reduction_full.begin(); it != reduction_full.end(); it++)
        sizeOfReductions += it->second.size();

    // if(reduction_full.empty())
    //     std::cerr<<"reduction full inside calculateposterior is empty!\n";
    // if (sizeOfReductions == 0)
    // {
    //     // int a;
    //     // std::cerr << "size of reduction is zero!" << std::endl;
    //     // cin >> a;
    //     return std::numeric_limits<double>::lowest();
    // }

    //return negative infinity if AOF is a unigram with only one children
    if(AOF.or_children_.size() == 1 && AOF.or_children_[0].size() == 1)
    {
        // std::cerr<<"bad AOF!\n";
        return std::numeric_limits<double>::lowest();
    }

    if(reduction_full.empty() && !is_delete)
    {
        // std::cerr<<"empty reduction while adding an aof!\n";
        return std::numeric_limits<double> :: lowest();
    }
    int orChildren = 0;
    //Factor 1:
    int n = AOF.or_children_.size();

    // if (n == 1)
    // {
    //     Symbolic_State<StateType> c2b(2860);
    //     Symbolic_State<StateType> c2a(3062);
    //     if ((AOF.or_children_[0][0] == c2a && AOF.or_children_[0][1] == c2b) || (AOF.or_children_[0][0] == c2b && AOF.or_children_[0][1] == c2a))
    //     {
    //         std::cerr << "found\n";
    //     }
    // }

    //This is old factor 1
    for (int i = 0; i < n; i++)
    {
        orChildren += AOF.or_children_[i].size();
        // if (AOF.or_children_[i].size() == 1)
        //     orChildren += 1;
        // else
        //     orChildren += 2*AOF.or_children_[i].size();
        for (int j = 0; j < AOF.or_children_[i].size(); j++)
        {   
            //find how many time node a_ij appears in i-th position of all reductions
            int RD_i = 0;
            for(auto && it = reduction_full.begin(); it != reduction_full.end(); it++)
                if(it->first[i] == AOF.or_children_[i][j])
                    RD_i += it->second.size();

            if (RD_i != 0)
                factor1_log += RD_i * log(RD_i);

        }
    }

    double factor1 = 0;
    

    // factor1 = factor1_log - n * sizeOfReductions * log(sizeOfReductions);
    factor1 = CalculateNewFactor1(AOF,reduction_full);
    assert(factor1 == 1 ? (reduction_full.size() == 1) : true);

    //Calculate Factor2 
    // iterate CM matrix by context
    int deleted_node = 0;
    int reduced_node = 0;
    double prior;

    auto & context_index = boost::multi_index::get<0>(this->CM_);
    std::unordered_set<ContextType<StateType>> all_contexts;

    if(OFFLINE)
    {
        // std::cerr<<"using offline!\n";
        // int a;
        // std::cin >> a;
       for(auto && it = context_index.begin(); it != context_index.end();++it)
        {
            //if the context is not visited before
            if(all_contexts.find(it->context_) == all_contexts.end())
            {
                std::vector<int> counts;
                all_contexts.insert(it->context_);
                auto same_contxt = context_index.equal_range(it->context_);
                int merged_count = 0;
                for(auto && iter = same_contxt.first; iter != same_contxt.second; ++iter)
                {
                    assert(iter->count_ != 0);
                    counts.push_back(iter->count_);
                    ++merged_count;
                }
                --merged_count;

                //total number of CM entries with same contexts
                double sum = std::accumulate(counts.begin(),counts.end(),0);
                for (int i = 0; i < counts.size(); i++)
                    factor2_log += counts[i]*log(sum/counts[i]);
                

                //calculate the context size            
                int total_context_size = 0;
                for(int i = 0; i < it->context_.size(); ++i)
                    total_context_size += it->context_[i].size();
                
                //calculate the size of all configurations that can be replaced
                int total_config_size = it->configuration_.size() * n;
                deleted_node += merged_count * (total_context_size + total_config_size + 1);
                
                reduced_node += it->configuration_.size() * (n-1);
            }
        }

    }
    
    else
    {
        for(auto && it = context_index.begin(); it != context_index.end();++it)
        {
            //if the context is not visited before
            if(all_contexts.find(it->context_) == all_contexts.end())
            {
                std::vector<int> counts;
                all_contexts.insert(it->context_);
                auto same_contxt = context_index.equal_range(it->context_);
                for(auto && iter = same_contxt.first; iter != same_contxt.second; ++iter)
                {
                    assert(iter->count_ != 0);
                    counts.push_back(iter->count_);
                }

                //total number of CM entries with same contexts
                double sum = std::accumulate(counts.begin(),counts.end(),0);
                for (int i = 0; i < counts.size(); i++)
                    factor2_log += counts[i]*log(sum/counts[i]);
            }
        }


         //calculate reduced node and deleted node using the CM generated on thirdlevel
        //prior:    
        

        //save the true CM
        Context_Matrix<StateType> true_CM = this->CM_;

        //Gnerate CM according to third level
        std::unordered_map<std::vector<Symbolic_State<StateType> >, double> third_level = this->GetThirdLevel();
        ConfigBuffer<StateType> third_level_seq;
        for(auto && iter = third_level.begin(); iter != third_level.end(); ++iter)
            third_level_seq.push_back(iter->first);

        assert(!AOF.sorted_configs_.empty());
        GenerateCMAndReductions(AOF.sorted_configs_, third_level_seq,true);

        //Calculate reduced and deleted node according to the new CM
        auto& third_level_context_index = boost::multi_index::get<0>(this->CM_);
        std::unordered_set<ContextType<StateType>> third_level_all_contexts;

        for(auto && it = third_level_context_index.begin(); it != third_level_context_index.end();++it)
        {
            //if the context is not visited before
            if(third_level_all_contexts.find(it->context_) == third_level_all_contexts.end())
            {
                third_level_all_contexts.insert(it->context_);
                auto same_contxt = third_level_context_index.equal_range(it->context_);
                int merged_count = 0;
                for(auto && iter = same_contxt.first; iter != same_contxt.second; ++iter)
                {
                    assert(iter->count_ != 0);
                    ++merged_count;
                }
                --merged_count;

                //calculate the context size            
                int total_context_size = 0;
                for(int i = 0; i < it->context_.size(); ++i)
                    total_context_size += it->context_[i].size();
                
                //calculate the size of all configurations that can be replaced
                int total_config_size = it->configuration_.size() * n;
                deleted_node += merged_count * (total_context_size + total_config_size + 1);
                
                reduced_node += it->configuration_.size() * (n-1);
            }
        }


        //restore the original CM
        this->CM_ = true_CM;


    }
    

    double factor2 = factor2_log;
    double likelihood = ALPHA_LIKELIHOOD * (factor1 + factor2);

   


    if(use_beta)
        prior = -ALPHA_PRIOR * (AOF.AOFGrammarSize() * BETA  - deleted_node * GAMMA - reduced_node * REDUCED);
    else
        prior = -ALPHA_PRIOR * (AOF.AOFGrammarSize(1) - deleted_node - reduced_node);
        // std::cerr<<"prior gain returned for Generate(): "<<prior<<std::endl;
        // std::cerr<<"log likelihood gain returned for Generate(): "<< likelihood<<std::endl;
        
        // std::cerr<<"Posterior gain for testing: "<< likelihood+prior<<std::endl<<"prior: "<< prior<< "\nGrammarSize: " << AOF.AOFGrammarSize()
        // <<"AOFGrammar size: "<<AOF.AOFGrammarSize()<<"\n"
        // <<"deleted_node: "<<deleted_node<<"\n"
        // <<"reduced_node: "<<reduced_node<<"\n";
        // std::cerr <<"factor1: "<<factor1<< "\n" <<"factor2: "<<factor2<<std::endl<<std::endl
        // <<"n: "<<n<<"\n"
        // <<"context size: "<<all_contexts.size()<<std::endl;

        // std::cerr<<"AOF or_children: \n";
        // for(auto or_node : AOF.or_children_)
        // {
        //     for(auto state : or_node)
        //     {
        //         std::cerr<<state.GetContent()<<"_";                
        //     }
        //     std::cerr<<std::endl;
        // }

        // std::cerr << "AOF weights: \n";
        // for(auto or_node : AOF.weights_)
        // {
        //     for(auto weight : or_node)
        //         std::cerr <<weight << "_";
        //     std::cerr<<std::endl;
        // }

        // std::cerr << "ContextMatrix Size is" << CM_.size() << std::endl; 
        // auto & context_index = boost::multi_index::get<0>(CM_);
        // std::cerr<<"Context Size: "<<context_index.size()<<std::endl;
        // std::unordered_set<ContextType<StateType> > explored_contexts;
        // int i = 1;
        // for (auto it = context_index.begin(); it != context_index.end(); it++)
        // {
        //     if(explored_contexts.find(it->context_) == explored_contexts.end())
        //     {
        //         explored_contexts.insert(it->context_);
        //         auto same_context = context_index.equal_range(it->context_);
        //         for(auto iter = same_context.first; iter != same_context.second;++iter)
        //         {
        //             std::cerr<<i<<std::endl;
        //             ++i;
        //             std::cerr << "Config:";
        //             for (int i = 0; i < iter->configuration_.size(); i++)
        //             {
        //                 std::cerr << "("<<iter->configuration_[i].GetContent() << " ,"<<iter->configuration_[i].GetId()<<") ";
                        
        //             }
        //             std::cerr << "|| Context:";

        //             for (int j = 0; j < iter->context_.first.size(); j++)
        //             {
        //                 std::cerr << "("<<iter->context_.first[j].GetContent() << " ,"<<iter->context_.first[j].GetId() << ") ";
        //             }
        //             std::cerr << "******";
        //             for (int j = 0; j < iter->context_.second.size(); j++)
        //             {                
        //                 std::cerr << "("<<iter->context_.second[j].GetContent() << " ,"<<iter->context_.second[j].GetId() << ") ";
        //             }
        //             std::cerr<<"count: "<<iter->count_;
        //             std::cerr << "\n";

        //             std::cerr << std::endl;
        //         }
        //     }
            
            // auto same_context = context_index.equal_range(it->context_);
            // for(auto iter = same_context.first; iter != same_context.second;++iter)
            // {
            //     std::cerr << "Config:";
            //     for (int i = 0; i < iter->configuration_.size(); i++)
            //     {
            //         std::cerr << "("<<iter->configuration_[i].GetContent() << " ,"<<iter->configuration_[i].GetId()<<") ";
                   
            //     }
            //     std::cerr << "|| Context:";

            //     for (int j = 0; j < iter->context_.first.size(); j++)
            //     {
            //         std::cerr << "("<<iter->context_.first[j].GetContent() << " ,"<<iter->context_.first[j].GetId() << ") ";
            //     }
            //     std::cerr << "******";
            //     for (int j = 0; j < iter->context_.second.size(); j++)
            //     {                
            //         std::cerr << "("<<iter->context_.second[j].GetContent() << " ,"<<iter->context_.second[j].GetId() << ") ";
            //     }
            //     std::cerr<<"count: "<<iter->count_;
            //     std::cerr << "\n";

            //     std::cerr << std::endl;
            // }
            
    //     }
    // }

    // if (factor2 != 0 && n == 1)
    //     factor2 = factor2;
    
    // int a;
    // if(this->posterior_ + likelihood + prior > 0)
    // {
    //     std::cerr<<"Posterior gain for testing: "<< likelihood+prior<<std::endl<<"prior: "<< prior<< "\nGrammarSize: " << AOF.AOFGrammarSize()
    //     <<"AOFGrammar size: "<<AOF.AOFGrammarSize()<<"\n"
    //     <<"deleted_node: "<<deleted_node<<"\n"
    //     <<"reduced_node: "<<reduced_node<<"\n";
    //     std::cerr <<"factor1: "<<factor1<< "\n" <<"factor2: "<<factor2<<std::endl<<std::endl
    //     <<"n: "<<n<<"\n"
    //     <<"context size: "<<all_contexts.size()<<std::endl;

    //     std::cerr<<"AOF or_children: \n";
    //     for(auto or_node : AOF.or_children_)
    //     {
    //         for(auto state : or_node)
    //         {
    //             std::cerr<<state.GetContent()<<"_";                
    //         }
    //         std::cerr<<std::endl;
    //     }

    //     std::cerr << "AOF weights: \n";
    //     for(auto or_node : AOF.weights_)
    //     {
    //         for(auto weight : or_node)
    //             std::cerr <<weight << "_";
    //         std::cerr<<std::endl;
    //     }
    // }
    // std::cerr << "Size of Reducitons: " << sizeOfReductions << std::endl;
    // for (auto it = reduction_full.begin(); it != reduction_full.end(); it++)
    // {
        
    //     std::cerr << "Reduciton: ";
    //     for (auto iter = (*it).first.begin(); iter != (*it).first.end(); iter++)
    //     {
    //         std::cerr << iter->GetContent() << "_";
    //     }
    //     std::cerr << std::endl;

        //print positions:
        // for (auto iter = (*it).second.begin(); iter != (*it).second.end(); iter++)
        // {
        //     std::cerr << "Pos: "
        //                 << "(" << (*iter).first << "), (" << (*iter).second << ")" << std::endl;
        // }
    // }
    // std::cerr << "before cin\n";
    // std::cin >> a;
    
    
    // if(likelihood + prior == -10.6)
    // {
    //     std::cerr <<"the current AOF: \n";
    //     for(auto or_node : AOF.or_children_)
    //     {
    //         for(auto state : or_node)
    //             std::cerr <<"("<<state.GetContent() << ", "<<state.GetId() << ")";
    //         std::cerr << std::endl;
    //     }
    // // }
    if(!use_beta)
    {
        std::cerr<<"Posterior gain for testing: "<< likelihood+prior<<std::endl<<"prior: "<< prior<< "\nAOF GrammarSize: " << AOF.AOFGrammarSize(1)
        << "\nDeleted_node: " << deleted_node << "\nReduced_node: " << reduced_node << std::endl
        << "factor1: "<<factor1<< "\n" <<"factor2: "<<factor2<<std::endl
        <<"n: "<<n<<"\n"
        <<"context size: "<<all_contexts.size()<<std::endl<<std::endl;
    }
    

    return likelihood+prior;

}


template <class StateType>
double Online_Learner<StateType>::Checkpoint::CalculateNewFactor1(const AOFStruct<StateType> &AOF,
                                                                      const RdBuffer<StateType> &reduction_full)
{
    //normalize weight 
    std::vector<std::vector<double>> probs;
    for(const auto & or_node : AOF.weights_)
        probs.push_back(std::vector<double>(or_node.size(),0.0));
    
    for(int i =0; i < AOF.weights_.size(); ++i)
    {
        double sum_of_elems = std::accumulate(AOF.weights_[i].begin(), AOF.weights_[i].end(), 0);
        for(int j =0; j < AOF.weights_[i].size();++j)
        {
            probs[i][j] = AOF.weights_[i][j]/ sum_of_elems;
            // if(probs[i][j] == 0)
            // {
            //     std::cerr<<"probs is  0!\n";
            //     throw std::exception();
            // }
        }
    }

    
    double total_factor1_log = 0;
    for(const auto & reduction : reduction_full)
    {
        double factor1_log = 0;
        int reduction_size = reduction.second.size();
        SequenceType<StateType> config = reduction.first;
        for(int i = 0; i < config.size();++i)
        {
            for(int j = 0; j < AOF.or_children_[i].size();++j)
            {
                if(config[i] == AOF.or_children_[i][j])
                {
                    factor1_log += log(probs[i][j]);
                    break;
                }
            }
        }
        factor1_log *= reduction_size;
        total_factor1_log += factor1_log;

    }
    
    assert(total_factor1_log <= 0);
    return total_factor1_log;
}                                                                                                                               
                                                                      

template<class StateType>
SN_Ptr Online_Learner<StateType>::Checkpoint::Select(double Cp, bool skip_expand)
{
    // If the graph reaches certain threshold complexity
    // we consider it as terminal node, return nullptr
    // TODO: deal with terminal node
	/** std::cerr << "MCTS Search: inside Select" << std::endl; */

    // If it is not a terminal node, and does not have a child
    //      Expand itself
    // If it is not a termnial node, and has children
    //      with PROBABILITY% we expand itself
    //      with 1-PROBABILITY% we find children with max UCT

    // Don't have children, expand itself
    std::cerr << "[Select()]\n";

    std::vector<SN_Ptr> children = this->GetChildren();
    if (children.empty())
    {
    	std::cerr << "\tInside Select:no children, return itself\n"; 
        return this->shared_from_this();
    }

    // std::cerr<<"selected node's AOF: "<<std::endl;
    // for(auto or_node : this->AOFragment_.or_children_)
    // {
    //     for(auto state : or_node )
    //     {
    //         auto content = state.GetContent();
    //         if(content == "")
    //             std::cout << state.GetId()<<"_";
    //         else
    //             std::cout << content<<"_";
            
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;

    // Have children, PROBABILITY% expand itself
    double weight = EPSILON_MIN + (EPSILON_START - EPSILON_MIN) * exp(-1.0 * this->GetChildren().size() / EPSILON_DECAY);
    std::cout << "weight: "<<weight<<std::endl;
    // bool do_expand = rand()*1.0/RAND_MAX < weight;
    double prob = gen()/(double)(gen.max() - gen.min());
    std::cout <<"expand prob" << prob<<std::endl;
    // int a;
    // std::cin >> a;
    bool do_expand = prob  < weight;

    if (!skip_expand && do_expand)
    {
        // std:cerr<<"selected checkpoint's AOF and_node id:" << this->AOFragment_.and_node_.GetId()<<std::endl;
        return this->shared_from_this();

    }

    // Have children, 1-PROBABILITY% find max UCT
    double max_score = children[0]->GetValue();
    SN_Ptr max_node = children[0];

    int pos = 0;
    int max_pos = 0;
    for (const auto &child : children)
    {
        // Every presented child should have at least VisitCount 1
        if (!child->GetVisitCount())
        {
            std::cerr << "\tSelect() find child with 0 VisitCount\n";
            throw std::exception();
        }
        else
        {
            unsigned N_parent = this->GetVisitCount();
            unsigned N_child = child->GetVisitCount();
            double Q_child = child->GetValue();
            double score = Q_child + Cp * sqrt(2 * log(N_parent) / N_child);
            std::cerr << "\tAt level ["<< this->GetLevel()  << "], Children [" << pos << "] score: " << score << "\n";
	    std::cerr << "Q value: "<<Q_child << "  visit count term: "<< Cp * sqrt(2 * log(N_parent) / N_child)<< "\n";
            // find the child with max UCT
            if (score > max_score)
            {
                max_score = score;
                max_node = child;
                max_pos = pos;   
            }
        }
        pos++;
    }
    std::cerr << "\tSelected #" << max_pos << " children\n";
    auto max_cp_node = static_cast<Checkpoint*> (max_node.get());
    std::cerr<<"The selected node's AOF:\n";
    std::cerr<<"and node: "<< "("<<max_cp_node -> AOFragment_.and_node_.GetContent()<<", "<< max_cp_node -> AOFragment_.and_node_.GetId()<<")\n";
    for(auto or_node : max_cp_node->AOFragment_.or_children_)
      {
    	for(auto state : or_node)
    	  {
    	    std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    	  }
    	std::cerr<<std::endl;
      }
    return max_node->Select(Cp,false);
}


template<class StateType>
SN_Ptr Online_Learner<StateType>::Checkpoint::Expand()
{
    
	/** std::cerr << "MCTS Search: Expand:" << std::endl; */
    
    // Create a copy of current node
    // Connect current node as a children of the new node
    // Delete current node from its parent
    std::cerr << "[Expand()]\n";
    int num_of_children = this->GetChildren().size();
    // std::cerr << "Number of Children: " << num_of_children << "\n";

    // std::cerr<<"before expand, parent's AOF:\n";
    // for(auto or_node : this->AOFragment_.or_children_)
    // {
    //     for(auto state : or_node)
    //     {
    //         std::cerr<<state.GetContent()<<"_";
    //     }
    //     std::cerr<<std::endl;
    // }

    std::shared_ptr<Checkpoint> expanded_node = std::make_shared<Checkpoint>(*this);
    //clear the node's content
    expanded_node->Clean();
    expanded_node->SetLevel(this->GetLevel()+1);    
    expanded_node->SetParent(this->shared_from_this());
    expanded_node -> num_of_neg_post_gain_ = 0;
    // If can further increase posterior, do expand
    // std::cout << "Before generate fragment\n";
    bool success;
    double delete_prob = gen()/(double)(gen.max() - gen.min());
    // std::cerr<<"delete prob in expand():"<<delete_prob<<std::endl;
    // int c;
    // std ::cin >> c;
    if(delete_prob <= DELETE_PROB)
    {
        // std::cerr<<"prepare to delete in expand()\n";
        success = expanded_node->Delete();
        if(!success)
        {
            std::cerr<<"delete failed, prepare to generate()\n";
            success = expanded_node->Generate(EXPAND_AOF_ITER, EXPAND_WITH_STOCHASTIC);
            

        }
    }
    else  
    {
       
        success = expanded_node->Generate(EXPAND_AOF_ITER, EXPAND_WITH_STOCHASTIC);

    }
    // std::cerr<<"after expand, parent's AOF:\n";
    // for(auto or_node : this->AOFragment_.or_children_)
    // {
    //     for(auto state : or_node)
    //     {
    //         std::cerr<<state.GetContent()<<"_";
    //     }
    //     std::cerr<<std::endl;
    // }
    // std::cerr<<std::endl;
   

    //  std::cerr<<"after expand, expanded node's parent's AOF:\n";
    // for(auto or_node : static_cast<Checkpoint *>(expanded_node->GetParent().get())->AOFragment_.or_children_)
    // {
    //     for(auto state : or_node)
    //     {
    //         std::cerr<<state.GetContent()<<"_";
    //     }
    //     std::cerr<<std::endl;
    // }
    // std::cerr<<"Expand finds AOF\n";
    /** std::cerr<<"out of generateFragment!\n"; */
    if (success)
    {
        /** std::cerr << "MCTS Search: Expand: Before update" << std::endl;
        // update expanded node's private parameters*/
        // std::cerr<<"Before updateGraphandBuffer in Expand()\n"; 
        // std::cerr<<"Expand finishes updating graph and buffer\n";
        /** std::cerr<<"After updateGraphandBuffer in Expand()\n";

        std::cerr << "MCTS Search: Expand: Before return" << std::endl; */
        //clear the children of the expanded node
        // std::cerr<<"after expand, expanded node's AOF:\n";
        // for(auto or_node : expanded_node->AOFragment_.or_children_)
        // {
        //     for(auto state : or_node)
        //     {
        //         std::cerr<<state.GetContent()<<"_";
        //     }
        //     std::cerr<<std::endl;
        // }
        // std::cerr<<std::endl;
        // std::cerr<<"posterior after expand:"<< expanded_node->posterior_<<std::endl;

        this->AddChildren(expanded_node);
        // std::cout << "out of addchildren!" << std::endl;
        return expanded_node;
    }
    
	else
	    std::cerr << "\tCannot Expand, TERMINAL NODE\n";
    // If cannot further increase posterior, do not expand

    
    return this->shared_from_this();
}

template<class StateType>
SN_Ptr Online_Learner<StateType>::Checkpoint::Simulate()
{
    /** std::cerr << "MCTS Search: Simulate: inside" << std::endl; */
    std::cerr << "[Simulate()]\n";
    
    std::shared_ptr<Checkpoint> simulate_node = std::make_shared<Checkpoint>(*this);
    //clear simulated node's content
    //set simulation node's parent
    simulate_node->SetParent(0);
    //clear num of neg post gain
    simulate_node-> num_of_neg_post_gain_ = 0;


    // std::cerr<<"data buffer before simulate: \n";
    // for(auto data : this->data_buffer_)
    // {
    //     for(auto state : data)
    //         std::cerr<<state.GetContent()<<"_";
    //     std::cerr<<std::endl;
    // }

    double delete_prob = gen()/(double)(gen.max() - gen.min());
    std::cerr<<"delete prob: "<<delete_prob<<std::endl;
    // int a;
    // std::cin >> a;
    if(delete_prob <= DELETE_PROB)
    {
        std::cerr<<"doing delete!\n";
        simulate_node -> Delete();

    }
    
    simulate_node->SetLevel(simulate_node->GetLevel()+1);
    bool success = simulate_node->Generate(SIMULATE_AOF_ITER);
    

    unsigned simulate_round_count = 0;
    std::cerr << "\tSimulate from LEVEL [" << simulate_node->GetLevel() << "] to LEVEL ";

    while (success && (EARLY_STOP == 0 || simulate_round_count < EARLY_STOP) 
                    && (NEG_GAIN_COUNT == -1 || simulate_node->GetNegGainCount() <= NEG_GAIN_COUNT) )
    {
        ++simulate_round_count;
        std::cerr<<"the "<<simulate_round_count<<"th round in simulation\n";
        
        // std::cout << "AOF after simulate:"<<std::endl;
        // for(auto or_node : simulate_node->AOFragment_.or_children_)
        // {
        //     for(auto state : or_node )
        //     {
        //     std::cerr<<"("<<state.GetContent()<<","<<state.GetId()<<")";                    
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;


        delete_prob = gen()/(double)(gen.max() - gen.min());
        std::cerr<<"delete prob: "<<delete_prob<<std::endl;
        // int b;
        // std::cin >> b;
        if(delete_prob <= DELETE_PROB)
        {
            std::cerr<<"doing delete!\n";
            simulate_node -> Delete();
        }

        simulate_node->SetLevel(simulate_node->GetLevel()+1); 
        success = simulate_node->Generate(SIMULATE_AOF_ITER);
        // std::cout << "AOF after simulate:"<<std::endl;
        // for(auto or_node : simulate_node->AOFragment_.or_children_)
        // {
        //     for(auto state : or_node )
        //     {
        //         auto content = state.GetContent();
        //         if(content == "")
        //             std::cout << state.GetId()<<"_";
        //         else
        //             std::cout << content<<"_";
                
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;
        // std::cout << "Reduciton: ";
        // for(auto it = rd.begin(); it != rd.end();++it)
        // {
        //     for (auto iter = (*it).first.begin(); iter != (*it).first.end(); iter++)
        //         {
        //             std::cout << iter->GetContent() << "_";
        //         }
        //         std::cout << std::endl;

        //         //print positions:
        //     for (auto iter = (*it).second.begin(); iter != (*it).second.end(); iter++)
        //         {
        //             std::cout << "Pos: "
        //                         << "(" << (*iter).first << "), (" << (*iter).second << ")" << std::endl;
        //         }
        // }

    }
    std::cerr << "[" << simulate_node->GetLevel() << "]--" << "[" << simulate_round_count << "]" << " Rounds of Simulation\n";
    
    
    //add prior to the posterior
    // double prior = -ALPHA * simulate_node->graph_.GrammarSize(LEAF_BIAS);
    // simulate_node->posterior_ += prior;
    return simulate_node;
}


template<class StateType>
void Online_Learner<StateType>::Checkpoint:: UpdateWeightsWithEM(const std::vector<SequenceType<StateType> > &raw_memory)
{
    // std::cerr<<"inside EM!\n";
   
    //create a parser
    // std::cerr<<"before update EM\n";
    // this->PrintThirdLevel();
    //  std::cerr<<"the data buffer :\n";            
    // for(auto data : this->data_buffer_)
    // {
    //     for(auto state : data)
    //     {
    //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //     }
    //     std::cerr<<std::endl;
    // }
    // std::cerr<<"the raw data:\n";
    // for(auto data : raw_memory)
    // {
    //     for(auto state : data)
    //     {
    //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //     }
    //     std::cerr<<std::endl;
    // }
    if(NEG_GAIN_COUNT == -1)
    {
        std::vector<Symbolic_Rule<StateType> > all_rules = this->graph_.GetRules();
        std::vector<Symbolic_Rule<StateType> > rules_without_root;
        Symbolic_State<StateType> root_state = this->graph_.GetStateByVertexId(this->graph_.GetRoot());
        for(const  auto &rule : all_rules)
        {
            if(rule.GetSource() != root_state)
            {
                rules_without_root.push_back(rule);
            }
        }

        //store sufficient statistics for all or_rules met 
        std::unordered_map<Symbolic_Rule<StateType>, double> sufficient_statistics;
        std::unordered_map<Symbolic_State<StateType> ,std::unordered_set<Symbolic_Rule<StateType> > >  rules_under_or_node;
        std::unordered_map<VertexId, Symbolic_Rule<StateType> > dummy_to_rule;

   
        // this->graph_.OutputGraph("just_inside_em",PATH,true,SIMPLIFY);

        //parse each data sequence
        for(int i = 0; i < raw_memory.size(); ++i)
        {
            //first try full parse
            std::shared_ptr<grammar<StateType> > g = std::make_shared<grammar<StateType> >(all_rules, std::vector<Symbolic_State<StateType> >{root_state});                  
            this->parser_= std::make_shared<EarleyParser<StateType>>(*g);
            this->parseposition.Clear();

            double all_parse_tree_probs = 1;
            bool parsing_success = false;
            

            int pos = this->parser_->parse(raw_memory[i].begin(), raw_memory[i].end(), std::cout, parsing_success);

            //if cannot fully parse, try partial parse
            if(pos != raw_memory[i].size())
            {
                //recreate parser
                g = std::make_shared<grammar<StateType> >(rules_without_root, this->graph_.GetTopLevelStates());
                this->parser_ = std::make_shared<EarleyParser<StateType> >(*g);
                this->parseposition.Clear();

                this->PartialParseData(raw_memory[i],rules_without_root);
            }

            std::vector<std::vector<typename StateList<StateType>::statelist > > parsed_results = this->parser_->GetPartialParse();
            
            // //error check
            // if(parsed_results.size() != 1)
            // {            
            //     std::cerr<<"Parser one chart_ assumption fails\n the parsed_results size: "
            //             << parsed_results.size()<<std::endl;
            //     throw std::exception();
            // }
            std::vector<std::vector<double> > parse_tree_prob_chart;
            for(int idx = 0; idx < parsed_results.size(); ++idx)
            {
                parse_tree_prob_chart.push_back(std::vector<double>());
                int last_statelist_index = parsed_results[idx].size() - 1;
                double sub_parse_tree_prob = 0;
                std::unordered_map<Symbolic_Rule<StateType>, double> temp_rule_to_prob;
                for(state<StateType> st : parsed_results[idx][last_statelist_index])
                {
                    // find top most root
                    if (st.i_ == 0 && st.j_ != 0 && st.rule_->left().GetContent() == "$" && st.completed())
                    {
                        /** std::cerr << "The state is: [" << st.i_ << "," << st.j_ << "], FROM " << st.rule_->left().GetContent() << " TO ";
                        auto results =  st.rule_->right()[0];
                        for (auto result : results)
                        {
                            std::cerr << result.GetContent() << " ";
                        }
                        std::cerr << std::endl; */
                        // @param prob: probability of this particular possible parsing
                        double parse_tree_prob = 1;
                        std::vector<Symbolic_Rule<StateType> > or_rules_in_one_parse;
                        // @param q: queue used to backtrack all subparsing in one possible parsing
                        std::queue<state<StateType> > q;
                        // @param possible_parsing: possible parsing in one chart
                        Symbolic_State<StateType> possible_parsing_start;
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
                                /** std::cerr << "Back Pointer points to: [" << pointed_state.i_ << "," << pointed_state.j_ << "], FROM " << pointed_state.rule_->left().GetContent() << " TO ";
                                auto results =  pointed_state.rule_->right()[0];
                                for (auto result : results)
                                {
                                    std::cerr << result.GetContent() << " ";
                                }
                                std::cerr << std::endl; */
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
                            VertexId parentId = this->graph_.GetVertexIdByState(st.rule_->left());
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
                            //     // AOFStruct<StateType> peeker;
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
                            if (!this->graph_.GetVertexContent(parentId)->IsAnd())
                            {
                                /** std::cerr << "The rule's source is an Or-node" << std::endl; */
                                // find all out edge weights of the source
                                std::unordered_map<VertexId, double> outEdgeWeights =
                                        this->graph_.GetOutEdgeWeights(parentId, true);

                                // locate the outedge we are looking for
                                // find all right-hand-side vertexids
                                // @param right_hand_state_ids: the ids of right hand side states (e.g. S->NP VP. ids of NP, VP)
                                std::vector<VertexId> right_hand_state_ids;
                                for (Symbolic_State<StateType> right_hand_state : st.rule_->right()[st.right_])
                                    right_hand_state_ids.push_back(
                                            this->graph_.GetVertexIdByState(right_hand_state));
                                // find the children that has all the right-hand-side vertex,
                                // i.e. corresponds to the right-hand-rule
                                // @param dummy_vertices: the ids of the dummy vertices under the Or-node of parent (e.g. S->NP VP. dummy nodes under S)
                                std::vector<VertexId> dummy_vertices =
                                        this->graph_.ChildrenVertices(parentId);

                                for (VertexId id : dummy_vertices)
                                {
                                    // @param found_destination: flag indicating the (e.g. S->NP VP) branch is found
                                    // @param children_vertices: children vertices of the dummy node we are looking at
                                    bool found_destination = true;
                                    std::vector<VertexId> children_vertices =
                                            this->graph_.ChildrenVertices(id);

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
                                        Symbolic_Rule<StateType> cur_rule = {st.rule_->left(), st.rule_->right()[st.right_]};                                
                                        or_rules_in_one_parse.push_back(cur_rule);                        
                                        rules_under_or_node[st.rule_->left()].insert(cur_rule);

                                        //record the mapping from this dummy node to the rule
                                        dummy_to_rule.insert({id, cur_rule});
                                        double weight = outEdgeWeights[id];
                                        parse_tree_prob *= weight;
                                        // std::cout << "weight is: " << weight << ", probability is: " << prob << std::endl;
                                        break;
                                    }
                                    
                                }
                            
                            } 
                        }

                        for(Symbolic_Rule<StateType>& rule : or_rules_in_one_parse)
                        {
                            if(temp_rule_to_prob.find(rule) == temp_rule_to_prob.end())
                                temp_rule_to_prob[rule] = parse_tree_prob;
                            else
                                temp_rule_to_prob[rule] += parse_tree_prob;
                        }

                        parse_tree_prob_chart[idx].push_back(parse_tree_prob);
                        sub_parse_tree_prob += parse_tree_prob;
                    }

                }
                //error checking
                if(sub_parse_tree_prob == 0)
                {
                    std::cerr<<"the data is not partial parsable!\n";
                    throw std::exception();
                }
                for(auto iter = temp_rule_to_prob.begin(); iter != temp_rule_to_prob.end(); ++iter)
                {
                    
                    if(sufficient_statistics.find(iter->first) == sufficient_statistics.end())
                        sufficient_statistics[iter -> first] = iter -> second / sub_parse_tree_prob;
                    else
                        sufficient_statistics[iter -> first] += iter -> second / sub_parse_tree_prob;
                }

            }
            


        }

        std:: unordered_set<Symbolic_State<StateType> > src_states;   
        for(auto iter = sufficient_statistics.begin(); iter != sufficient_statistics.end(); ++iter)
        {
            Symbolic_State<StateType> src = iter->first.GetSource();
            
            //calculate new weights all rules with this source states
            if(src_states.find(src) == src_states.end())
            {
                // double denom = 0;
                src_states.insert(src);
                assert(rules_under_or_node.find(src) != rules_under_or_node.end());
                
                // //calculate new theta's denominator
                // for(Symbolic_Rule<StateType> rule : rules_under_or_node[src])
                // {
                //     assert(sufficient_statistics.find(rule) != sufficient_statistics.end());
                //     denom += sufficient_statistics[rule];
                // }

                VertexId src_id = this->graph_.GetVertexIdByState(src);
                //calcuate new weights
                auto cur_weights = this->graph_.GetOutEdgeWeights(src_id,false);
                for(auto iter = cur_weights.begin(); iter != cur_weights.end(); ++iter)
                {
                    //if this rule under or-node is not used to parse any data
                    if(dummy_to_rule.find(iter->first) == dummy_to_rule.end())
                        cur_weights[iter -> first] = OR_WEIGHT_LRATE * cur_weights[iter -> first];
                    else
                        cur_weights[iter -> first] =  OR_WEIGHT_LRATE * cur_weights[iter -> first] 
                                                    +  (1 - OR_WEIGHT_LRATE)* sufficient_statistics[dummy_to_rule.at(iter->first)];
                }

                //update new weights
                // std::cerr<<"all the weights to be updated in EM: \n";
                // for(auto iter = cur_weights.begin(); iter != cur_weights.end(); ++iter)
                // {
                //     std::cerr<<iter->second<<" ";
                //     if(iter->second == 0)
                //     {
                //         std::cerr<<"update weight in EM error!\n";
                //         throw std::exception();
                //     }
                // }
                // std::cerr<<std::endl;
                this->graph_.SetOutEdgeWeights(src_id, cur_weights);

            }
        }
       
    }
    else
    {
        std::vector<Symbolic_Rule<StateType> > all_rules = this->best_aog_.GetRules();
        std::vector<Symbolic_Rule<StateType> > rules_without_root;
        Symbolic_State<StateType> root_state = this->best_aog_.GetStateByVertexId(this->best_aog_.GetRoot());
        for(const  auto &rule : all_rules)
        {
            if(rule.GetSource() != root_state)
            {
                rules_without_root.push_back(rule);
            }
        }

        //store sufficient statistics for all or_rules met
        std::unordered_map<Symbolic_Rule<StateType>, double> sufficient_statistics;
        std::unordered_map<Symbolic_State<StateType> ,std::unordered_set<Symbolic_Rule<StateType> > >  rules_under_or_node;
        std::unordered_map<VertexId, Symbolic_Rule<StateType> > dummy_to_rule;
    
   
        // this->best_aog_.OutputGraph("just_inside_em",PATH,true,SIMPLIFY);

        //parse each data sequence
        for(int i = 0; i < raw_memory.size(); ++i)
        {
            //first try full parse
            std::shared_ptr<grammar<StateType> > g = std::make_shared<grammar<StateType> >(all_rules, std::vector<Symbolic_State<StateType> >{root_state});                  
            this->parser_= std::make_shared<EarleyParser<StateType>>(*g);
            this->parseposition.Clear();

            double all_parse_tree_probs = 1;
            bool parsing_success = false;
            

            int pos = this->parser_->parse(raw_memory[i].begin(), raw_memory[i].end(), std::cout, parsing_success);

            //if cannot fully parse, try partial parse
            if(pos != raw_memory[i].size())
            {
                //recreate parser
                g = std::make_shared<grammar<StateType> >(rules_without_root, this->best_aog_.GetTopLevelStates());
                this->parser_ = std::make_shared<EarleyParser<StateType> >(*g);
                this->parseposition.Clear();

                this->PartialParseData(raw_memory[i],rules_without_root);
            }

            std::vector<std::vector<typename StateList<StateType>::statelist > > parsed_results = this->parser_->GetPartialParse();
            
            // //error check
            // if(parsed_results.size() != 1)
            // {            
            //     std::cerr<<"Parser one chart_ assumption fails\n the parsed_results size: "
            //             << parsed_results.size()<<std::endl;
            //     throw std::exception();
            // }
            std::vector<std::vector<double> > parse_tree_prob_chart;
            for(int idx = 0; idx < parsed_results.size(); ++idx)
            {
                parse_tree_prob_chart.push_back(std::vector<double>());
                int last_statelist_index = parsed_results[idx].size() - 1;
                double sub_parse_tree_prob = 0;
                std::unordered_map<Symbolic_Rule<StateType>, double> temp_rule_to_prob;
                for(state<StateType> st : parsed_results[idx][last_statelist_index])
                {
                    // find top most root
                    if (st.i_ == 0 && st.j_ != 0 && st.rule_->left().GetContent() == "$" && st.completed())
                    {
                        /** std::cerr << "The state is: [" << st.i_ << "," << st.j_ << "], FROM " << st.rule_->left().GetContent() << " TO ";
                        auto results =  st.rule_->right()[0];
                        for (auto result : results)
                        {
                            std::cerr << result.GetContent() << " ";
                        }
                        std::cerr << std::endl; */
                        // @param prob: probability of this particular possible parsing
                        double parse_tree_prob = 1;
                        std::vector<Symbolic_Rule<StateType> > or_rules_in_one_parse;
                        // @param q: queue used to backtrack all subparsing in one possible parsing
                        std::queue<state<StateType> > q;
                        // @param possible_parsing: possible parsing in one chart
                        Symbolic_State<StateType> possible_parsing_start;
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
                                /** std::cerr << "Back Pointer points to: [" << pointed_state.i_ << "," << pointed_state.j_ << "], FROM " << pointed_state.rule_->left().GetContent() << " TO ";
                                auto results =  pointed_state.rule_->right()[0];
                                for (auto result : results)
                                {
                                    std::cerr << result.GetContent() << " ";
                                }
                                std::cerr << std::endl; */
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
                            VertexId parentId = this->best_aog_.GetVertexIdByState(st.rule_->left());
                            // try
                            // {
                            //     // auto all_states = this->best_aog_.GetStates();
                            //     /** std::cerr<<"size of states in this->best_aog_: "<<all_states.size()<<std::endl; */
                            //     parentId = this->best_aog_.GetVertexIdByState(st.rule_->left());
                                
                            // }
                            // catch(std::exception e)
                            // {
                                
                            //     // for(auto state : all_states)
                            //     //     std::cerr<<state.GetId()<<std::endl;
                            //     std::cerr<<"exception caught in EM!\n";
                            //     throw std::exception();
                            //     // AOFStruct<StateType> peeker;
                            //     /** std::cerr << "AOF id_: " << peeker.id_; */
                            //     // auto all_states = this->best_aog_.GetStates();
                            //     /** std::cerr<<"size of states in this->best_aog_: "<<all_states.size()<<std::endl;
                            //     for(auto state : all_states)
                            //     {
                            //         std::cerr<<"State Content: "<<state.GetContent()<<"_\n"<<"State ID: "<<state.GetId()<<std::endl;
                            //     }

                            //     std::cout << e.what() << std::endl; */
                            // }
                            // if it is Or-node, find the corresponding weight and update likelihood of this possible parsing
                            if (!this->best_aog_.GetVertexContent(parentId)->IsAnd())
                            {
                                /** std::cerr << "The rule's source is an Or-node" << std::endl; */
                                // find all out edge weights of the source
                                std::unordered_map<VertexId, double> outEdgeWeights =
                                        this->best_aog_.GetOutEdgeWeights(parentId, true);

                                // locate the outedge we are looking for
                                // find all right-hand-side vertexids
                                // @param right_hand_state_ids: the ids of right hand side states (e.g. S->NP VP. ids of NP, VP)
                                std::vector<VertexId> right_hand_state_ids;
                                for (Symbolic_State<StateType> right_hand_state : st.rule_->right()[st.right_])
                                    right_hand_state_ids.push_back(
                                            this->best_aog_.GetVertexIdByState(right_hand_state));
                                // find the children that has all the right-hand-side vertex,
                                // i.e. corresponds to the right-hand-rule
                                // @param dummy_vertices: the ids of the dummy vertices under the Or-node of parent (e.g. S->NP VP. dummy nodes under S)
                                std::vector<VertexId> dummy_vertices =
                                        this->best_aog_.ChildrenVertices(parentId);

                                for (VertexId id : dummy_vertices)
                                {
                                    // @param found_destination: flag indicating the (e.g. S->NP VP) branch is found
                                    // @param children_vertices: children vertices of the dummy node we are looking at
                                    bool found_destination = true;
                                    std::vector<VertexId> children_vertices =
                                            this->best_aog_.ChildrenVertices(id);

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
                                        Symbolic_Rule<StateType> cur_rule = {st.rule_->left(), st.rule_->right()[st.right_]};                                
                                        or_rules_in_one_parse.push_back(cur_rule);                        
                                        rules_under_or_node[st.rule_->left()].insert(cur_rule);

                                        //record the mapping from this dummy node to the rule
                                        dummy_to_rule.insert({id, cur_rule});
                                        double weight = outEdgeWeights[id];
                                        parse_tree_prob *= weight;
                                        // std::cout << "weight is: " << weight << ", probability is: " << prob << std::endl;
                                        break;
                                    }
                                    
                                }
                            
                            } 
                        }

                        for(Symbolic_Rule<StateType>& rule : or_rules_in_one_parse)
                        {
                            if(temp_rule_to_prob.find(rule) == temp_rule_to_prob.end())
                                temp_rule_to_prob[rule] = parse_tree_prob;
                            else
                                temp_rule_to_prob[rule] += parse_tree_prob;
                        }

                        parse_tree_prob_chart[idx].push_back(parse_tree_prob);
                        sub_parse_tree_prob += parse_tree_prob;
                    }

                }
                //error checking
                if(sub_parse_tree_prob == 0)
                {
                    std::cerr<<"the data is not partial parsable!\n";
                    throw std::exception();
                }
                for(auto iter = temp_rule_to_prob.begin(); iter != temp_rule_to_prob.end(); ++iter)
                {
                    
                    if(sufficient_statistics.find(iter->first) == sufficient_statistics.end())
                        sufficient_statistics[iter -> first] = iter -> second / sub_parse_tree_prob;
                    else
                        sufficient_statistics[iter -> first] += iter -> second / sub_parse_tree_prob;
                }

            }

            // //calculate probabilities of all parse trees of data x given learned grammar
            // for(int i = 0; i <  parse_tree_prob_chart.size(); ++i)
            // {
            //     double temp = 0;
            //     std::accumulate(parse_tree_prob_chart[i].begin(), parse_tree_prob_chart[i].end(),temp);
            //     all_parse_tree_probs *= temp;
            // }
            // if(all_parse_tree_probs >= 1)
            // {
            //     std::cerr<<"all parse tree prob for this data >= 1!\n";
            //     throw std::exception();
            // }
            
            // //calculate statistics for data x and update sufficient statistics
            // for(auto iter = x_ss.begin(); iter != x_ss.end(); ++iter)
            // {
            //     if(sufficient_statistics.find(iter->first) == sufficient_statistics.end())
            //         sufficient_statistics[iter->first] = iter->second / all_parse_tree_probs;
            //     else
            //         sufficient_statistics[iter->first] += iter -> second / all_parse_tree_probs;
            // }    
            


        }

        std:: unordered_set<Symbolic_State<StateType> > src_states;   
        for(auto iter = sufficient_statistics.begin(); iter != sufficient_statistics.end(); ++iter)
        {
            Symbolic_State<StateType> src = iter->first.GetSource();
            
            //calculate new weights all rules with this source states
            if(src_states.find(src) == src_states.end())
            {
                // double denom = 0;
                src_states.insert(src);
                assert(rules_under_or_node.find(src) != rules_under_or_node.end());
                
                // //calculate new theta's denominator
                // for(Symbolic_Rule<StateType> rule : rules_under_or_node[src])
                // {
                //     assert(sufficient_statistics.find(rule) != sufficient_statistics.end());
                //     denom += sufficient_statistics[rule];
                // }

                VertexId src_id = this->best_aog_.GetVertexIdByState(src);
                //calcuate new weights
                auto cur_weights = this->best_aog_.GetOutEdgeWeights(src_id,false);
                for(auto iter = cur_weights.begin(); iter != cur_weights.end(); ++iter)
                {
                    //if this rule under or-node is not used to parse any data
                    if(dummy_to_rule.find(iter->first) == dummy_to_rule.end())
                        cur_weights[iter -> first] = OR_WEIGHT_LRATE * cur_weights[iter -> first];
                    else
                        cur_weights[iter -> first] =  OR_WEIGHT_LRATE * cur_weights[iter -> first] 
                                                    +  (1 - OR_WEIGHT_LRATE)* sufficient_statistics[dummy_to_rule.at(iter->first)];
                }

                //update new weights
                // std::cerr<<"all the weights to be updated in EM: \n";
                // for(auto iter = cur_weights.begin(); iter != cur_weights.end(); ++iter)
                // {
                //     std::cerr<<iter->second<<" ";
                //     if(iter->second == 0)
                //     {
                //         std::cerr<<"update weight in EM error!\n";
                //         throw std::exception();
                //     }
                // }
                // std::cerr<<std::endl;
                this->best_aog_.SetOutEdgeWeights(src_id, cur_weights);

            }
        }
    }
    
    // std::cerr<<"after update EM\n";
    // this->PrintThirdLevel();
    // std::unordered_map<std::vector<Symbolic_State<StateType> >, double> third_level = this->GetThirdLevel();
    // for(auto iter = third_level.begin(); iter != third_level.end(); ++iter)
    // {
    //     std::cerr<<"the third level seq: \n";
    //     for(auto state : iter->first)
    //     {
    //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
    //     }
    //     std::cerr<<std::endl;
    //     std::cerr<<"this third level's weight: "<<iter->second<<std::endl;
    //     std::cerr<<"it's children or-nodes' weights: \n";
    //     for(auto state : iter->first)
    //     {
    //         if(!state.GetIsBasic() && !this->graph_.GetVertexContent( this->graph_.GetVertexIdByState(state))->IsAnd())
    //         {
    //             std::cerr<<"the or-node: ("<< state.GetContent()<<", "<<state.GetId()<<")\n";
    //             auto recur_weights = this->graph_.GetOutEdgeWeights(this->graph_.GetVertexIdByState(state),false);
    //             for(auto it = recur_weights.begin(); it != recur_weights.end(); ++it)
    //             {
    //                 std::cerr<<it->second <<" ";
    //             }
    //             std::cerr<<std::endl;
    //         }
    //     }
    // }
    // this->graph_.OutputGraph("problem_weight",PATH,true,SIMPLIFY);
    std::cerr<<"Out of EM!\n";
}

template<class StateType>
void Online_Learner<StateType>::Checkpoint::Backtrack(SN_Ptr termination_node)
{
    std::cerr << "[Backtrack]\n";
	Checkpoint* terminal_node = static_cast<Checkpoint *>(termination_node.get());
    // double grammar_size = terminal_node->graph_.GrammarSize(LEAF_BIAS);
	// double prior = -ALPHA * grammar_size;
    
	double reward = terminal_node->posterior_;
    if(NEG_GAIN_COUNT != -1)
    {
        reward = terminal_node -> best_posterior_;
    }
    reward -= this->root_posterior_;
	Search_Node *iter_ptr = this;
	// until reached the root (root's parent is nullptr)
    std::cerr << "\tReward updated is [" << reward << "] for [";
    int count = 0;    
	while (iter_ptr != nullptr)
	{
        count++;
		iter_ptr->Update(reward);
		iter_ptr = iter_ptr->GetParent().get();
	}
    std::cerr << count << "] levels.\n";
}

template<class StateType>
void Online_Learner<StateType>::Checkpoint::SetParent(SN_Ptr parent)
{
    //first create a copy
    // SN_Ptr grandpa = parent->GetParent();
    // std::vector<SN_Ptr> children = parent->GetChildren();
    // unsigned vist_cnt = parent->GetVisitCount();
    // double accum_val = parent->GetValue();
    // SN_Ptr new_par = std::shared_ptr<Checkpoint>(new Checkpoint(grandpa, children, vist_cnt, accum_val));
    // //pass the copy to setParent() in Search_node
    Search_Node::SetParent(parent);
}

//helper function to calculate n choose k

int Combination(int n, int k)
{
    //calculate numerator
    int numerator = 1;
    for(int i = 1; i <= k ; ++i)
        numerator *= n - (k-i);

    //calculate denominator
    int denom = 1;
    for(int i = 1; i <= k; ++i)
        denom *= i;
    
    return numerator/denom;

}
