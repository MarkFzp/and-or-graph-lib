
template <class StateType>
ConfigsMap<StateType> Online_Learner<StateType>::Checkpoint::GenerateKConfigs(unsigned k_gram)
{
	return GenerateKConfigs(k_gram, this->data_buffer_);
}

template <class StateType>
ConfigsMap<StateType> Online_Learner<StateType>::Checkpoint ::GenerateKConfigs(unsigned k_gram, ConfigBuffer<StateType> data_buffer)
{
	/** std::cerr<<"------------------------Inside GenerateKConfigs()---------------------\n"; */
	//k_gram should be at least 1
	if (k_gram < 1)
	{
		std::cerr << "And-Or Fragment should have at least two branches.\n";
		throw std::exception();
	}

	ConfigsMap<StateType> frags;

	/** std::cerr<<"Each data's size in databuffer:\n"; */
	for (auto data : data_buffer)
	{
		/** std::cerr<<data.size()<<std::endl;
		std::cerr<<"this data: ";
		for(auto state : data)
			std::cerr<<state.GetContent()<<"_";
		std::cerr<<std::endl; */

		//skip this example since it is too small to generate k-gram AOF
		if (k_gram > data.size())
			continue;

		//find all fragments with size k in data
		for (int i = 0; i < data.size() - k_gram + 1; i++)
		{
			std::vector<Symbolic_State<StateType>> frag(data.begin() + i,
														data.begin() + i + k_gram);

			if (frags.find(frag) == frags.end())
				frags[frag] = 1;
			else
				frags[frag]++;
		}
	}
	/** std::cerr<<"------------------------Out of GenerateKConfigs()---------------------\n"; */

	return frags;
}

template <class StateType>
AOFStruct<StateType> Online_Learner<StateType>::Checkpoint ::ChooseBigram()
{
	/** std::cerr<<"---------------------------Inside ChooseBigram()------------------------------------\n"; */

	// will be used to obtain a seed(true random numbers) for the random number engine
	std::random_device rd;
	// standard mersenne_twister_engine seeded with rd()
	// std::mt19937 gen(rd());
	double best_score = 0;
	double score = 0;
	std::vector<SequenceType<StateType>> best_configs;
	//use the hyperparameter k to find the best bigram

	AOFStruct<StateType> AOF;
	SequenceType<StateType> config1, config2, config3, config4;

	//first remove same data in data buffer
	std::unordered_set<SequenceType<StateType>> unique_data_set;
	for (auto data : this->data_buffer_)
		unique_data_set.insert(data);

	//turn set into a vector
	std::vector<SequenceType<StateType>> unique_data_vec;
	for (auto data : unique_data_set)
		unique_data_vec.push_back(data);

	//make a distribution
	ConfigsMap<StateType> two_configs_map = GenerateKConfigs(2, unique_data_vec);
	ConfigBuffer<StateType> configs;
	std::vector<int> counts;
	for (auto config : two_configs_map)
	{
		// for (auto state : config.first)
		// {
		// 	std::cerr << "(" << state.GetContent() << ", " << state.GetId() << ", " << config.second << ")";
		// }
		// std::cerr << std::endl;
		configs.push_back(config.first);
		counts.push_back(config.second);
	}

	// std::vector<Symbolic_State<std::string>> test_or_1 = {Symbolic_State<std::string>("Round", true),Symbolic_State<std::string>( "Square",true)};
	// std::vector<Symbolic_State<std::string>> test_or_2 = { Symbolic_State<std::string>( "Square",true),Symbolic_State<std::string>("Round", true)};
	// std::cerr << "the round square frequency: " << two_configs_map[test_or_1] << std::endl
	// 		  << "the square round frequency: " << two_configs_map[test_or_2] << std::endl;
	// int f;
	// std::cin >> f;

	//if no configs of size 2 can be generated, return an empty AOF
	if (configs.empty())
		return AOF;

	std::discrete_distribution<int> dis(counts.begin(), counts.end());
	std::vector<ConfigBuffer<StateType> > sampled_configs;
	std::vector<double> AOF_scores;
	for (int i = 0; i < this->sample_times_; ++i)
	{
		/** std::cerr<<i<<" round\n"; */
		int fst = dis(gen);
		int snd = dis(gen);
		/** std::cerr<<"num of configs passed in: "<<configs.size()<<std::endl;
		std::cerr<<"Each config passed in: \n"; 
		for(auto config : configs)
		{
			for(auto state : config)
			{
				
				std::cerr<<state.GetContent()<<"_";

			}
			std::cerr<<std::endl;
		}

		std::cerr<<"first config idx: "<<fst<<std::endl;
		std::cerr<<"second config idx: "<<snd<<std::endl; */

		config1 = configs[fst];
		config2 = configs[snd];
		
		/** std::cerr<<"configs generated from distribution:\n";
		std::cerr<<"config1: ";
		for(auto state : config1)
			std::cerr<<state.GetContent()<<"_";
		std::cerr<<std::endl;

		std::cerr<<"config2: ";
		for(auto state : config2)
			std::cerr<<state.GetContent()<<"_";
		std::cerr<<std::endl; */

		//if cannot generate a bigram, sampling again
		if (config1[0] == config2[0] && config1[1] == config2[1])
		{
			//--i;
			continue;
		}
		config3 = {config1[0], config2[1]};
		config4 = {config2[0], config1[1]};
		ConfigBuffer<StateType> AOF_configs = {config1, config2, config3, config4};
		sampled_configs.push_back(AOF_configs);
		//update the best bigram, use frequency as metrics for now
		score = BigramScore(AOF_configs, two_configs_map);
		AOF_scores.push_back(score);
		// if (score > best_score)
		// {
		// 	best_score = score;
		// 	best_configs = AOF_configs;
		// }
		// else if(score == best_score)
		// {
		// // 	std::cerr<<"find equal score!\n";
		// // 	std::cerr<<"Current Best Configs Generated: \n";
		// // 	for(auto config : best_configs)
		// // 	{
		// // 		for(auto state : config)
		// // 		{
		// // 			std::cerr<<state.GetContent()<<"_";
		// // 		}
		// // 		std::cerr<<std::endl;
		// // 	} 

		// // 	std::cerr<<"the equal configs: \n";
		// // 	for(auto config : AOF_configs)
		// // 	{
		// // 		for(auto state : config)
		// // 		{
		// // 			std::cerr<<state.GetContent()<<"_";
		// // 		}
		// // 		std::cerr<<std::endl;
		// // 	} 
		// // 	int a;
		// // 	std:: cin >> a;
		// 	if(gen() % 2)
		// 	{
				
		// 		best_score = score;
		// 		best_configs = AOF_configs;
		// 	}
		// }
	}

	//if cannot generate a bigram return an uninitialized AOF
	if (AOF_scores.empty())
	{
		return AOF;
	}

	//sample from the sampled configs to choose one bigram
	std::discrete_distribution<int> score_dis(AOF_scores.begin(), AOF_scores.end());
	best_configs = sampled_configs[score_dis(gen)];

	//initialize AOF
	AOF.and_node_ = Symbolic_State<StateType>(AOF.id_);
	++(AOF.id_);

	/** std::cerr<<"Best Configs Generated: \n";
	for(auto config : best_configs)
	{
		for(auto state : config)
		{
			std::cerr<<state.GetContent()<<"_";
		}
		std::cerr<<std::endl;
	} */

	SequenceType<StateType> or_node_1 = {best_configs[0][0]};
	SequenceType<StateType> or_node_2 = {best_configs[0][1]};

	if (best_configs[0][0] != best_configs[1][0])
		or_node_1.push_back(best_configs[1][0]);

	if (best_configs[0][1] != best_configs[1][1])
		or_node_2.push_back(best_configs[1][1]);

	//set AOF structure's member variables
	AOF.or_children_ = {or_node_1, or_node_2};
	AOF.weights_.push_back(std::vector<double>(or_node_1.size(), 0));
	AOF.weights_.push_back(std::vector<double>(or_node_2.size(), 0));

	// AOFStruct<std::string> test_bigram;

	// test_bigram.or_children_.push_back(test_or_1);
	// test_bigram.or_children_.push_back(test_or_2);
	// std::cerr << "chosen bigram:\n";
	// for (auto or_node : AOF.or_children_)
	// {
	// 	for (auto state : or_node)
	// 		std::cerr << "(" << state.GetContent() << ", " << state.GetId() << ")";
	// 	std::cerr << std::endl;
	// }
	// std::cerr<<"score: "<<best_score<<std::endl;
	// int g;
	// std:: cin >>g;
	// if(test_bigram == AOF)
	// {
	// 	std::cerr<<"find you!\n";
	// 	int d;
	// 	std::cin >> d;
	// }

	
	// int q;
	// std::cin >>q;
	// std::cerr<<"Weights of Each Or node:\n";
	// for(auto or_node : AOF.weights_)
	// {
	// 	for(auto weight : or_node)
	// 		std::cerr<<weight<<"_";
	// 	std::cerr<<std::endl;
	// }

	// std::cerr << "------------------------------Out of ChooseBigram()------------------------------------------\n";
	return AOF;
}

template <class StateType>
ConfigBuffer<StateType> AOFStruct<StateType>::GenerateAllConfigsWrapper()
{
	ConfigBuffer<StateType> all_configs;
	SequenceType<StateType> accum;
	GenerateAllConfigs(all_configs, 0, accum);
	return all_configs;
}

template <class StateType>
ConfigBuffer<StateType> AOFStruct<StateType>::GenerateAllConfigsAndWeights(const ConfigsMap<StateType> &config_map)
{
	ConfigBuffer<StateType> all_configs;
	SequenceType<StateType> accum;
	std::vector<unsigned> counts;
	GenerateAllConfigs(all_configs, 0, accum);

	//clear weights
	for (int i = 0; i < this->weights_.size(); ++i)
		for (int j = 0; j < this->weights_[i].size(); ++j)
			this->weights_[i][j] = 0;

	//find the number of times each config appears in data buffer and update weights
	for (auto config : all_configs)
	{
		auto result = config_map.find(config);
		if (result != config_map.end())
		{
			counts.push_back(result->second);
			for (int i = 0; i < config.size(); ++i)
			{
				//find the idx of the state
				for (int j = 0; j < this->or_children_[i].size(); ++j)
				{
					if (this->or_children_[i][j] == config[i])
					{
						this->weights_[i][j] += result->second;
						break;
					}
				}
			}
		}
		else
			counts.push_back(0);
	}

	//create an index array
	std::vector<int> index;
	for (int i = 0; i < all_configs.size(); ++i)
		index.push_back(i);

	//sort the index vector
	std::sort(index.begin(), index.end(),
			  [&](const int &a, const int &b) { return counts[a] > counts[b]; });

	//reconstruct result vector
	ConfigBuffer<StateType> sorted_result;
	for (auto idx : index)
	{
		sorted_result.push_back(all_configs[idx]);
	}

	this->sorted_configs_ = sorted_result;
	return sorted_result;
}

//recursive function to find all possible configurations
template <class StateType>
void AOFStruct<StateType>::GenerateAllConfigs(ConfigBuffer<StateType> &result,
											  int i,
											  const SequenceType<StateType> &accum)
{
	if (i == this->or_children_.size())
	{
		result.push_back(accum);
		return;
	}

	auto row = this->or_children_[i];

	for (int j = 0; j < row.size(); ++j)
	{
		auto tmp(accum);
		tmp.push_back(row[j]);
		GenerateAllConfigs(result, i + 1, tmp);
	}
}

template <class StateType>
AOFRules<StateType> AOFStruct<StateType>::GetAllRules()
{
	if (!this->or_children_.size())
	{
		AOFRules<StateType> empty_res;
		return empty_res;
	}

	if (this->or_children_.size() == 1 && this->or_children_[0].size() == 1)
	{
		std::cerr << "meaningless AOF!\n";
		throw std::exception();
	}
	//generte or-nodes
	if (this->all_rules_.empty())
	{

		std::vector<Symbolic_State<StateType>> or_nodes;
		for (int i = 0; i < this->or_children_.size(); ++i)
		{
			// if or_children[i] size is 1, skip or_node
			if (this->or_children_[i].size() > 1)
			{
				or_nodes.push_back(Symbolic_State<StateType>(this->id_));
				++AOFStruct::id_;
			}
			else
				or_nodes.push_back(this->or_children_[i][0]);
		}

		//if the AOF has more than one Or-node, add and -> or rule
		if (this->or_children_.size() != 1)
		{
			Symbolic_Rule<StateType> and_to_or(this->and_node_, or_nodes);
			this->all_rules_.push_back(and_to_or);
		}

		//add or to leaf rule
		for (int i = 0; i < this->or_children_.size(); ++i)
		{
			for (int j = 0; this->or_children_[i].size() > 1 && j < this->or_children_[i].size(); ++j)
			{
				SequenceType<StateType> or_child = {this->or_children_[i][j]};
				this->all_rules_.push_back(Symbolic_Rule<StateType>(or_nodes[i], or_child));
			}
		}
	}
	return this->all_rules_;
}

template <class StateType>
double AOFStruct<StateType>::AOFGrammarSize(double leaf_bias) const
{

	double grammar_size = 0;

	if (this->or_children_.size() == 1)
	{
		for (int j = 0; j < this->or_children_[0].size(); ++j)
		{
			if (this->or_children_[0][j].GetId() == -1)
				grammar_size += 1 + leaf_bias;
			else
				grammar_size += 2;
		}
	}
	else
	{
		++grammar_size;
		for (int i = 0; i < this->or_children_.size(); ++i)
		{
			int num_of_or_children = this->or_children_[i].size();
			if (num_of_or_children > 1)
			{
				++grammar_size;
				for (int j = 0; j < num_of_or_children; ++j)
				{
					if (this->or_children_[i][j].GetId() == -1)
						grammar_size += 1 + leaf_bias;
					else
						grammar_size += 2;
				}
			}
			else
			{
				if (this->or_children_[i][0].GetId() == -1)
					grammar_size += leaf_bias;
				else
					++grammar_size;
			}
		}
	}

	return grammar_size;
}

template <class StateType>
bool AOFStruct<StateType>::operator==(const AOFStruct<StateType> &rhs) const
{
	if (this->or_children_.size() != rhs.or_children_.size())
		return false;

	for (int i = 0; i < this->or_children_.size(); ++i)
	{
		if(this->or_children_[i].size() != rhs.or_children_[i].size())
			return false;
		bool equal = std::is_permutation(this->or_children_[i].begin(),
										 this->or_children_[i].end(),
										 rhs.or_children_[i].begin());
		if (!equal)
			return false;
	}

	return true;
}

template <class StateType>
bool AOFStruct<StateType>::operator!=(const AOFStruct<StateType> &rhs) const
{
	return !(*this == rhs);
}
template <class StateType>
int Online_Learner<StateType>::Checkpoint ::FindOrNodeToAdd(const ConfigBuffer<StateType> &config_base,
															const SequenceType<StateType> &frag)
{

	//if only one msimatch, return the mismatch position, else return -1
	int or_size = config_base.size();
	int mis_cnt = 0;
	int idx = -1;
	int j = 0;
	for (int i = 0; i < frag.size(); ++i)
	{
		//if all the prev parts of frag can be found in config_base, return the last state's index
		if (j == or_size)
			return j;

		//if frag[i] is not found in or_config, it should be the only one that is mismatched
		auto result = std::find(config_base[j].begin(), config_base[j].end(), frag[i]);
		if (result == config_base[j].end())
		{

			++mis_cnt;
			idx = i;
			--j;
			//if there is more than one states in frag that cannot be found in config_base, return -1
			if (mis_cnt > 1)
				return -1;
		}
		++j;
	}
	return idx;
}

template <class StateType>
int Online_Learner<StateType>::Checkpoint ::FindStateToAdd(const ConfigBuffer<StateType> &config_base,
														   const SequenceType<StateType> &frag)
{
	int cnt = 0;
	int idx = -1;
	// std::cerr<<"size of config_base: "<<config_base.size()<<std::endl;
	// std::cerr<<"size of frag:"<<frag.size()<<std::endl;
	for (int i = 0; i < frag.size(); ++i)
	{
		auto result = std::find(config_base[i].begin(), config_base[i].end(), frag[i]);
		//if the current state is not found in config_base, record the state's position
		if (result == config_base[i].end())
			idx = i;
		else
			++cnt;
	}

	if (cnt == frag.size() - 1)
		return idx;
	else
		return -1;
}

template <class StateType>
AOFStruct<StateType> Online_Learner<StateType>::Checkpoint ::DeleteOr(AOFStruct<StateType> AOF,
																	  RdBuffer<StateType> &best_rd,
																	  const ConfigsMap<StateType> &configs_map,
																	  double &best_posterior_gain)
{
	// std::cerr<<"-----------------------inside deleteor -----------------------------------\n";
	// std::cerr<<"the passed in AOF:\n";
	// for(auto or_node : AOF.or_children_)
	// {
	// 	for(auto state : or_node)
	// 	{
	// 		std::cerr<<"("<<state.GetContent()<<","<<state.GetId()<<")";
	// 	}
	// 	std::cerr<<std::endl;
	// }
	// std::cerr<<"passed in best_posterior_gain: "<<best_posterior_gain<<std::endl;
	/** std::cerr<<"Inside DeleteOr()\n"; */
	std::random_device rd;
	// std::mt19937 gen(100);
	int num_of_or_nodes = AOF.or_children_.size();
	double posterior_gain;
	AOFStruct<StateType> best_AOF = AOF;

	//refuse to delete an or_node for a bigram
	if (num_of_or_nodes == 1)
		return AOF;

	for (int i = num_of_or_nodes - 1; i >= 0; --i)
	{
		// ++NUM_OF_VARIATION_EACH_ROUND;

		//delete each or_node
		SequenceType<StateType> temp = AOF.or_children_[i];
		AOF.or_children_[i] = AOF.or_children_[num_of_or_nodes - 1];
		AOF.or_children_.pop_back();

		//save original weights
		std::vector<double> temp_row_weight = AOF.weights_[i];
		AOF.weights_[i] = AOF.weights_[num_of_or_nodes - 1];
		AOF.weights_[num_of_or_nodes - 1] = temp_row_weight;
		std::vector<std::vector<double>> temp_weights = AOF.weights_;
		AOF.weights_.pop_back();

		//find all configurations and update weights
		ConfigBuffer<StateType> all_configs = AOF.GenerateAllConfigsAndWeights(configs_map);

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
				Cache<StateType> cache = this->AOF_cache_map_[AOF.or_children_];
				reductions = cache.rd_cache_;
				this->CM_ = cache.cm_cache_;
			}
			else
			{
				reductions = GenerateCMAndReductions(all_configs, this->data_buffer_);
				//cache the generated reductions and CM
				Cache<StateType> cache = {this->CM_, reductions};
				this->AOF_cache_map_[AOF.or_children_] = cache;

			}
		}
		else	
			reductions = GenerateCMAndReductions(all_configs, this->data_buffer_);


		posterior_gain = CalculatePosterior(AOF, reductions);

		// if(posterior_gain + this->posterior_ >  0)
		// {
		// 	std::cerr<<"find you in delete or before end!\n";
		// 	std::cerr<<"current posterior: "<<this->posterior_<<std::endl;
		// 	std::cerr<<"the post gain: "<< posterior_gain<<std::endl;
		// 	int a;
		// 	std:: cin >> a;
		// }

		if (posterior_gain > best_posterior_gain)
		{
			best_posterior_gain = posterior_gain;
			best_AOF = AOF;
			best_rd = reductions;
		}
		else if (posterior_gain == best_posterior_gain)
		{
			if (gen() % 2)
			{

				best_posterior_gain = posterior_gain;
				best_AOF = AOF;
				best_rd = reductions;
			}
		}
		//restore  the config_base for the next round
		AOF.or_children_.push_back(temp);
		AOF.weights_ = temp_weights;
	}
	
	// std::cerr<<"-----------------------Outside DeleteOr -----------------------------------\n";

	return best_AOF;
}

template <class StateType>
AOFStruct<StateType> Online_Learner<StateType>::Checkpoint ::AddOr(AOFStruct<StateType> AOF,
																   const ConfigBuffer<StateType> &k_frag,
																   RdBuffer<StateType> &best_rd,
																   const ConfigsMap<StateType> &configs_map,
																   const ConfigsMap<StateType> &two_config_map,
																   double &best_posterior_gain)
{
	// std::cerr<<"-----------------------inside Addor -----------------------------------\n";
	// std::random_device rd;
	// std::mt19937 gen(417);
	// std::cerr<<"the passed in AOF:\n";
	// for(auto or_node : AOF.or_children_)
	// {
	// 	for(auto state : or_node)
	// 	{
	// 		std::cerr<<"("<<state.GetContent()<<","<<state.GetId()<<")";
	// 	}
	// 	std::cerr<<std::endl;
	// }
	// std::cerr<<"passed in best_posterior_gain: "<<best_posterior_gain<<std::endl;

	typedef std::unordered_multimap<Symbolic_State<StateType>, int> TestMap;
	double posterior_gain;
	std::vector<double> added_weight = {0};
	AOFStruct<StateType> best_AOF = AOF;
	std::vector<int> indexes;
	TestMap tested_states_pos;
	/** std::cerr<<"Inside AddOr()\n"; */
	for (auto frag : k_frag)
	{

		// idx = FindOrNodeToAdd(AOF.or_children_, frag);
		/** std::cerr<<"The idx to be added: "<<idx<<std::endl; */
		//if we can add an or node to make this configuration appear in AOF

		assert(frag.size() - 1 == AOF.or_children_.size() );
		//find if the frag's start or end state can be added as an or-node to AOF
		bool add_at_begin = true;
		bool add_at_end = true;
		for(int i = 1; i < frag.size(); ++i)
		{
			if(std::find(AOF.or_children_[i-1].begin(), AOF.or_children_[i-1].end(),frag[i] ) == AOF.or_children_[i-1].end())
			{
				add_at_begin = false;
				break;
			}
		}
		//check if it will result in weight 0
		if(add_at_begin)
		{
			for(int i = 0; i < AOF.or_children_[0].size(); ++i)
			{
				if(two_config_map.find({frag[0],AOF.or_children_[0][i]}) == two_config_map.end()) 
				{
					add_at_begin = false;
					break;
				}
			}
		}

		for(int i = 0; i < AOF.or_children_.size(); ++i)
		{
			if(std::find(AOF.or_children_[i].begin(), AOF.or_children_[i].end(),frag[i] ) == AOF.or_children_[i].end())
			{
				add_at_end = false;
				break;
			}
		}
		
		if(add_at_end)
		{
			for(int i = 0; i < AOF.or_children_.back().size(); ++i)
			{
				if(two_config_map.find({AOF.or_children_.back()[i], frag.back()}) == two_config_map.end()) 
				{
					add_at_end = false;
					break;
				}
			}
		}


		if(add_at_begin)
			indexes.push_back(0);
		if(add_at_end)
			indexes.push_back(frag.size()-1);


		for(auto idx : indexes)
		{
			//check if the state to be added already been tested at the same pos before
			auto range = tested_states_pos.equal_range(frag[idx]);
			bool tested = std::any_of(range.first, range.second,
									  [idx](typename TestMap::value_type &x) { return x.second == idx; });
			if (tested)
				continue;

			// ++NUM_OF_VARIATION_EACH_ROUND;
			tested_states_pos.insert({frag[idx], idx});

			//std::cerr<<"I am in!\n";
			//save original weights
			std::vector<std::vector<double>> temp_weights = AOF.weights_;

			/** std::cerr<<"Before Insert new or node\n"; */
			//generate all configurations
			SequenceType<StateType> new_or_node = {frag[idx]};
			if (idx < AOF.or_children_.size())
			{
				AOF.or_children_.insert(AOF.or_children_.begin() + idx, new_or_node);
				AOF.weights_.insert(AOF.weights_.begin() + idx, added_weight);
			}
			else
			{
				AOF.or_children_.push_back(new_or_node);
				AOF.weights_.push_back(added_weight);
			}

			ConfigBuffer<StateType> all_configs = AOF.GenerateAllConfigsAndWeights(configs_map);

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
					Cache<StateType> cache = this->AOF_cache_map_[AOF.or_children_];
					reductions = cache.rd_cache_;
					this->CM_ = cache.cm_cache_;
				}
				else
				{
					reductions = GenerateCMAndReductions(all_configs, this->data_buffer_);
					//cache the generated reductions and CM
					Cache<StateType> cache = {this->CM_, reductions};
					this->AOF_cache_map_[AOF.or_children_] = cache;

				}
			}
			else
				reductions = GenerateCMAndReductions(all_configs, this->data_buffer_);

			
			/** std::cerr<<"After update Weight\n"; */
			posterior_gain = CalculatePosterior(AOF, reductions);

			/** std::cerr<<"After calculate likelihood\n"; */
			if (posterior_gain > best_posterior_gain)
			{

				best_posterior_gain = posterior_gain;
				best_AOF = AOF;
				best_rd = reductions;
			}
			else if (posterior_gain == best_posterior_gain)
			{
				if (gen() % 2)
				{

					best_posterior_gain = posterior_gain;
					best_AOF = AOF;
					best_rd = reductions;
				}
			}
			//restore AOF
			/** std::cerr<<"num of or node before erase: "<<AOF.or_children_.size()<<std::endl; */
			AOF.or_children_.erase(AOF.or_children_.begin() + idx);
			AOF.weights_ = temp_weights;
		}
		/** std::cerr<<"prepare for testing next frag\n"; */
	}
	// std::cerr<<"-----------------------Outside Addor ----------------------------------\n";
	
	return best_AOF;
}

template <class StateType>
AOFStruct<StateType> Online_Learner<StateType>::Checkpoint ::AddLeaf(AOFStruct<StateType> AOF,
																	 const ConfigBuffer<StateType> &k_frag,
																	 RdBuffer<StateType> &best_rd,
																	 const ConfigsMap<StateType> &configs_map,
																	 double &best_posterior_gain)
{
	// std::cerr<<"-----------------------inside AddLeaf -----------------------------------\n";
	// std::cout <<" AOF passed in: \n";
	// for(auto or_node : AOF.or_children_)
	// {
	// 	for(auto state : or_node)
	// 	{
	// 		auto content = state.GetContent();
	// 			std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
	// 	}
	// 	std::cout << std::endl;
	// }
	// std::cout << std::endl;
	// std::mt19937 gen(100);

	// std::cerr<<"the passed in AOF:\n";
	// for(auto or_node : AOF.or_children_)
	// {
	// 	for(auto state : or_node)
	// 	{
	// 		std::cerr<<"("<<state.GetContent()<<","<<state.GetId()<<")";
	// 	}
	// 	std::cerr<<std::endl;
	// }
	// std::cerr<<"passed in best_posterior_gain: "<<best_posterior_gain<<std::endl;

	double posterior_gain;
	AOFStruct<StateType> best_AOF = AOF;
	bool tested = false;
	//maintain a map that keeps track of already tested states and its position added to
	std::unordered_multimap<Symbolic_State<StateType>, int> tested_states;

	/** std::cerr<<"Inside AddLeaf()\n"; */
	// check each configuration in configs against config_base
	for (auto frag : k_frag)
	{
		tested = false;
		//find the index of the state in the current frag that can be added to the AOF
		int idx = FindStateToAdd(AOF.or_children_, frag);

		//try every possible states that can be added to the AOF
		if (idx != -1)
		{

			// std::cout <<" AOF before operation: \n";
			// for(auto or_node : AOF.or_children_)
			// {
			// 	for(auto state : or_node)
			// 	{
			// 		auto content = state.GetContent();
			//             std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
			// 	}
			// 	std::cout << std::endl;
			// }
			// std::cout << std::endl;

			//find if the state has already been tested as a children of an or-node
			auto range = tested_states.equal_range(frag[idx]);
			for (auto iter = range.first; iter != range.second; ++iter)
			{
				if (iter->second == idx)
				{
					tested = true;
					break;
				}
			}

			if (tested)
				continue;

			// ++NUM_OF_VARIATION_EACH_ROUND;
			//save original weights
			std::vector<std::vector<double>> temp_weights = AOF.weights_;
			AOF.weights_[idx].push_back(0);

			AOF.or_children_[idx].push_back(frag[idx]);

			// std::cout <<" AOF after push back: \n";
			// for(auto or_node : AOF.or_children_)
			// {
			// 	for(auto state : or_node)
			// 	{
			// 		auto content = state.GetContent();
			//             std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
			// 	}
			// 	std::cout << std::endl;
			// }

			// std::cout << std::endl;
			ConfigBuffer<StateType> all_configs = AOF.GenerateAllConfigsAndWeights(configs_map);
			
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
					Cache<StateType> cache = this->AOF_cache_map_[AOF.or_children_];
					reductions = cache.rd_cache_;
					this->CM_ = cache.cm_cache_;
				}
				else
				{
					reductions = GenerateCMAndReductions(all_configs, this->data_buffer_);
					//cache the generated reductions and CM
					Cache<StateType> cache = {this->CM_, reductions};
					this->AOF_cache_map_[AOF.or_children_] = cache;

				}
			}
			else	
				reductions = GenerateCMAndReductions(all_configs, this->data_buffer_);

			
					
			posterior_gain = CalculatePosterior(AOF, reductions);

			//update best AOF found
			if (posterior_gain > best_posterior_gain)
			{
				best_posterior_gain = posterior_gain;
				best_AOF = AOF;
				best_rd = reductions;
				// std::cout << "likelihood_gain > best_likelihood_gain!\n";
				// if(best_AOF.or_children_[0].size() == 1)
				// 		std::cout << "find you inside equal!\n";
			}
			else if (posterior_gain == best_posterior_gain)
			{
				if (gen() % 2)
				{

					best_posterior_gain = posterior_gain;
					best_AOF = AOF;
					best_rd = reductions;
					// std::cout << "likelihood_gain == best_likelihood_gain!\n";
					// if(best_AOF.or_children_[0].size() == 1)
					// 	std::cout << "find you inside equal!\n";
				}
			}
			//restore the state of configs and configs_base for the next round
			AOF.or_children_[idx].pop_back();
			AOF.weights_ = temp_weights;
			// std::cout <<"restored AOF: \n";
			// for(auto or_node : AOF.or_children_)
			// {
			// 	for(auto state : or_node)
			// 	{
			// 		auto content = state.GetContent();
			//             std::cout << "(" << state.GetContent() << "," << state.GetId() << ")";
			// 	}
			// 	std::cout << std::endl;
			// }
			// std::cout << std::endl;
			//store the state as tested
			std::pair<Symbolic_State<StateType>, int> state(frag[idx], idx);
			tested_states.insert(state);
		}
	}

	// std::cerr<<"this round addleaf's best posterior gain: " <<best_posterior_gain<<std::endl;
	// if(best_AOF.or_children_[0].size() == 1)
	// 	std::cout << "find you before return!\n";
	
	// std::cerr<<"-----------------------Outside AddLeaf -----------------------------------\n";

	return best_AOF;
}

template <class StateType>

AOFStruct<StateType> Online_Learner<StateType>::Checkpoint ::DeleteLeaf(AOFStruct<StateType> AOF, RdBuffer<StateType> &best_rd,
																		const ConfigsMap<StateType> &configs_map,
																		double &best_posterior_gain)
{
	// std::cerr<<"-----------------------inside deleteLeaf -----------------------------------\n";
	// std::cerr<<"the passed in AOF:\n";
	// for(auto or_node : AOF.or_children_)
	// {
	// 	for(auto state : or_node)
	// 	{
	// 		std::cerr<<"("<<state.GetContent()<<","<<state.GetId()<<")";
	// 	}
	// 	std::cerr<<std::endl;
	// }
	// std::cerr<<"passed in best_posterior_gain: "<<best_posterior_gain<<std::endl;
	// std::random_device rd;
	// std::mt19937 gen(10);
	int num_of_nodes = AOF.or_children_.size();
	int num_of_children;
	double posterior_gain;
	AOFStruct<StateType> best_AOF = AOF;
	AOFStruct<StateType> test_aof;
	test_aof.or_children_.push_back({Symbolic_State<StateType>("B_t2",true)});
	test_aof.or_children_.push_back({Symbolic_State<StateType>("L_t1",true)});
	/** std::cerr<<"Inside DeleteLeaf()\n"; */

	//a fast way to delete a leaf node and find its likelihood, O(n)
	for (int i = 0; i < num_of_nodes; ++i)
	{
		//if the current rule is the only child of the or-node, leave it to DeleteOrNode()
		num_of_children = AOF.or_children_[i].size();
		if (num_of_children == 1)
			continue;

		for (int j = num_of_children - 1; j >= 0; --j)
		{
			// ++NUM_OF_VARIATION_EACH_ROUND;
			//delete the ith rule in AOF
			Symbolic_State<StateType> temp_child = AOF.or_children_[i][j];
			AOF.or_children_[i][j] = AOF.or_children_[i][num_of_children - 1];
			AOF.or_children_[i].pop_back();

			//store the original weights to be restored for next round
			int temp_weight = AOF.weights_[i][j];
			AOF.weights_[i][j] = AOF.weights_[i][num_of_children - 1];
			AOF.weights_[i][num_of_children - 1] = temp_weight;
			std::vector<std::vector<double>> temp_weights = AOF.weights_;

			//make the weight matrix's shape align with the or_children matrix
			AOF.weights_[i].pop_back();

			//find all configurations and update weights after deleting the leaf nodes
			ConfigBuffer<StateType> all_configs = AOF.GenerateAllConfigsAndWeights(configs_map);

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
					Cache<StateType> cache = this->AOF_cache_map_[AOF.or_children_];
					reductions = cache.rd_cache_;
					this->CM_ = cache.cm_cache_;
				}
				else
				{
					reductions = GenerateCMAndReductions(all_configs, this->data_buffer_);
					//cache the generated reductions and CM
					Cache<StateType> cache = {this->CM_, reductions};
					this->AOF_cache_map_[AOF.or_children_] = cache;

				}
			}
			else
				reductions = GenerateCMAndReductions(all_configs, this->data_buffer_);

			
			
			posterior_gain = CalculatePosterior(AOF, reductions);
			if( test_aof == AOF)
			{
				std::cerr<<"current best post gain: " << best_posterior_gain << std::endl;
				std::cerr <<"the reduction size: " << reductions.size()<<std::endl;
				std::cerr << "the calculated posterior gain: " << posterior_gain << std::endl;
				for(auto reduction : reductions)
				{
					for(auto state : reduction.first)
						std::cerr << "(" << state.GetContent() << ", "<<state.GetId() << ")"<<std::endl;
					std::cerr << std::endl;
				}
			}
			// if(best_posterior_gain + this->posterior_ >  0)
			// {
			// 	std::cerr<<"find you in delete leaf before end!\n";
			// 	std::cerr<<"current posterior: "<<this->posterior_<<std::endl;
			// 	std::cerr<<"the post gain: "<< posterior_gain<<std::endl;
			// 	int a;
			// 	std:: cin >> a;
			// }
			if (posterior_gain > best_posterior_gain)
			{
				best_posterior_gain = posterior_gain;
				best_AOF = AOF;
				best_rd = reductions;
			}
			else if (posterior_gain == best_posterior_gain)
			{
				if (gen() % 2)
				{

					best_posterior_gain = posterior_gain;
					best_AOF = AOF;
					best_rd = reductions;
				}
			}
			//restore the AOF for the next round
			AOF.or_children_[i].push_back(temp_child);
			AOF.weights_ = temp_weights;
		}
	}
	// std::cerr<<"-----------------------Outside DeleteLeaf -----------------------------------\n";
	
	return best_AOF;
}

template <class StateType>
double Online_Learner<StateType>::Checkpoint::BigramScore(const ConfigBuffer<StateType> &configs,
													   const ConfigsMap<StateType> &configs_map)
{

	std::unordered_set<SequenceType<StateType> > unique_config;
	
	double score = 0;
	for (auto config : configs)
	{
		if(unique_config.find(config) == unique_config.end())
		{
			unique_config.insert(config);
			auto result = configs_map.find(config);
			if (result == configs_map.end())
				continue;
			else
				score += result->second;
		}
		
	}
	return score / unique_config.size();
}

template <class StateType>
long Online_Learner<StateType>::Checkpoint::FindSubBuffer(const SequenceType<StateType> &buf,
														  const SequenceType<StateType> &seq)
{
	if (buf.size() < seq.size())
		return -1;

	long start = 0;
	bool isFound = false;
	// abcde abc
	while (start <= buf.size() - seq.size() - 1)
	{
		int numOfConsecutive = 0;
		if (buf[start] == seq[0])
		{
			numOfConsecutive++;
			for (int i = 0; i < seq.size() - 1; i++)
			{
				if (seq[i] != buf[start + i])
				{
					break;
				}
				else
				{
					numOfConsecutive++;
				}
			}
		}
		if (numOfConsecutive == seq.size())
		{
			isFound = true;
			break;
		}
		else
		{
			start++;
		}
	}
	// long start = 0, b = 0, s = 0;

	// while (s < seq.size())
	// {
	// 	if (buf[b] == seq[s]) b++, s++;
	// 	else
	// 	{
	// 		start++;
	// 		while (buf[start] != seq[0]) start++;
	// 		s = 0, b = start;
	// 		if (start > buf.size() - seq.size() - 1)
	// 		{
	// 			start = -1;
	// 			break;
	// 		}
	// 	}
	// }
	if (isFound)
	{
		return start;
	}
	else
	{
		return -1;
	}
}

template <class StateType>
std::unordered_map<std::vector<Symbolic_State<StateType>>, double> Online_Learner<StateType>::Checkpoint ::GetThirdLevel() const
{
	// root is And
	std::unordered_map<std::vector<Symbolic_State<StateType>>, double> result;
	if (this->graph_.GetVertexContent(this->graph_.GetRoot())->IsAnd())
	{
		std::vector<VertexId> children_vertices_id = this->graph_.ChildrenVertices(this->graph_.GetRoot());
		SequenceType<StateType> seq;
		for (VertexId child_id : children_vertices_id)
			seq.push_back(this->graph_.GetStateByVertexId(child_id));
		result.insert(std::make_pair(seq, 1));
		return result;
	}

	// Root node is or
	std::unordered_map<VertexId, double> weights = this->graph_.GetOutEdgeWeights(this->graph_.GetRoot(), true);
	std::vector<VertexId> dummy_children_vertices_id = this->graph_.ChildrenVertices(this->graph_.GetRoot());
	for (VertexId child_id : dummy_children_vertices_id)
	{
		std::vector<VertexId> third_level_children_id = this->graph_.ChildrenVertices(child_id);

		std::vector<Symbolic_State<StateType>> seq;
		for (int i = 0; i < third_level_children_id.size(); i++)
			seq.push_back(this->graph_.GetStateByVertexId(third_level_children_id[i]));

		result.insert(std::make_pair(seq, weights[child_id]));
	}

	return result;
}

template <class StateType>
void Online_Learner<StateType>::Checkpoint::PrintThirdLevel() const
{
	std::cerr << "------------------------ThirdLevel Summary:------------------------\n";
	if(NEG_GAIN_COUNT == -1)
	{
		// if (this->graph_.GetVertexContent(this->graph_.GetRoot())->IsAnd())
		// {
		// 	std::vector<VertexId> children_vertices_id = this->graph_.ChildrenVertices(this->graph_.GetRoot());
		// 	for (VertexId child_id : children_vertices_id)
		// 	{
		// 		Symbolic_State<StateType> state = this->graph_.GetStateByVertexId(child_id);
		// 		if (state.GetId() == -1)
		// 			std::cerr << state.GetContent() << " ";
		// 		else
		// 			std::cerr << state.GetId() << " ";
		// 	}
		// 	std::cout << "weight: " << this->data_buffer_.size() << "\n";
		// 	return;
		// }

		std::vector<VertexId> dummy_children_vertices_id = this->graph_.ChildrenVertices(this->graph_.GetRoot());

		// Root node is or
		// std::unordered_map<SequenceType<StateType>, double> weights = this->GetThirdLevel();
		// for (VertexId child_id : dummy_children_vertices_id)
		// {
		// 	std::vector<VertexId> third_level_children_id = this->graph_.ChildrenVertices(child_id);
		// 	SequenceType<StateType> seq;
		// 	for (int i = 0; i < third_level_children_id.size(); i++)
		// 	{
		// 		Symbolic_State<StateType> state = this->graph_.GetStateByVertexId(third_level_children_id[i]);
		// 		if (state.GetId() == -1)
		// 			std::cerr << state.GetContent() << " ";
		// 		else
		// 			std::cerr << state.GetId() << " ";
		// 		seq.push_back(state);
		// 	}   
		// 	std::cerr << "weight:" << weights[seq] <<"\n";
		// }
		std::cerr << "data size is: " << dummy_children_vertices_id.size() << "\n";
	}
	else
	{
		// if (this->best_aog_.GetVertexContent(this->best_aog_.GetRoot())->IsAnd())
		// {
		// 	std::vector<VertexId> children_vertices_id = this->best_aog_.ChildrenVertices(this->best_aog_.GetRoot());
		// 	for (VertexId child_id : children_vertices_id)
		// 	{
		// 		Symbolic_State<StateType> state = this->best_aog_.GetStateByVertexId(child_id);
		// 		if (state.GetId() == -1)
		// 			std::cerr << state.GetContent() << " ";
		// 		else
		// 			std::cerr << state.GetId() << " ";
		// 	}
		// 	std::cout << "weight: " << this->data_buffer_.size() << "\n";
		// 	return;
		// }

		std::vector<VertexId> dummy_children_vertices_id = this->best_aog_.ChildrenVertices(this->best_aog_.GetRoot());

		// // Root node is or
		// std::unordered_map<SequenceType<StateType>, double> weights = this->GetThirdLevel();
		// for (VertexId child_id : dummy_children_vertices_id)
		// {
		// 	std::vector<VertexId> third_level_children_id = this->best_aog_.ChildrenVertices(child_id);
		// 	SequenceType<StateType> seq;
		// 	for (int i = 0; i < third_level_children_id.size(); i++)
		// 	{
		// 		Symbolic_State<StateType> state = this->best_aog_.GetStateByVertexId(third_level_children_id[i]);
		// 		if (state.GetId() == -1)
		// 			std::cerr << state.GetContent() << " ";
		// 		else
		// 			std::cerr << state.GetId() << " ";
		// 		seq.push_back(state);
		// 	}   
		// 	std::cerr << "weight:" << weights[seq] <<"\n";
		// }
		std::cerr << "data size is: " << dummy_children_vertices_id.size() << "\n";
	}
	
}


template<class StateType>    
AOFStruct<StateType> Online_Learner<StateType>::Checkpoint::FindAOFFromGraph(Symbolic_State<StateType> aof_root_node)
{
	AOFStruct<StateType> aof;
	aof.and_node_ = aof_root_node;
	VertexId aof_root_id = this->graph_.GetVertexIdByState(aof_root_node);

	//if aof_root is an And-node
	if(this->graph_.GetVertexContent(aof_root_id)->IsAnd())
	{
		std::vector<VertexId> root_children = this->graph_.ChildrenVertices(aof_root_id);
		std::vector<Symbolic_State<StateType> > rule_targets;
		
		//create the first and-> or rule
		for(VertexId child : root_children)
		{
			Symbolic_State<StateType> root_child_state = this->graph_.GetStateByVertexId(child);
			rule_targets.push_back(root_child_state);
		}
		aof.all_rules_.push_back(Symbolic_Rule<StateType>(aof_root_node,rule_targets));

		

		for(int i = 0; i < root_children.size();++i)
		{
			VertexId child_id = root_children[i];
			Symbolic_State<StateType> child_state =this->graph_.GetStateByVertexId(child_id);
			
			//if this child is an And-node, treat it as a Or-node with only one child in AOF
			if(this->graph_.GetVertexContent(child_id)->IsAnd())
			{
				std::vector<Symbolic_State<StateType> > or_child = {child_state};
				aof.or_children_.push_back(or_child);
				aof.weights_.push_back(std::vector<double>(1,1.0));
				
			}

			//if the children is an or-node, first check if it is shared, then find the children under dummy and add it to AOF
			else
			{
				//check if this or-node is shared by other nodes
				if(this->graph_.ParentsVertices(child_id).size() > 1)
				{
					aof.or_children_.push_back({child_state});
					aof.weights_.push_back(std::vector<double>(1,1.0));
					continue;
				}
				
				//add children under dummy to aof
				std::vector<VertexId> dummys = this->graph_.ChildrenVertices(child_id);
				std::unordered_map<VertexId,double> weights = this->graph_.GetOutEdgeWeights(child_id, false);
				aof.or_children_.push_back(std::vector<Symbolic_State<StateType> >());
				aof.weights_.push_back(std::vector<double>());
				for(VertexId dummy : dummys)
				{
					aof.weights_[i].push_back(weights[dummy]);
					std::vector<VertexId> dummy_children = this->graph_.ChildrenVertices(dummy);
					if(dummy_children.size() != 1)
					{
						std::cerr<<"dummy children size: "<<dummy_children.size()<<std::endl;
						for(auto child : dummy_children)
						{
							auto state = this->graph_.GetStateByVertexId(child);
							std::cerr<<state.GetContent()<<"_";			
						}
						std::cerr<<std::endl;
						std::cerr<<"aof root is an and node, Dummy has more than one children\n"<<std::endl;
						throw std::exception();
					}
					Symbolic_State<StateType> dummy_child_state = this->graph_.GetStateByVertexId(dummy_children[0]);
					aof.or_children_[i].push_back(dummy_child_state);
					Symbolic_Rule<StateType> or_to_leaf(child_state, {dummy_child_state});
					aof.all_rules_.push_back(or_to_leaf);
				}
			}
		}
	}
	//if aof_root is an Or-node
	else
	{
		aof.weights_.push_back(std::vector<double>());
		aof.or_children_.push_back(SequenceType<StateType>());
		std::vector<VertexId> dummys = this->graph_.ChildrenVertices(aof_root_id);
		std::unordered_map<VertexId,double> weights = this->graph_.GetOutEdgeWeights(aof_root_id, false);
		
		for(VertexId dummy : dummys)
		{
			std::vector<VertexId> dummy_children = this->graph_.ChildrenVertices(dummy);
			if(dummy_children.size() != 1)
			{
				std::cerr<<"dummy children size: "<<dummy_children.size()<<std::endl;
				for(auto child : dummy_children)
				{
					auto state = this->graph_.GetStateByVertexId(child);
					std::cerr<<"("<<state.GetContent()<<","<<state.GetId()<<")";			
				}
				std::cerr<<std::endl;
				std::cerr<<"aof root is an or node, Dummy has more than one children\n"<<std::endl;
				throw std::exception();
			}
			Symbolic_State<StateType> child_state = this->graph_.GetStateByVertexId(dummy_children[0]);
			aof.or_children_[0].push_back(child_state);
			aof.weights_[0].push_back(weights[dummy]);
			Symbolic_Rule<StateType> rule(aof_root_node,{child_state});
			aof.all_rules_.push_back(rule);
		}
	}
	return aof;
}

// template<class StateType>
// RdBuffer<StateType> Online_Learner<StateType>::Checkpoint::FindReductionCore(const ConfigBuffer<StateType>& sorted_configs,
// 																		const ConfigBuffer<StateType>& reduction_space) const
// {
// 	struct CMP //this is a comparator with a boolean map with the same size as the data buffer
// 	{
// 		std::vector<std::vector<bool> > search_map_;

// 		CMP(const std::vector<SequenceType<StateType> >& related_data)
// 		{
// 			this->search_map_.resize(0);
// 			for(auto data: related_data)
// 			{
// 				this->search_map_.push_back(std::vector<bool>(data.size()));
// 			}
// 		}
// 		//udpate boolean map so that no two configs overlap
// 		void Update(unsigned row, unsigned col, unsigned len)
// 		{
// 			for(unsigned i = 0; i < len; i++)
// 			this->search_map_[row][col + i] = true;
// 		}
// 		//check the boolean map so that no two configs overlap
// 		bool Check(unsigned row, unsigned col, unsigned len)
// 		{
// 			for(unsigned i = 0; i < len; i++)
// 				if(this->search_map_[row][col + i])
// 					return false;
// 			return true;
// 		}
// 	};
// 	CMP cmp(reduction_space);

// 	RdBuffer<StateType> reductions;

// 	//classic nested for loops for each config check all the data
// 	for(auto& config: sorted_configs)
// 	{
// 	//   std::cout << "config is:\n";
// 	//   for (int i = 0; i < config.size(); i++)
// 	//   {
// 	// 	  std::cout << "(" << config[i].GetContent() << "," << config[i].GetId() << ")";
// 	//   }
// 	//   std::cout << "\n";
// 		for(unsigned i = 0; i < reduction_space.size(); i++)
// 		{
// 			auto search_result = std::search(reduction_space[i].begin(),
// 											reduction_space[i].end(),
// 											config.begin(), config.end());
// 			while(search_result != reduction_space[i].end())
// 			{
// 			//   std::cout << "Found\n";
// 				unsigned pos = search_result - reduction_space[i].begin();
// 				if(cmp.Check(i, pos, config.size()))
// 				{
// 					cmp.Update(i, pos, config.size());
// 					try
// 					{
// 						reductions.at(config);
// 					}
// 					catch(const std::out_of_range& oor)
// 					{
// 						std::vector<std::pair<unsigned, unsigned> > positions(0);
// 						reductions[config] = positions;
// 					}
// 					reductions.at(config).push_back(std::make_pair(i, pos));
// 				}
// 				search_result = std::search(search_result + config.size(),
// 											reduction_space[i].end(),
// 											config.begin(), config.end());
// 			}
// 		}
// 	}
// 	return reductions;
// }

// template<class StateType>
// RdBuffer<StateType> Online_Learner<StateType>::Checkpoint::FindReduction(const ConfigBuffer<StateType>& sorted_configs)
// {
// 	return FindReductionCore(sorted_configs, this->data_buffer_);
// }

template <class StateType>
RRMap<StateType> Online_Learner<StateType>::Checkpoint::FindRRMap(AOFStruct<StateType> &AOF)
{
	// std::cout << "------ Inside FindRRMap ------\n";
	RRMap<StateType> rrmap;

	// std::cout << "-----Inside FindRRMap-----\n";
	auto all_rules = AOF.GetAllRules();
	// std::cout<<"all rules from GetAllRules() from RRMap of the chosen aof: \n";
	// for(auto rule : all_rules)
	// {
	// 	AOG_LIB::Symbolic_State<StateType> source = rule.GetSource();
    //     std::vector<AOG_LIB::Symbolic_State<StateType> > results = rule.GetResults();
    //     for (AOG_LIB::Symbolic_State<StateType> result : results)
    //     {
    //         std::cout << "\t" << "(" << source.GetContent() << "," << source.GetId() << ")" 
    //             << " -> " << "(" << result.GetContent() << "," << result.GetId() << ")" << "\n";
    //     }
	// }
	Symbolic_State<StateType> aof_head = all_rules[0].GetSource();
	// std::cout << "AOF_HEAD:" << "(" << aof_head.GetContent() << "," << aof_head.GetId() << ")\n";

	// Find third level children of AOG.
	std::vector<VertexId> dummy_children = this->graph_.ChildrenVertices(this->graph_.GetRoot());
	std::unordered_map<VertexId, double> weights = this->graph_.GetOutEdgeWeights(this->graph_.GetRoot(), false);
	// for(auto iter = weights.begin(); iter != weights.end(); ++iter)
	// {
	// 	std::cerr<<iter->second<<std::endl;
	// }
	ConfigBuffer<StateType> all_third_level_sequence;

	// std::cout << "All third level sequence:\n";
	for (VertexId dummy_child : dummy_children)
	{
		std::vector<VertexId> third_level_children = this->graph_.ChildrenVertices(dummy_child);
		SequenceType<StateType> third_level_seq;
		for (VertexId third_level_child : third_level_children)
			third_level_seq.push_back(this->graph_.GetStateByVertexId(third_level_child));
		// push the number of third level child corresponding to its weight
		// std::cout << "Sequence weight: " << weights[dummy_child] << "\n\t";
		// for (int i = 0; i < third_level_seq.size(); i++)
		// 	std::cout << "(" << third_level_seq[i].GetContent() << "," << third_level_seq[i].GetId() << ")";
		// std::cout << "\n";
		// std::cout << "weights:" << weights[dummy_child] << "||";
		// for (int j = 0; j < third_level_seq.size(); j++)
		// {
		// 	std::cout << "(" << third_level_seq[j].GetContent() << "," << third_level_seq[j].GetId() << ")";
		// }
		// std::cout << "\n";

		for (int i = 0; i < weights[dummy_child]; i++)
		{
			all_third_level_sequence.push_back(third_level_seq);
		}
	}

	bool is_delete = false;
	for (SequenceType<StateType> third_level_seq : all_third_level_sequence)
	{
		// find aof in AOG
		if (std::find(third_level_seq.begin(), third_level_seq.end(), aof_head) != third_level_seq.end())
		{
			is_delete = true;
			break;
		}
	}

	std::vector<std::vector<int>> adjusted_pos;
	for (int i = 0; i < all_third_level_sequence.size(); i++)
		adjusted_pos.push_back(std::vector<int>(all_third_level_sequence[i].size(), 0));

	std::vector<bool> modified_pos(all_third_level_sequence.size(), false);

	ConfigBuffer<StateType> old_third_level_sequence = all_third_level_sequence;
	RdBuffer<StateType> reduction;

	rrmap.is_delete_ = is_delete;
	// if this is add
	if (!is_delete)
	{
		// std::cout << "--------is add-------\n";
		assert(!AOF.sorted_configs_.empty());
		ConfigBuffer<StateType> all_config = AOF.sorted_configs_;
		reduction = this->GenerateCMAndReductions(all_config, all_third_level_sequence, false);
		// std::cerr<<"all third level seq:\n";
		// for(auto seq : all_third_level_sequence)
		// {
		// 	for(auto state : seq)
		// 	{
		// 		std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
		// 	}
		// 	std::cerr<<std::endl;
		// }

		// std::cerr<<"all configs:\n";	
		// for(auto config : all_config)
		// {
		// 	for(auto state : config)
		// 	{
		// 		std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
		// 	}
		// 	std::cerr<<std::endl;
		// }
		// std::cerr<<"all reductions: \n";
		// for(auto entry : reduction)
		// {
		// 	for(auto state : entry.first)
		// 	{
		// 		std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
		// 	}
		// 	std::cerr << std::endl;
		// 	for(auto pos : entry.second)
		// 	{
				
		// 			std::cerr<<"("<<pos.first<<", "<<pos.second<<")";
				
				
		// 	}
		// 	std::cerr<<std::endl;
		// }
		for (auto it = reduction.begin(); it != reduction.end(); it++)
		{
			std::vector<std::pair<unsigned, unsigned>> positions_in_data_buffer = it->second;
			SequenceType<StateType> sequence_to_be_replaced = it->first;

			for (std::pair<unsigned, unsigned> position : positions_in_data_buffer)
			{
				// abcdehh f->cd, pos_x=0, pos_y=2
				// hijskhh g->ehh, pos_x=0, pos_y=4
				int pos_x = position.first;
				int pos_y = position.second - adjusted_pos[pos_x][position.second];

				modified_pos[pos_x] = true;

				// substitute_seq = ab
				SequenceType<StateType> substitute_seq(all_third_level_sequence[pos_x].begin(), all_third_level_sequence[pos_x].begin() + pos_y);

				// std::cout << "First update:";
				// for (int i = 0; i < substitute_seq.size(); i++)
				// 	std::cout << "(" << substitute_seq[i].GetContent() << "," << substitute_seq[i].GetId() << ")";
				// std::cout << "\n";
				// substitute_seq = ab f
				substitute_seq.push_back(aof_head);

				// std::cout << "Second update:";
				// for (int i = 0; i < substitute_seq.size(); i++)
				// 	std::cout << "(" << substitute_seq[i].GetContent() << "," << substitute_seq[i].GetId() << ")";
				// std::cout << "\n";

				// substitute_seq = ab f ehh
				substitute_seq.insert(substitute_seq.end(), all_third_level_sequence[pos_x].begin() + pos_y + sequence_to_be_replaced.size(), all_third_level_sequence[pos_x].end());

				// std::cout << "Third update:";
				// for (int i = 0; i < substitute_seq.size(); i++)
				//     std::cout << "(" << substitute_seq[i].GetContent() << "," << substitute_seq[i].GetId() << ")";
				// std::cout << "\n";

				for (int i = position.second; i < adjusted_pos[pos_x].size(); i++)
					adjusted_pos[pos_x][i] += sequence_to_be_replaced.size() - 1;

				// update data buffer
				all_third_level_sequence[pos_x] = substitute_seq;
			}
		}

	} // If this is delete AOF
	else
	{
		// std::cout << "--------is delete-------\n";
		// for (auto or_child : AOF.or_children_)
		// {
		// 	for (int i = 0; i < or_child.size(); i++)
		// 	{
		// 		std::cout << "(" << or_child[i].GetContent() << "," << or_child[i].GetId() << ") ";
		// 	}
		// 	std::cout << "\n";
		// }
		// std::cout << "\n";

		// find all configurations and its corresponding weight
		std::vector<int> indices(AOF.or_children_.size(), 0);
		ConfigBuffer<StateType> all_config;
		std::vector<double> probability;
		while (1)
		{
			SequenceType<StateType> config;
			double weight = 1.0;
			// print current combination
			for (int i = 0; i < AOF.or_children_.size(); i++)
			{
				config.push_back(AOF.or_children_[i][indices[i]]);
				// if the or_children is Or-node not and
				if (AOF.weights_[i].size() != 1)
					weight *= AOF.weights_[i][indices[i]];
				// if it is And-node, just multiply by 1
				//might be problematic
			}
			all_config.push_back(config);
			probability.push_back(weight);

			// if(AOF.and_node_.GetId() == 73981)
			// {
			// std::cerr<<"all configs: \n";
			// for(auto config : all_config)
			// {
			// 	for(auto state : config)
			// 		std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
			// 	std::cerr<<std::endl;
			// }

			// std::cerr<<"all probs: \n";
			// for(auto prob : probability)
			// {
			// 	std::cerr<<prob<<std::endl;
			// }
			// }
			// find the rightmost array that has more
			// elements left after the current element
			// in that array
			int next = AOF.or_children_.size() - 1;
			while (next >= 0 &&
				   (indices[next] + 1 >= AOF.or_children_[next].size()))
				next--;

			// no such array is found so no more
			// combinations left
			if (next < 0)
				break;

			// if found move to next element in that
			// array
			indices[next]++;

			// for all arrays to the right of this
			// array current index again points to
			// first element
			for (int i = next + 1; i < AOF.or_children_.size(); i++)
				indices[i] = 0;
		}

		// sample enginee
		std::discrete_distribution<int> dis(probability.begin(), probability.end());

		// find all position in AOG that contains the AOF to be deleted
		reduction = GenerateCMAndReductions({{aof_head}}, all_third_level_sequence, false);
		// std::cerr<<"all reductions: \n";
		// for(auto entry : reduction)
		// {
		// 	for(auto state : entry.first)
		// 	{
		// 		std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
		// 	}
		// 	std::cerr << std::endl;
		// 	for(auto pos : entry.second)
		// 	{
				
		// 			std::cerr<<"("<<pos.first<<", "<<pos.second<<")";
				
				
		// 	}
		// 	std::cerr<<std::endl;
		// }
		assert(reduction.size() == 1);

		for (auto it = reduction.begin(); it != reduction.end(); it++)
		{
			std::vector<std::pair<unsigned, unsigned>> positions_in_data_buffer = it->second;
			for (std::pair<unsigned, unsigned> position : positions_in_data_buffer)
			{
				// ab f ehh
				int pos_x = position.first;
				int pos_y = position.second + adjusted_pos[pos_x][position.second];

				modified_pos[pos_x] = true;

				int sample_idx = dis(gen);
				// std::cerr<<"sampled idx: "<<sample_idx<<std::endl;
				
			
				// substitute_seq = ab
				SequenceType<StateType> substitute_seq(all_third_level_sequence[pos_x].begin(), all_third_level_sequence[pos_x].begin() + pos_y);
				// substitute_seq = ab xxx
				substitute_seq.insert(substitute_seq.end(), all_config[sample_idx].begin(), all_config[sample_idx].end());
				// substitute_seq = ab xxx ehh
				substitute_seq.insert(substitute_seq.end(), all_third_level_sequence[pos_x].begin() + pos_y + 1, all_third_level_sequence[pos_x].end());

				for (int i = position.second; i < adjusted_pos[pos_x].size(); i++)
					adjusted_pos[pos_x][i] += all_config[sample_idx].size() - 1;

				all_third_level_sequence[pos_x] = substitute_seq;
			}
		}
	}

	// std::cout << "Modified pos by reduction is:\n";
	// create RRMap with old and new third level sequence
	for (int i = 0; i < all_third_level_sequence.size(); i++)
	{
		// if this data is modified, it forms a new rule in rrmap
		if (modified_pos[i])
		{
			// std::cout << i << "\n";
			auto find_match = rrmap.map_.find(old_third_level_sequence[i]);
			if (find_match == rrmap.map_.end())
			{
				ConfigBuffer<StateType> new_third_level_sequence = {all_third_level_sequence[i]};
				rrmap.map_.insert(std::make_pair(old_third_level_sequence[i], new_third_level_sequence));
			}
			else
				rrmap.map_[old_third_level_sequence[i]].push_back(all_third_level_sequence[i]);
		}
	}
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
	// std::cout << "------ Outside FindRRMap ------\n";

	return rrmap;
}
