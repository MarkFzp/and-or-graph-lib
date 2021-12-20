//
// Created by Luyao Yuan on 17/10/25.
//

#ifndef AOG_LIB_ONLINE_LEARNER_H
#define AOG_LIB_ONLINE_LEARNER_H


#include <utility>
#include <math.h> // log, sqrt
#include <stdlib.h> // exception
#include <algorithm> // find
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include<cassert>
#include "Learner.h"
#include "Context_Matrix.h"
#include "../Utils/MCTS.h"
#include "../Utils/Earley_Parser.h"
#include <limits>
#include <cstdio>
#include <string>

unsigned MEMORY_LIMIT;
unsigned SEARCH_RESOURCE;
unsigned DATA_SET_SIZE;
unsigned BEAM_WIDTH;
unsigned START_EXPAND;
unsigned RANDOM_SEED;
unsigned EXPAND_AOF_ITER;
unsigned SIMULATE_AOF_ITER;
double ALPHA_PRIOR;
double ALPHA_LIKELIHOOD;
double EPSILON_MIN;
double EPSILON_START;
double EPSILON_DECAY;
double Cp;
double THRESHOLD;
double LEAF_BIAS;
double BETA;
double GAMMA;
double REDUCED;
double DELETE_PROB;
double OR_WEIGHT_LRATE;
unsigned EARLY_STOP; // 0 indicates early stopping is not applied to simulation
int NEG_GAIN_COUNT;
bool EXPAND_WITH_STOCHASTIC;
bool TRUNCATE;
bool SIMPLIFY;
bool OFFLINE;
bool CACHE;
bool EM;

std::string PATH;
std::mt19937 gen;

// std::vector<int> NUM_OF_BIGRAM;
// std::vector<double> NUM_OF_VARIATION;
// int TOTAL_VARIATIONS = 0;
// int NUM_OF_VARIATION_EACH_ROUND = 0;
// std :: vector<double> BEST_POSTERIOR;
// std::vector<bool> KL_DIVERGENCE;
int SAMPLING = 1000;

template<class T>
using ConfigBuffer = std::vector<SequenceType<T> >;
template<class T>
using AOFRules = std::vector<AOG_LIB::Symbolic_Rule<T> >;

template<class T>
using RdBuffer = std::unordered_map<SequenceType<T>, std::vector<std::pair<unsigned, unsigned> > >;
template<class T>
using ConfigsMap = std::unordered_map<SequenceType<T>, unsigned>;

// template<class T>
// using RRMap = std::unordered_map<SequenceType<T>, std::vector<SequenceType<T> > >;

// void OutputToFile(std::string filename, const AOG_LIB::T_AOG<std::string> &t_aog, AOG_LIB::VertexId sample_root,int sample_time);

/********************************************************************
 * Current assumption about the AOG representation in this project:
 * 		Root is an or_node with at least one dummy nodes as children 
 * 		even if there is only one data associated with root
 *********************************************************************/
namespace AOG_LIB
{
	class Checkpoint;

	template <class StateType>
	struct RRMap
	{
		std::unordered_map<SequenceType<StateType>, std::vector<SequenceType<StateType> > > map_;
		bool is_delete_;
	};

	template <class StateType>
	struct AOFStruct
	{
		ConfigBuffer<StateType> or_children_;
		std::vector<std::vector<double> > weights_;
		AOG_LIB::Symbolic_State<StateType> and_node_;
		AOFRules<StateType> all_rules_;
		ConfigBuffer<StateType> sorted_configs_;
		static int id_;
		
		//member functions used to find all configurations of this AOF
		ConfigBuffer<StateType> GenerateAllConfigsAndWeights(const ConfigsMap<StateType>& config_map);
		ConfigBuffer<StateType> GenerateAllConfigsWrapper();
		AOFRules<StateType> GetAllRules();
		unsigned NumberOfOrNodes() const;
		double AOFGrammarSize(double leaf_bias=LEAF_BIAS) const;
		bool operator== (const AOFStruct<StateType> &rhs) const; 
		bool operator!= (const AOFStruct<StateType> &rhs) const;
		

		private:
			void GenerateAllConfigs(ConfigBuffer<StateType> &result,
									int i,
		                        const SequenceType<StateType> &accum);
	};

	template <class StateType>
	struct Replacement
	{
		RdBuffer<StateType> reductions_;
		RRMap<StateType> rrmap_;
		bool is_delete_;
	};

	template <class StateType>
	struct Cache
	{
		Context_Matrix<StateType> cm_cache_;
		RdBuffer<StateType> rd_cache_;
	};

	struct ParsePositionStruct
	{
		// position of state that cannot be parsed in a sequence
		// first: the size of left input seq, to calculate position in seq, use seq.size - first
		// second: the number of parsed sequence before this unparsed one
		std::vector<std::pair<int, int> > unparsed_position;

		// position of current successful parsing number
		int parsed_number = 0;

		void Clear()
		{
			this->parsed_number = 0;
			this->unparsed_position.clear();
		}
	};

	template <class StateType>
	int AOFStruct<StateType>::id_ = 1;

	template<class StateType>
	class Online_Learner : public Learner<StateType>
	{
		//checkpoint to be used in MCTS nodes
		class Checkpoint : public Search_Node, public std::enable_shared_from_this<Checkpoint>
		{
		  	std::shared_ptr<EarleyParser<StateType>> parser_;
		  	ParsePositionStruct parseposition;
		  	std::vector<SequenceType<StateType>> data_buffer_;
		  	AOG_LIB::T_AOG<StateType> graph_;
			double posterior_;
			double root_posterior_;
			//current and-or fragment
			//AOFBuffer<StateType> AOFragment_;
			AOFStruct<StateType> AOFragment_;

			//context matrix generated from data buffer and current and-or fragment
			Context_Matrix<StateType> CM_;

			std::unordered_map<ConfigBuffer<StateType>, Cache<StateType> > AOF_cache_map_;

			int num_of_neg_post_gain_;
			double  best_posterior_;
			T_AOG<StateType> best_aog_;
			std::vector<SequenceType<StateType> > best_buffer_;

			//this checkpoint's graph_


			//contain same data as the learner's memory, but updated
			// to be consistent with this checkpoint's graph


			// A hyperparameter used to generate the initial bigram for constructing a AOF
			unsigned sample_times_ = 50;

	

			/*
			 * Data Editing functions
			 * irrelevant to MCTS only related to our task.
			 * Yet, can be called in MCTS related functions
			 *
			 *
			 * Please be notified, these functions changes the checkpoint member variables
			 * If used in expand(), a copy of this checkpoint should already be made.
			 * If used in simulate(), it should be fine to call directly.
			 *
			 *
			 * For a checkpoint, graph, data buffer and likelihood should always be consistent.
			 * Similarly, the and-or fragment and context matrix are always consistent.
			 * Please make sure in all public functions, generate functions and update functions
			 * are called iteratively. Hence, in our code, the function calling order should be:
			 *      Generate();
			 *      ();
			 *      if expand():
			 *          make a copy and be a child of that copy
			 *      UpdateGraph();
			 *      UpdateBuffer();
			 *      CalculateLikelihood();
			 * Potentially wrap the GenerateXXX and UpdateXXX into two functions, Generate(), Update()
			 */
			//TODO select new and-or fragment

			
			AOFStruct<StateType> ChooseBigram();

			ConfigsMap<StateType> GenerateKConfigs(unsigned k_gram,
													ConfigBuffer<StateType> data_buffer);

			ConfigsMap<StateType> GenerateKConfigs(unsigned k_gram);

			bool Generate(int iterations, bool use_stochastic = true);

			RdBuffer<StateType> GenerateCMAndReductions(const ConfigBuffer<StateType> &configurations, 
														const ConfigBuffer<StateType>& reduction_space, 
														bool generate_cm = true);

			long FindSubBuffer(const SequenceType<StateType> &, const SequenceType<StateType> &);

      		//find the subset of configurations generated by AOF that can be used to replace nodes in T_AOG
      		// RdBuffer<StateType> FindReduction(const ConfigBuffer<StateType>& configs /*, bool fst_data = false*/);

			// RdBuffer<StateType> FindReductionCore(const ConfigBuffer<StateType>& configs, const ConfigBuffer<StateType>& reduction_space) const;

			RRMap<StateType> FindRRMap(AOFStruct<StateType> &AOF);

			// update data_buffer and graph, merge graph if needed
			void UpdateGraphAndBuffer(const Replacement<StateType> &replacement);

			double CalculatePosterior(const AOFStruct<StateType> &AOF,
									  const RdBuffer<StateType> &reduction_full,bool use_beta=true, bool is_delete=false);
			
			double CalculateNewFactor1(const AOFStruct<StateType> &AOF,
										const RdBuffer<StateType> &reduction_full);
			
			int CountGraphReplacement(const RdBuffer<StateType> &reductions);

			AOFStruct<StateType> FindAOFFromGraph(Symbolic_State<StateType> aof_root_node);


			/*****************************************************
			  helper functions for selecting best AOF to expand
			 ******************************************************/

			

      		//generate all possible configurations from AOF
			// ConfigBuffer<StateType> GenerateAllConfigs_Wrapper(const AOFBuffer<StateType>& AOF);



			//metric used to evaluate the goodness of the initial bigram
			//now use the sum of the 4 frequency counts to evaluate the bigram
			double BigramScore(const ConfigBuffer<StateType>& configs,
	                        const ConfigsMap<StateType> &configs_map);


			AOFStruct<StateType> DeleteLeaf(AOFStruct<StateType> AOF,
											RdBuffer<StateType>& best_rd,
											const ConfigsMap<StateType> &configs_map,
											double &best_likelihood);

			AOFStruct<StateType> AddLeaf(AOFStruct<StateType> AOF,
										const ConfigBuffer<StateType> &k_frag,
										RdBuffer<StateType> &best_rd,
										const ConfigsMap<StateType> &configs_map,
										double &best_likelihood_gain);

			AOFStruct<StateType> DeleteOr(AOFStruct<StateType> AOF,
										RdBuffer<StateType> &best_rd,
										const ConfigsMap<StateType> &configs_map,
	        							double &best_likelihood_gain);

			AOFStruct<StateType> AddOr(AOFStruct<StateType> AOF,
									const ConfigBuffer<StateType> &k_frag,
									RdBuffer<StateType>& best_rd,
									const ConfigsMap<StateType> &configs_map,
									const ConfigsMap<StateType> &two_config_map,
									double &best_likelihood_gain);


			int FindStateToAdd(const ConfigBuffer<StateType> &config_base,
			                   const SequenceType<StateType> &frag);


			int FindOrNodeToAdd(const ConfigBuffer<StateType> &config_base,
			                    const SequenceType<StateType> &frag);

			void Backtrack(SN_Ptr termination_node);

		public:

			// Checkpoint() : Search_Node()
		  Checkpoint(std::shared_ptr<T_AOG<StateType>> graph_ptr, ConfigBuffer<StateType> data) : Search_Node(), graph_(*graph_ptr), data_buffer_(data), posterior_(0),root_posterior_(0),num_of_neg_post_gain_(0)
		  																		
		  {
			  if(NEG_GAIN_COUNT != -1)
				this->best_posterior_ = std::numeric_limits<double>::lowest();
		  }

			// Checkpoint(SN_Ptr parent, std::vector<SN_Ptr> children, unsigned visit_cnt, double accum_val)
			// 		: Search_Node(parent, children, visit_cnt, accum_val)
			// {}


			Checkpoint(const Checkpoint &other) : Search_Node(other.GetParent(),
			                                            other.GetChildren(),
			                                            other.GetVisitCount(),
			                                            other.GetValue(),
														other.GetLevel()), 
														posterior_(other.posterior_),
														root_posterior_(other.root_posterior_),
														AOFragment_(other.AOFragment_),
														CM_(other.CM_),														
														data_buffer_(other.data_buffer_),
														num_of_neg_post_gain_(other.num_of_neg_post_gain_)														
													
			{
				this->graph_ = T_AOG<StateType>(other.graph_);
				if(NEG_GAIN_COUNT != -1)
				{
					this->best_posterior_ = other.best_posterior_;
					this->best_aog_ = T_AOG<StateType>(other.best_aog_);
					if(!OFFLINE)
						this->best_buffer_ = other.best_buffer_;
				}
			}

			void SetPosterior(double posterior) { this->posterior_ = posterior; }
			void SetBestPosterior(double posterior){this->best_posterior_ = posterior;}
			void SetRootPosterior(double root_posterior){this->root_posterior_ = root_posterior;}
			double GetPosterior() const{ return this->posterior_; };
			double GetBestPosterior() const{return this->best_posterior_;}
			double GetRootPosterior() const{return this->root_posterior_;}

			std::shared_ptr<T_AOG<StateType> > GetGraphPtr() const
			{
				std::shared_ptr<T_AOG<StateType> > p(new T_AOG<StateType>(this->graph_));
				return p;
			}
			std::shared_ptr<T_AOG<StateType> > GetBestGraphPtr() 			
			{ 	
				std::shared_ptr<T_AOG<StateType> > p(new T_AOG<StateType>(this->best_aog_));
				return p;
			}

			void SetBestAOG(std::shared_ptr<T_AOG<StateType> > graph_ptr)
			{
				this->best_aog_ = T_AOG<StateType>(*graph_ptr);
			}
			std::shared_ptr<EarleyParser<StateType> > GetParser() const
			{
				return this->parser_;
			}

			void SetParser(std::shared_ptr<EarleyParser<StateType>> parser)
			{
				this->parser_ = parser;
			}

			ParsePositionStruct GetParsePosition() const
			{
				return this->parseposition;
			}

			bool Delete();

			virtual void SetParent(SN_Ptr);

		
			// Checkpoint& operator= (const Checkpoint&);


			//TODO expand, select and simulate method of MCTS,
			//Heuristic goes here
			/*
			 * Make a copy of self to fill in current place
			 * update self contents
			 * link self to the copy just made
			*/
			SN_Ptr Expand();

			//Potentially random
			/* Roll out from the node just expanded
			 * Should be called by the node just created
			 * Nodes simulated before ending don't need to be stored
			 * but the final node with reward is useful and should be returned
			 * When simulation ends, also do backtrack to update related nodes before return
			*/
			SN_Ptr Simulate();

			//Use UCT algorithm to choose from children
			SN_Ptr Select(double,bool);

			std::vector<SequenceType<StateType> > GetDataBuffer() { return this->data_buffer_; }
			std::vector<SequenceType<StateType> > GetBestDataBuffer() { return this->best_buffer_; }

			void UpdateDataBuffer(const std::vector<SequenceType<StateType> > & memory) {this->data_buffer_ = memory;}
			void UpdateBestDataBuffer(const std::vector<SequenceType<StateType> > & memory) {this->best_buffer_ = memory;}

			AOG_LIB::T_AOG<StateType> GetGraph() { return this->graph_; }

			void SetGraph(const AOG_LIB::T_AOG<StateType> &new_graph) {this -> graph_ = new_graph;}
			std::vector<Symbolic_State<StateType> >
				MergeNewParsingWithGraph(const std::vector<Symbolic_State<StateType> > &);

			void PartialParseData(const std::vector<Symbolic_State<StateType> > &, const std::vector<Symbolic_Rule<StateType> > &);

			std::vector<Symbolic_State<StateType> > GetBestParsing(double &best_probability);

			int EarleyParserWrapper(const std::vector<Symbolic_State<StateType> > &, 
									std::vector<Symbolic_Rule<StateType> >, bool & parsing_success);

			std::unordered_map<SequenceType<StateType>, double> GetThirdLevel() const;	
			void PrintThirdLevel() const;

			SequenceType<StateType> ParseDelete(const SequenceType<StateType> &reparsed_seq, double &reduced_likelihood);

			void UpdateWeightsWithEM(const std::vector<SequenceType<StateType> > &raw_memory);
			// void ClearNegPostGain(){this->num_of_neg_post_gain_ = 0;}
			int GetNegGainCount(){return this->num_of_neg_post_gain_;}
		};

		MCTS searcher_;
		//define the maximum number of CM entries a learner can have
		unsigned memory_limit_;
		double likelihood_;
		std::unordered_set<SN_Ptr> particle_pool;
		//same data as in data buffer but in different representation
		//each checkpoint has its own data buffer, because data can be represented
		//in different ways given different rules(AOG)
		std::vector<SequenceType<StateType> > memory_;
		std::vector<SequenceType<StateType> > raw_memory_;

	public:
		explicit Online_Learner(std::shared_ptr<T_AOG<StateType> >, std::string);

		/*
		 * Given one data sequence:
		 *    TODO partial parse the sequence maybe use Earley parser
		 *    update online learner memory
		 *    attach sequence to the root/merge directly
		 *    --> select and-or fragment
		 *    |   generate context matrix
		 *    |   update likelihood
		 *    |_  update data buffer
		 * @param:
		 *      new_seq: a vector of symbolic state sequence should learn from
		 */
		void Learn(std::vector<Symbolic_State<StateType> > &);

		/* This GetThirdLevel() function returns the third level node's appearance number
         * @param:
         *      void
         * @return:
         *      a map from one dummy node sequence to its appearance number 
         */

		/* This GetThirdLevel() function prints the third level nodes
         * @param:
         *      void
         * @return:
         *      void 
         */

		

	};

	//TODO online learner constructor
	template<class StateType>
	Online_Learner<StateType>::Online_Learner(std::shared_ptr<T_AOG<StateType> > graph_ptr, std::string filename)
			:Learner<StateType>(graph_ptr)
	{
		//initialize AOFStruct's id according to passed in graph
		
		if(!graph_ptr)
			AOFStruct<StateType>::id_ = 0;
		else
		{
			SequenceType<StateType> all_states = this->related_graph_ptr_->GetStates();
			int max_state_id = 0;
			for(auto state : all_states)
			{
				if(state.GetId() > max_state_id)
					max_state_id = state.GetId();
			}

			AOFStruct<StateType> :: id_ = ++max_state_id;
		}
		
    	std::ifstream f(filename, std::ifstream::in);
    	std::string s = "";
		std::unordered_map<std::string, double> config_params;
		while (getline(f, s))
		{
			if (s[0] == '#')
				continue;
			std::stringstream ss(s);
			std::string item;
			std::vector<std::string> tokens;
			while (getline(ss, item, '='))
				tokens.push_back(item);

			// empty line, skip
			if (tokens.size() == 0)
				continue;

			if (tokens.size() != 2)
			{
				std::cerr << "Input Hyper Parameters doesn't match\n";
				throw std::exception();
			}
			config_params.insert(std::make_pair(tokens[0],std::stod(tokens[1])));
		}

		bool find_alpha = false;

		SEARCH_RESOURCE = unsigned(config_params["SEARCH_RESOURCE"]);
		MEMORY_LIMIT = unsigned(config_params["MEMORY_LIMIT"]);
		BEAM_WIDTH = unsigned(config_params["BEAM_WIDTH"]);
		EPSILON_MIN = config_params["EPSILON_MIN"];
		EPSILON_START = config_params["EPSILON_START"];
		EPSILON_DECAY = config_params["EPSILON_DECAY"];
		Cp = config_params["Cp"];
		THRESHOLD = config_params["THRESHOLD"];
		LEAF_BIAS = config_params["LEAF_BIAS"];
		BETA = config_params["BETA"];
		

		if(config_params.find("ALPHA") != config_params.end())
		{
			find_alpha = true;
			ALPHA_PRIOR = config_params["ALPHA"];
			ALPHA_LIKELIHOOD = 1;
		}
		
		if( config_params.find("ALPHA_PRIOR") == config_params.end())
		{
			if(!find_alpha)
				ALPHA_PRIOR = 1.0;
		}
		else
			ALPHA_PRIOR = config_params["ALPHA_PRIOR"];

		if(config_params.find("ALPHA_LIKELIHOOD") == config_params.end())
		{
			ALPHA_LIKELIHOOD = 1;
		}
		else
			ALPHA_LIKELIHOOD = config_params["ALPHA_LIKELIHOOD"];
		
		

		if(config_params.find("GAMMA") == config_params.end())
		{
			GAMMA = 1.0;
		}
		else
			GAMMA = config_params["GAMMA"];

		if(config_params.find("REDUCED") == config_params.end())
		  {
		    REDUCED = 1.0;
		  }
		else
		  REDUCED = config_params["REDUCED"];
		
		if(config_params.find("RANDOM_SEED") == config_params.end())
		{
			std::random_device rd;
			RANDOM_SEED = rd();
		}
		else
			RANDOM_SEED = unsigned(config_params["RANDOM_SEED"]);
		
		if(config_params.find("TRUNCATE") == config_params.end())
		{
			TRUNCATE = true;
		}
		else
			TRUNCATE = bool(config_params["TRUNCATE"]);
		
		if(config_params.find("SIMPLIFY") == config_params.end())
		{
			SIMPLIFY = true;
		}
		else
			SIMPLIFY = bool(config_params["SIMPLIFY"]);
		
		if(SIMPLIFY){
			TRUNCATE = true;
		}
		if(config_params.find("OFFLINE") == config_params.end())
			OFFLINE = true;
		else
			OFFLINE = bool(config_params["OFFLINE"]);
		
		if(config_params.find("CACHE") == config_params.end())
			CACHE = false;
		else
			CACHE = bool(config_params["CACHE"]);

		if(config_params.find("EXPAND_AOF_ITER") == config_params.end())
			EXPAND_AOF_ITER = 20;
		else
			EXPAND_AOF_ITER = unsigned(config_params["EXPAND_AOF_ITER"]);
		
		if(config_params.find("SIMULATE_AOF_ITER") == config_params.end())
			SIMULATE_AOF_ITER = 5;
		else
			SIMULATE_AOF_ITER = unsigned(config_params["SIMULATE_AOF_ITER"]);
		
		if(config_params.find("OR_WEIGHT_LRATE") == config_params.end())
			OR_WEIGHT_LRATE = 0.9;
		else
			OR_WEIGHT_LRATE = double(config_params["OR_WEIGHT_LRATE"]);
		
		if(config_params.find("EARLY_STOP") == config_params.end())
			EARLY_STOP = 0;
		else{
			EARLY_STOP = unsigned(config_params["EARLY_STOP"]);}
			// if(EARLY_STOP == 0){
			// 	std::cerr << "EARLY_STOP cannot be assigned 0!!" << std::endl;
			// 	throw std::exception();
			// }
		
		
		if(config_params.find("NEG_GAIN_COUNT") == config_params.end())
			NEG_GAIN_COUNT = -1;
		else
			NEG_GAIN_COUNT= int(config_params["NEG_GAIN_COUNT"]);

			
		if(config_params.find("EM") == config_params.end())
			EM = true;
		else
			EM = bool(config_params["EM"]);


		DELETE_PROB = config_params["DELETE_PROB"];
		EXPAND_WITH_STOCHASTIC = bool(config_params["EXPAND_WITH_STOCHASTIC"]);
		START_EXPAND = std::min(unsigned(config_params["START_EXPAND"]), MEMORY_LIMIT);

		gen.seed(RANDOM_SEED);
		srand(RANDOM_SEED);
		
		std::cerr << "ALPHA_PRIOR:" << ALPHA_PRIOR <<  "ALPHA_LIKELIHOOD:" << ALPHA_LIKELIHOOD << "," << "SEARCH_RESOURCE:" << SEARCH_RESOURCE << "," << "MEMORY_LIMIT:" << MEMORY_LIMIT << ","
		<< "BEAM_WIDTH:" << BEAM_WIDTH << ",EPSILON_MIN:" << EPSILON_MIN << ",EPSILON_START:" << EPSILON_START << ",EPSILON_DECAY:" << EPSILON_DECAY
		<< ",Cp:" << Cp << ",BETA:" << BETA << ",RANDOM_SEED:" << RANDOM_SEED << "DELETE_PROB:" << DELETE_PROB << "EXPAND_WITH_STOCHASTIC:" << EXPAND_WITH_STOCHASTIC << "\n"; 
	}

	//TODO learn function for the learner
	template<class StateType>
	void Online_Learner<StateType>::Learn(SequenceType<StateType> &new_seq)
	{
		static unsigned data_index = 0;
		data_index++;

		std::cout<<"the passed in new seq:\n";
		for(auto state : new_seq)
		{
			std::cout<<state.GetContent()<<"_";
		}
		std::cout<<std::endl;
		
		if(!this->related_graph_ptr_->HasRoot())
		{
			std::cerr<<"taog does not have a root!\n";
			Symbolic_State<StateType> source_state(0);
			std::shared_ptr<AOG_Vertex<StateType> > source_ptr(new AOG_Vertex<StateType>(source_state, true, true));
			this->related_graph_ptr_->AddVertex(source_ptr);
		}

		// create checkpoint
		std::shared_ptr<Checkpoint> checkpoint_root = std::make_shared<Checkpoint>(this->related_graph_ptr_,
		                                                                           this->memory_);

		//Print the graph of this checkpoint:
		// checkpoint_root->GetGraph().OutputGraph("checkpoint", PATH, false);
		
		SequenceType<StateType> updated_seq = checkpoint_root->MergeNewParsingWithGraph(new_seq);
		
		this->related_graph_ptr_ = std::make_shared<AOG_LIB::T_AOG<StateType> >(checkpoint_root->GetGraph());
		
		//reached memory limit 
		if (this->memory_.size() == MEMORY_LIMIT)
		{
			assert(this->raw_memory_.size() == MEMORY_LIMIT);
			this->memory_.erase(this->memory_.begin());
			this->raw_memory_.erase(this->raw_memory_.begin());
		}
		
		this->memory_.push_back(updated_seq);
		this->raw_memory_.push_back(new_seq);
		
		// calculate likelihood of Checkpoint root
		double log_initial_probability = 0;
		std::unordered_map<SequenceType<StateType>, double> weights =
				checkpoint_root->GetThirdLevel();
		for (int i = 0; i < memory_.size(); i++)
			log_initial_probability += log(weights[memory_[i]]);

		log_initial_probability *= ALPHA_LIKELIHOOD;

		checkpoint_root->UpdateDataBuffer(this->memory_);
		checkpoint_root->UpdateBestDataBuffer(this->memory_);
		checkpoint_root->SetBestAOG(checkpoint_root->GetGraphPtr());

		double grammar_size = checkpoint_root->GetGraph().GrammarSize(1);
		std::cout << "grammar size: "<<grammar_size<<std::endl;
		double prior = -ALPHA_PRIOR * grammar_size;
		std::cerr<<"initial prior: "<<prior<<std::endl;
		std::cerr<<"initial log likelihood: " << log_initial_probability<<std::endl;

		checkpoint_root->SetPosterior(prior + log_initial_probability);
		checkpoint_root -> SetBestPosterior(prior + log_initial_probability);
		checkpoint_root -> SetRootPosterior(prior+log_initial_probability);
		this->searcher_.DeleteRoot();
		this->searcher_.UpdateRoot(checkpoint_root);
		this->searcher_.UpdateResource(SEARCH_RESOURCE);
		
		//Use offline learning to test MCTS
		if(this->memory_.size() >= START_EXPAND)
		{
			// record the current posterior
			std::cerr<<"initial posterior: "<< checkpoint_root->GetPosterior()<<std::endl;
			std::cerr<<"number of rules in the initial graph: "<<this->related_graph_ptr_->NumOfRules()<<std::endl;
			
			// BEST_POSTERIOR.push_back(checkpoint_root->GetPosterior());
       		// OutputToFile("../AOG_Sample_Results/sample_"+std::to_string(BEST_POSTERIOR.size()-1),checkpoint_root->GetGraph(),checkpoint_root->GetGraph().GetRoot(),SAMPLING);                    
			// KL_DIVERGENCE.push_back(true);
       		// OutputToFile("../AOG_Sample_Results/sample_"+std::to_string(BEST_POSTERIOR.size()-1),checkpoint_root->GetGraph(),checkpoint_root->GetGraph().GetRoot(),SAMPLING);                    


			std::vector<SN_Ptr> search_results = this->searcher_.Search();

			//find the best searching result
			Checkpoint *current_best_cpt = checkpoint_root.get();
			double current_best_likelihood = checkpoint_root->GetPosterior();

			std::shared_ptr<T_AOG<StateType> > cpt_graph = current_best_cpt->GetGraphPtr();
			std::vector<Symbolic_Rule<StateType>> rules = cpt_graph->GetRules();

			if(OFFLINE){
				std::ofstream file_pos;

				file_pos.open(PATH + "all_posteriors.txt", std::ofstream::out | std::ofstream::app);

				if(file_pos.is_open())
				{
					file_pos << "Root Posterior: " << checkpoint_root->GetPosterior() << "  Root Prior: " << -ALPHA_PRIOR * (cpt_graph->GrammarSize(1))<< std::endl;
				}
				else{
					std::cout << "Unable to open all_posteriors.txt" << std::endl;
				}

				file_pos.close();
			}

			
			std::cerr<<"start posterior is: "<<current_best_likelihood<<std::endl;

			//if no early neg count stop
			if(NEG_GAIN_COUNT == -1)
			{
				for (auto search_result: search_results)
				{
					Checkpoint *cpt_result = static_cast<Checkpoint *>(search_result.get());
					std::cerr<<"posterior of search result without neg gain: "<<cpt_result->GetPosterior()<<std::endl; 
					if (cpt_result->GetPosterior() > current_best_likelihood)
					{						
						/* std::cerr<<"Inside if statement in Learn()\n"; */
						current_best_cpt = cpt_result;
						current_best_likelihood = cpt_result->GetPosterior();
					}
				}
			}
			else
			{
				for (auto search_result: search_results)
				{
					Checkpoint *cpt_result = static_cast<Checkpoint *>(search_result.get());
					std::cerr<<"posterior of search result: "<<cpt_result->GetBestPosterior()<<std::endl; 
					if (cpt_result->GetBestPosterior() > current_best_likelihood)
					{						
						std::cerr<<"Inside if statement in Learn()\n"; 
						std::cerr<<"current cpt result: "<< cpt_result -> GetBestGraphPtr()->GetRoot()<<std::endl;
						current_best_cpt = cpt_result;
						current_best_likelihood = cpt_result->GetBestPosterior();
					}
				}
			}
			
			if(current_best_cpt)
			{
				
				if(!OFFLINE && EM)
				{
					//update the weight using EM if incremental learning
					current_best_cpt -> UpdateWeightsWithEM(this->raw_memory_);
				}
				
				if(NEG_GAIN_COUNT == -1)
				{
					// std::cout << "picked best posterior: " << current_best_cpt->GetPosterior() << std::endl;
					this->related_graph_ptr_ = current_best_cpt->GetGraphPtr();
					// std::cout << "grammar size: " << this->related_graph_ptr_->GrammarSize(1) << std::endl;
					this->memory_ = current_best_cpt->GetDataBuffer();
				}
				else
				{
					this->related_graph_ptr_ = current_best_cpt -> GetBestGraphPtr();
					if(!OFFLINE)
						this->memory_ = current_best_cpt -> GetBestDataBuffer();
					else
						this->memory_ = current_best_cpt -> GetDataBuffer();
				}
			}

			// Print out the learned_tree for all the search results:
			std::cerr << "Start to print learned tree...\n";
			int search_results_index = search_results.size();

			// from last resource to first
			std::ofstream file_pos;
			file_pos.open(PATH + "all_posteriors.txt", std::ofstream::out | std::ofstream::app);
			if(!file_pos.is_open()){
				std::cerr << "Unable to open all_posteriors.txt" << std::endl;
			}
			if(OFFLINE){
				for(auto search_result : search_results)
				{
					Checkpoint *cpt_result = static_cast<Checkpoint *>(search_result.get());
					std::shared_ptr<T_AOG<StateType> > cpt_graph = cpt_result->GetGraphPtr();
					cpt_graph -> OutputLearnedTree(PATH, std::to_string(search_results_index));

					file_pos << "Resource#: " << std::to_string(search_results_index) << " Posterior: " << cpt_result->GetPosterior() << "  Prior: " << -ALPHA_PRIOR * (cpt_graph->GrammarSize(1)) << std::endl;
					// cpt_graph->PrintOrNodeAndChildrenCount(file_pos);

					search_results_index--;
				}
			}
			else{
				std::shared_ptr<T_AOG<StateType> > current_best_graph = current_best_cpt->GetGraphPtr();
				
				current_best_graph -> OutputLearnedTree(PATH, std::to_string(1 + data_index - START_EXPAND));

				file_pos << "Resource#: " << std::to_string(1 + data_index - START_EXPAND) << " Posterior: " << current_best_cpt->GetPosterior() << "  Prior: " << -ALPHA_PRIOR * (current_best_graph->GrammarSize(1)) << std::endl;
				// current_best_graph->PrintOrNodeAndChildrenCount(file_pos);

			}
			file_pos.close();


			current_best_cpt->PrintThirdLevel();
			
		}

		
	}
	/********************************************/
	/********* Partial Parser functions *********/
	/********************************************/

	#include "Partial_Parser_Utils.h"

	/****************************************/
	/********* CheckPoint functions *********/
	/****************************************/
	#include "Checkpoint_Utils.h"
	
	#include "Checkpoint_Functions.h"

}
namespace std
{
	template <class StateType>
    struct hash<std::pair<SequenceType<StateType>, int> >
    {
        size_t operator()(const pair<SequenceType<StateType>, int> &rd) const noexcept
        {
            return (hash<SequenceType<StateType> >{}(rd.first) ^ (hash<int>{}(rd.second)<<1));
        }
    };

	template<class StateType>
	struct hash<AOG_LIB::AOFStruct<StateType> >
	{
		size_t operator()(const AOG_LIB::AOFStruct<StateType> &aof) const noexcept
		{
			size_t hash_value = 0;
			for(const auto& or_node : aof.or_children_)
			{
				size_t temp = 0;
				for(const auto& state : or_node)
					temp += hash<AOG_LIB::Symbolic_State<StateType> >{}(state);
				hash_value = (hash_value << 1) ^ temp;
			}
			return hash_value;
		}
	};
}

//output the sample result to a file
void OutputToFile(std::string filename, const AOG_LIB::T_AOG<std::string> &t_aog, AOG_LIB::VertexId sample_root,int sample_time)
{
	//open file
	std::ofstream outFile;
	outFile.open(filename);
	if (outFile.is_open())
	{
		
		while (sample_time)
		{
			std::vector<AOG_LIB::VertexId> res;
			double prob;
			t_aog.Sample(sample_root, res, prob);
			sample_time--;
			for (AOG_LIB::VertexId id : res)
			{
				outFile << t_aog.GetStateByVertexId(id).GetContent() << " ";
			}
			outFile << "\n";
		}
	}
	else
	{
		std::cerr << "Error opening file!\n";
		throw std::exception();
	}
	outFile.close();
}

#endif //AOG_LIB_ONLINE_LEARNER_H
