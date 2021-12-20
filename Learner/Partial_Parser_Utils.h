#include "../Utils/Earley_Parser.h"


// parse from the position where earley parser stops parsing
// i.e. s1,s2,s3,s4, if only s1 and s2 can be parsed, next time start from s3 as a new input_seq
template<class StateType>
void Online_Learner<StateType>::Checkpoint::PartialParseData(const SequenceType<StateType> &input_seq, const std::vector<Symbolic_Rule<StateType> > &rules_to_parse)
{
    /** std::cerr << "PartialParseData:Inside PartialParseData" << std::endl; */

    if (input_seq.size() == 0)
        return;
    
    /** std::cerr << "PartialParseData:Before EarleyParserWrapper" << std::endl; */
    // parse the sequence, record the position of the first unparsed data
    bool parsing_success = false;
    int pos = this->EarleyParserWrapper(input_seq, rules_to_parse, parsing_success);
    
    //  std::cerr << "PartialParseData:After EarleyParserWrapper" << std::endl;
    // std::cerr << "PartialParseData:position is: " << pos << std::endl;
    // std::cerr<<"parsing_success: "<<parsing_success<<std::endl;

    // std::cerr<<"all rules used to parse: \n";
    // auto all_rules = this->graph_.GetRules();
    // for(auto rule : all_rules)
    // {
    //     std::cerr<<"("<<rule.GetSource().GetContent()<<", "<< rule.GetSource().GetId()<<")->";
    //     for(auto target : rule.GetResults())
    //     {
    //         std::cerr<<"("<<target.GetContent()<<", "<<target.GetId()<<")";
    //     }
    //     std::cerr<<std::endl;
    // }
    // if it cannot be parsed with this start, record the unparsed position
    // bool and_or_node_id_found = false;
    // if(std::find(and_or_node_ids.begin(), and_or_node_ids.end(), input_seq[0].GetId()) != and_or_node_ids.end())
    //     and_or_node_id_found = true;
    // std::cerr << "The value of pos returned: " << pos << "\n";
    if (pos == 1 && !parsing_success)//&& !and_or_node_id_found)
    {
        // std::cerr << "Before unparsed position is pushed back\n";
        this->parseposition.unparsed_position.push_back(std::make_pair(input_seq.size(), this->parseposition.parsed_number));

    }

    // else increment the parsed number
    else
        this->parseposition.parsed_number++;

    // if finish parsing, return
    if (pos == input_seq.size())
    {
        /** std::cerr << "PartialParseData:return" << std::endl; */
        return;
    }
    // else, parse from next unparsed position
    std::vector<Symbolic_State<StateType> > unparsed_input_seq(input_seq.begin() + pos, input_seq.end());
    this->PartialParseData(unparsed_input_seq, rules_to_parse);
}

// for each parsing input, try the largest partial parse by reducing the size of rules
// i.e. S -> VP NP, try S first, if cannot parse S, try from next level (VP/NP)
template<class StateType>
int Online_Learner<StateType>::Checkpoint::EarleyParserWrapper(const SequenceType<StateType> &input_seq,
                                                    std::vector<Symbolic_Rule<StateType> > all_rules, bool & parsing_success)
{

    /** std::cerr << "EarleyParserWrapper:Inside EarleyParserWrapper" << std::endl;
    std::cerr << "EarleyParserWrapper:Before parser.parse" << std::endl;
    // TODO: debugging purpose, later change to the second line
    std::cerr << "EarleyParserWrapper:size of input_seq: " << input_seq.size() << std::endl; */
    // if all rules are eliminated, no rule can parse the input sequence
    // move to the next input position since there is no parsing for starting at this position
    if (all_rules.size() == 0)
    {
        /** std::cerr << "EarleyParserWrapper:NO RULE TO PARSE THIS INPUT" << std::endl; */
        return 1;
    }

    this->parser_->update(all_rules, input_seq);
    int pos = this->parser_->parse(input_seq.begin(), input_seq.end(),std::cout, parsing_success);
    // int pos = this->parser_->parse(input_seq.begin(), input_seq.end());
    //  std::cerr << "EarleyParserWrapper:After parser.parse" << std::endl;
    // std::cerr << "EarleyParserWrapper:Position is:" << pos << std::endl; 
    // std::cerr << "Parse success: "<<parsing_success<<std::endl;
    // earley parser parse with rules & states
    if (pos == 0)
    {
        // std::cerr << "0 returned, parse with reduced rules, else return position \n";
        /** std::cerr << "EarleyParserWrapper:Before deleting rules" << std::endl; */
        // if cannot parse from start, erase all top level rules
        SequenceType<StateType> top_level_rules = this->parser_->get_top_level_rules(all_rules);
        // if cannot further erase top level rules, no rule can parse the input sequence
        if (top_level_rules.empty())
        {
            std::cerr << "NO WAY TO FURTHER TRUNCATE THE RULES" << std::endl; 
            return 1;
        }
        /** std::cerr << "Top Level Rules start from"; 
        for (auto top_state : top_level_rules)
        {
            std::cerr << " " << top_state.GetContent();
        }
        std::cerr << std::endl; */

        for (typename std::vector<Symbolic_Rule<StateType> >::iterator it = all_rules.begin();
                it < all_rules.end(); it++)
        {
            if (std::find(top_level_rules.begin(), top_level_rules.end(), it->GetSource()) != top_level_rules.end())
            {
                /** std::cerr << "EarleyParserWrapper:Erased rule is From " << it->GetSource().GetContent() << " To ";
                for (auto result : it->GetResults())
                {
                    std::cerr << result.GetContent() << " ";
                }
                std::cerr << std::endl; */
                it = all_rules.erase(it);
                it--;               
            }
        }
        //  std::cerr << "EarleyParserWrapper:rules left is: " << all_rules.size() << std::endl; 
        // std::cerr << "EarleyParserWrapper:rules are: " << std::endl; 
        // for (auto rule : all_rules)
        // {
        //     std::cerr << "From " << rule.GetSource().GetContent() << " To ";
        //     auto results =  rule.GetResults();
        //     for (auto result : results)
        //     {
        //         std::cerr << result.GetContent() << " ";
        //     }
        //     std::cerr << std::endl;
        // } 
        return this->EarleyParserWrapper(input_seq, all_rules, parsing_success);
    }

    // else, return the position that has been already parsed
    return pos;
}

template<class StateType>
std::vector<Symbolic_State<StateType> > Online_Learner<StateType>
    ::Checkpoint::MergeNewParsingWithGraph(const SequenceType<StateType> &input_seq)
{
    std::vector<Symbolic_Rule<StateType> > all_rules = this->graph_.GetRules();
     std::vector<Symbolic_Rule<StateType> > rules_without_root;
    Symbolic_State<StateType> root_state = this->graph_.GetStateByVertexId(this->graph_.GetRoot());
    for (const auto &rule : all_rules)
    {
        if (rule.GetSource() != root_state)
        {
           rules_without_root.push_back(rule);
        }
    }
    std::shared_ptr<grammar<StateType> > g = std::make_shared<grammar<StateType> >(rules_without_root,
		                                                        this->graph_.GetTopLevelStates());
	this->parser_ = std::make_shared<EarleyParser<StateType> >(*g);

    // clear position record
    this->parseposition.Clear();

    this->PartialParseData(input_seq, rules_without_root);
    double reduced_likelihood;
    // @param states_under_dummy: the Symbolic_States right under the dummy node (i.e. no other rule points to them)
    std::vector<Symbolic_State<StateType> > states_under_dummy = this->GetBestParsing(reduced_likelihood);
    


   
    // std::cout << "Check: after getting bestParsing, the reduced_likelihood is: " << reduced_likelihood << "\n";
    for (int i = 0; i < this->parseposition.unparsed_position.size(); i++)
    {
        // std::cerr<<"inside unparsed position\n";
        // already parsed + unparsed before itself
        int already_parsed_start = this->parseposition.unparsed_position[i].second;
        int unparsed_position_in_seq = input_seq.size() - this->parseposition.unparsed_position[i].first;
        //std::cout << "Print out the inserted symbolic state. ID: " << input_seq[unparsed_position_in_seq].GetId() <<
        //"  Content: " << input_seq[unparsed_position_in_seq].GetContent() << "\n";
        states_under_dummy.insert(states_under_dummy.begin() + already_parsed_start + i, 
                                                input_seq[unparsed_position_in_seq]);
        
    }

   
//    std::cout << "Best Parsing is:" << std::endl;
//    for (auto rule : best_parsing)
//    {
//        std::cout << "From (" << rule.GetSource().GetContent()<< ","<< rule.GetSource().GetId()<<") To";
//        for (auto state : rule.GetResults())
//    	    std::cout << "(" <<  state.GetContent()<<","<<state.GetId()<<") ";
//        std::cout << std::endl;
//    }
    /** std::cerr << "MergeNewParsingWithGraph:best_parsing size is:" << best_parsing.size() << std::endl; */

    
   //    std::cout << "before:\n";
//    for (auto state: states_under_dummy)
//        std::cout << "(" << state.GetContent() << "," << state.GetId() << ")" << " ";
//    std::cout << "\n";

    // @param candidate_dummy_nodes: dummy nodes right below the root of related_graph_ptr
    std::vector<VertexId> candidate_dummy_nodes = this->graph_.ChildrenVertices(
            this->graph_.GetRoot());

    for (VertexId candidate_dummy_node : candidate_dummy_nodes)
    {
        // @param candidate_dummy_children_id: ids of children of one specific dummy node
        std::vector<VertexId> candidate_dummy_children_id = this->graph_.ChildrenVertices(
                candidate_dummy_node);
        std::vector<Symbolic_State<StateType> > candidate_dummy_children;
        // @param candidate_dummy_children: all children of dummy node
        for (VertexId candidate_dummy_child_id: candidate_dummy_children_id)
            candidate_dummy_children.push_back(
                    this->graph_.GetStateByVertexId(candidate_dummy_child_id));

        bool is_match = true;
        // no match
        if (candidate_dummy_children.size() != states_under_dummy.size())
            continue;
        for (int i = 0; i < candidate_dummy_children.size(); i++)
        {
            if (candidate_dummy_children[i] != states_under_dummy[i])
            {
                is_match = false;
                break;
            }
        }
        if (is_match)
        {
            VertexId root_id = this->graph_.GetRoot();
            std::unordered_map<VertexId, double> weights = this->graph_.GetOutEdgeWeights(root_id,
                                                                                                        false);
            weights[candidate_dummy_node] += 1;
            // std::cerr<<"all the weights to be updated in merge new parsing with graph: \n";
            // for(auto iter = weights.begin(); iter != weights.end(); ++iter)
            // {
            //     std::cerr<<iter->second<<" ";
            //     if(iter->second == 0)
            //     {
            //         std::cerr<<"update weight in mergenewparsingwithgraph error!\n";
            //         throw std::exception();
            //     }
            // }
            // std::cerr<<std::endl;
            this->graph_.SetOutEdgeWeights(root_id, weights);
            std::cout << "Exact match in parsing\n";
            return states_under_dummy;
        }
    }

    // // case when there is no match in previous parse, connect the partial parse to the root.
    // for (Symbolic_Rule<StateType> rule: best_parsing)
    // {
    //     /** std::cerr << "AddRule from: " << rule.GetSource().GetContent() << " TO ";
    //     for (auto result : rule.GetResults())
    //         std::cerr << result.GetContent() << " ";
    //     std::cerr << std::endl; */
       
    //     this->graph_.AddRule(rule);
        

    // }

    //reset the first data's weight after changing the root to or node
    VertexId root_id = this->graph_.GetRoot();
    root_state = this->graph_.GetStateByVertexId(root_id);
    Symbolic_Rule<StateType> rule_from_root_to_new_parsing(root_state, states_under_dummy);
    this->graph_.AddRule(rule_from_root_to_new_parsing);
    
    //if it is the first data, make the root to be an or-node by adding and deleting one dummy rule
    if(this->graph_.GetVertexContent(root_id)->IsAnd())
    {
        Symbolic_State<StateType> dummy_state(-2);
        Symbolic_Rule<StateType> dummy_rule(root_state,{dummy_state});
        this->graph_.AddRule(dummy_rule);
        this->graph_.DeleteRule(dummy_rule);
    }
    
			


    //chage the weight of the first rule, only occur when root has two data after merging
    // if(change_weight)
    // {
    //     std::unordered_map<VertexId, double> new_weights;
    //     VertexId curr_root_id = this->graph_.GetRoot();
    //     std::vector<VertexId> dummys = this->graph_.ChildrenVertices(curr_root_id);
    //     //check if my implementation is correct
        
    //     std::cerr << dummys.size()<<std::endl;
    //     assert(dummys.size() == 2);
    //     if(dummys[0] == prev_root_id)
    //         new_weights.insert({{dummys[0],count},{dummys[1],1}});
    //     else   
    //         new_weights.insert({{dummys[0],1},{dummys[1],count}});

    //     this->graph_.SetOutEdgeWeights(curr_root_id,new_weights);  
    // }

    /** std::cerr << "AddRule from: " << rule_from_root_to_new_parsing.GetSource().GetContent() << " TO ";
    for (auto result : rule_from_root_to_new_parsing.GetResults())
        std::cerr << result.GetContent() << " ";
    std::cerr << std::endl;

    std::cerr << "MergeNewParsingWithGraph: No matching dummy node in previous parsing, next level rule number is:" << states_under_dummy.size() << std::endl; 
    for (auto state: states_under_dummy)
        std::cerr << state.GetContent() << " ";
    std::cerr << std::endl;     */
//    std::cout << "after:\n";
//    for (auto state: states_under_dummy)
//        std::cout << "(" << state.GetContent() << "," << state.GetId() << ")" << " ";
//    std::cout << "\n";
    // for (int i = 0; i < states_under_dummy.size(); i++)
    //     std::cout << "(" << states_under_dummy[i].GetContent() << "," << states_under_dummy[i].GetId() << ") ";

    // std::cout << std::endl;
    return states_under_dummy;
}

template<class StateType>
std::vector<Symbolic_State<StateType> > Online_Learner<StateType>::Checkpoint::GetBestParsing(double &best_probabability)
{
    std::vector<Symbolic_State<StateType> > best_parsing;
    // std::vector<typename StateList<StateType>::statelist > : each chart
    /** std::cerr << "-------------------GetBestParsing()-------------------" << std::endl; */
    std::vector<std::vector<typename StateList<StateType>::statelist> > all_partial_rules = this->parser_->GetPartialParse();
    /** std::cerr << "GetBestParsing():size of all partial rules " << all_partial_rules.size() << std::endl; */
    
    // for each run from earley parser
    // @ param chart: parsed chart for each run
    for (std::vector<typename StateList<StateType>::statelist> chart : all_partial_rules)
    {
        // best_parsing_in_chart: best parsing in each run of the parser
        Symbolic_State<StateType> best_parsing_start;
        double best_prob = -1;
        int last_statelist_index = chart.size() - 1;
        // find the top most root in the last statelist (i.e. all possible parsing for this seq)
        for (state<StateType> st : chart[last_statelist_index])
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
                double prob = 1;
                
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
                        auto pointed_state = chart[back_pointer.first][back_pointer.second];
                        /** std::cerr << "Back Pointer points to: [" << pointed_state.i_ << "," << pointed_state.j_ << "], FROM " << pointed_state.rule_->left().GetContent() << " TO ";
                        auto results =  pointed_state.rule_->right()[0];
                        for (auto result : results)
                        {
                            std::cerr << result.GetContent() << " ";
                        }
                        std::cerr << std::endl; */
                        q.push(chart[back_pointer.first][back_pointer.second]);
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
                    //     /** std::cerr<<"size of states in AOG: "<<all_states.size()<<std::endl; */
                    //     parentId = this->graph_.GetVertexIdByState(st.rule_->left());
                        
                    // }
                    // catch(std::exception e)
                    // {
                        
                        // for(auto state : all_states)
                        //     std::cerr<<state.GetId()<<std::endl;
                        
                        // AOFStruct<StateType> peeker;
                        /** std::cerr << "AOF id_: " << peeker.id_; */
                        // auto all_states = this->graph_.GetStates();
                        /** std::cerr<<"size of states in AOG: "<<all_states.size()<<std::endl;
                        for(auto state : all_states)
                        {
                            std::cerr<<"State Content: "<<state.GetContent()<<"_\n"<<"State ID: "<<state.GetId()<<std::endl;
                        }

                        std::cout << e.what() << std::endl; */
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
                                double weight = outEdgeWeights[id];
                                if(weight == 0)
                                {
                                    std::cerr<<"the parent vertices: parent vertice: "<<parentId<<std::endl;
                                    auto parent_state = this->graph_.GetStateByVertexId(parentId);
                                    std::cerr<<"the parent vertice state: "<<"("<< parent_state.GetContent() <<", "<<parent_state.GetId()<<")"<< std::endl;

                                    std::cerr<<"the rule's target id: \n";
                                    for(auto target : right_hand_state_ids)
                                        std::cerr<<target<<" ";
                                    std::cerr << std::endl;
                                    
                                    std::cerr<<"the rule's target state: \n";
                                    for(auto target : right_hand_state_ids)
                                        std::cerr<<"("<<this->graph_.GetStateByVertexId(target).GetContent()<<", "<< this->graph_.GetStateByVertexId(target).GetId()<<")";
                                    std::cerr << std::endl;
                                    // this->graph_.OutputGraph("weight_0",PATH,true, SIMPLIFY);
                                    std::cerr<<"weight 0 found!\n";
                                    throw std::exception();
                                }
                                prob *= weight;
                                // std::cout << "weight is: " << weight << ", probability is: " << prob << std::endl;
                                break;
                            }
                        }
                    } 
                }

                // else find the largest probability parsing
                if (prob > best_prob)
                {
                    best_prob = prob;
                    best_parsing_start = possible_parsing_start;
                }
            }
        }
        // std::cout << "Best Prob is: " << best_prob << "\n";
        best_probabability *= best_prob;
        assert((best_parsing_start.GetId() == -1 && best_parsing_start.GetContent() != "") || best_parsing_start.GetId() != -1);
        best_parsing.push_back(best_parsing_start);
    }
    // std::cerr << "Print the content in the best_parsing: \n";
    // for (auto it = best_parsing.begin(); it != best_parsing.end(); it++)
    // {
    //     std::cerr << it->GetId() << " ";
    // }

    return best_parsing;
}
