//
// Created by yuanluyao on 10/23/17.
//

#ifndef AOG_LIB_T_AOG_H
#define AOG_LIB_T_AOG_H

#include <memory>
#include <random>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <stack>
#include <utility>
#include <algorithm>

#include "../Core/Graph.hpp"
#include "AOG_Vertex.h"
#include "AOG_Edge.h"
#include "Symbolic_Rule.h"

template <class T>
using AOG_Graph = AOG_LIB::Graph<AOG_LIB::AOG_Vertex<T>, AOG_LIB::AOG_Edge>;

template <class T> class Online_Learner; //forward declaration
namespace AOG_LIB
{
    /* This class defines an AOG graph with basic operations
     * An AOG graph consists of And-nodes and Or-nodes with their own semantic meanings and edges with different weights.
     * Each vertex in the AOG graph has its own unique id.
     */
    template<class StateType>
    class T_AOG : public AOG_Graph<StateType>
    {
        friend class Online_Learner<StateType>;
        // all the rules stored in a hash table for easy look up
        std::unordered_set<Symbolic_Rule<StateType> > all_rules_;

        //all leaf states stored in a hash table for easy look up
        std::unordered_set<Symbolic_State<StateType> > all_leaf_states_;

        //a map from state to graph vertex id
        std::unordered_map<Symbolic_State<StateType>, VertexId> state_to_vertex_;
        // the vertex id of the root node
        VertexId root_id_;
        //an auxiliary boolean variable for error checking when adding a vertex
        bool has_root_;

    public:
        //default constructor
        T_AOG();
    
        //construct from a set of rules. Currently it assigns the weight of 1.0 to all edges
        T_AOG(const std::vector<Symbolic_Rule<StateType> >&);

        //T_AOG copy constructor
        T_AOG(const T_AOG &other);

        //construct from a set of leaf states
        explicit T_AOG(const std::vector<Symbolic_State<StateType> > &);
        //return the number of rules the T_AOG contains
        unsigned NumOfRules() const;
        //return the number of symbolic states in the T_AOG
        unsigned NumOfStates() const;
        //return the number of leaf states in the T_AOG
        unsigned NumOfLeafStates() const;

        double GrammarSize(double leaf_bias, bool simplified_graph=false);

        bool HasRoot() const{return has_root_;}
        /* This function adds a vertex to the AOG graph
         * @param: 
         *      aog_vertex: the vertex to be added to the AOG graph
         * @throw: 
         *      throw an exception when user tries to add a second root node to the AOG graph
         * @return:
         *      an unsigned integer that indicates the id of the newly added vertex(if the vertex alreadly exists, return the existing id).
         */    
        virtual VertexId AddVertex(const std::shared_ptr<AOG_Vertex<StateType> >);

        /* This function adds an edge between two vertices.      
         * Cannot add an edge that starts and ends at same vertex
         * @params:
         *     source: the vertex id of the source node
         *     target: the vertex id of the target node
         *     aog_edge: the edge to be added from the source node to the target node
         *     multi_edge: true if the AOG graph allows multiple edges between two vertices
         * @return:
         *     true if the edge is added successfully
         */        
        virtual bool AddEdge(const VertexId, const VertexId, const std::shared_ptr<AOG_Edge>,bool=true, int=-1);

        /* This function deletes a vertex and all related edges from the AOG graph
         * @params: 
         *     vid: the id of the vertex to be deleted
         * @return:
         *     true if the vertex is successfully
         */
        virtual bool DeleteVertex(const VertexId);

        /* This function returns all the rules used to construct the AOG graph
         * @param:
         *     void
         * @return:
         *     a vector that contains all the symbolic
         */
        std::vector<Symbolic_Rule<StateType> > GetRules() const;

        /* This function returns all the rules at a certain level. Level 1 returns the root rule.
         * @param:
         *     level: the level user wants to request
         * @return:
         *     a vector of rules at a given level
         */
        std::vector<Symbolic_Rule<StateType> > GetRulesOfLevel(long);

        /* This function returns all the leaf states in the AOG graph
         * @param:
         *     void
         * @return:
         *     a vector that contains all the leaf states
         */
        std::vector<Symbolic_State<StateType> > GetLeafStates() const;

       /* This function returns all the symbolic states in the AOG graph
        * @param:
        *     void
        * @return:
        *     a vector that contains all the states in vertices in AOG graph
        */
        std::vector<Symbolic_State<StateType> > GetStates() const;

        /* @param: 
         *     source: the id of the vertex containing the desired state
         * @return:
         *     a Symbolic_State object that represents the state inside the given vertex
         */
        Symbolic_State<StateType> GetStateByVertexId(VertexId) const;

        /* @param:
         *     void
         * @return:
         *     the vertex id of the root vertex if there exist one, else throw an exception
         */
        VertexId GetRoot() const;

        /* @param:
         *     state: the symbolic state the user want to query about
         * @return:
         *     an unsigned integer that indicates the id of the vertex containing the given symbolic state
         */
        VertexId GetVertexIdByState(const Symbolic_State<StateType> &) const;
    
        /* This function modifies the symbolic state of a given vertex.
         * If the given state already exists in the AOG graph, merge the given vertex with the vertex that already has that state.
         * Else change the state of the given vertex to new state
         * @params:
         *     curr_id: the id of the vertex whose state needs to be changed
         *     state: the symbolic state to be changed to
         * @return:
         *     if the given new state already exists, return the id of the merged vertex, else the id of the vertex passed in
         */
        /* VertexId SetState(VertexId, const Symbolic_State<StateType> &); */

        /* This function changes a given node to an and/or-node
         * @params:
         *      source_id: the id of the vertex to be changed
         *      is_and: true if trying to set the vertex to an and-node
         * @return:
         *      void
         */
        void SetIsAnd(VertexId, const bool);

        /* This function reassigns weights to a given vertex's outedges
         * @params:
         *      source: the source vertex of all the outedges to be modified
         *      weight: an unordered_map that contains the mapping from the id of outedges to the weights each outedge need to be set to
         * @return:
         *      true if the weights are set successfully
         */
        bool SetOutEdgeWeights(VertexId, const std::unordered_map<VertexId ,double>&);

        /* This function gets weights of all the outedges from a given source vertex
         * @params:
         *      source: the id of the source vertex
         *      is_normalized: a boolean value. If true, return the normalized weights of all the outedges from the target vertex.
         * @return:
         *      return an unordered_map that maps all the outedges' target vertices' ids to their corresponding weights
         */
        std::unordered_map<VertexId,double> GetOutEdgeWeights(VertexId,bool) const;

        /* This function normalizes all the weights from a given source vertex
         * @param:
         *      src_id: the id of the source vertex whose outedges will be normalized
         * @return:
         *      an unordered_map that maps ids of all target vertices of outedges to the normalized weights.
         */
        std::unordered_map<VertexId,double> Normalize(VertexId);

        /* This function checks a certain rule for its existence
         * A->BC and A->CB are considered different rules
         * @param:
         *      rule: the rule to be checked
         * @return:
         *      true if the rule already exists in T_AOG, false otherwise
         */
        bool ExistRule(const Symbolic_Rule<StateType>&);

        /* This function adds a rule to the T_AOG
         * @param:
         *      rule: the rule to be added
         * @return:
         *      bool: whether the add is successful
         */
        bool AddRule(const Symbolic_Rule<StateType>&);

        /* This function deletes corresponding edges and vertices in the AOG graph given the rule to be deleted
         * @param:
         *      rule: the rule to be deleted
         * @return:
         *      false if the given rule does not exist, true otherwise
         */
        bool DeleteRule(const Symbolic_Rule<StateType>&);

        /* This function takes an And-Or fragment represented by a vector of rules, and searches for parts that can
         * induce the same leaves as the fragment. If any is found, replace it with the fragment
         * The second parameter is a map indicating which rule should be replaced by which configuration
         * @param:
         *      fragment: a vector of rules representing the fragment
         *      match: a map indicating which configuration replaces which rule
         * @return:
         *      void
         */
        void ReplaceFragment(const std::vector<Symbolic_Rule<StateType> >&,
                             const std::unordered_map<Symbolic_State<StateType>, std::vector<Symbolic_State<StateType> > >&);

        /* This Sample() function eliminates all Or-nodes by sampling over edge weights
         * @param:
         *      root: The vertex to sample from
         * @return:
         *      an unordered map containing parent VertexId (key) and a vector of child VertexId (value)
         */
        std::shared_ptr<std::vector<VertexId>> Sample(VertexId, std::vector<VertexId> &, double &) const;

        /* This Visualize function outputs a file that can be used for visualization
         * @param:
         *      filename: the name of the visualization input file; if not given, default is "aog.txt"
         *      dir: the target directory; if not given, default is the current working directory
         * @return:
         *      void
         * @output:
         *      a text file containing the following information:
         *
         *          x -> y[weight="z"]: an edge pointing from node x to node y with weight z
         *          x[label="x", name="x", AO="O"]: an Or-node with id x
         *          y[label="y", name="y", AO="A"]: an And-node with id y
         *          z[label="z", name="z", AO="X"]: an auxiliary root generated automatically when required
         *
         *      Example:
         *
         *          digraph g{
         *              0 -> 1[weight="0.4"];
         *              0 -> 2[weight="0.6"];
         *              0[label="0", name="0", AO="O"];
         *              1[label="1", name="1", AO="A"];
         *              2[label="2", name="2", AO="A"];
         *          }
         */
        void Visualize(std::string="aog", std::string="");

        /* This GetTopLevelStates() function returns the top most level state(s) of the T-AOG 
         * @param:
         *      void
         * @return:
         *      a vector containing the top most level state(s) of the T-AOG 
         */
        std::vector<Symbolic_State<StateType> > GetTopLevelStates();

        /* This GetTopLevelRules() function returns the top most level rule(s) of the T-AOG 
         * @param:
         *      void
         * @return:
         *      a vector containing the top most level rule(s) of the T-AOG 
         */
        std::vector<Symbolic_Rule<StateType> > GetTopLevelRules();
        
        void DeleteNoParentRules(const Symbolic_State<StateType>& state);

        void OutputGraph(std::string="graph", std::string="", bool=true, bool=true);

        void TruncateGraph();

        void SimplifyGraph();

        void OutputLearnedTree(std::string path, std::string appendix="");

        virtual unsigned NumberOfVertices() const
        {return state_to_vertex_.size();}

        void SetRoot(Symbolic_State<StateType> root);
        void PrintOrNodeAndChildrenCount(std::ofstream &file);
    };

    template <class StateType>
    void T_AOG<StateType>::SetRoot(Symbolic_State<StateType> root)
    {
        this->has_root_ = true;
        this->root_id_ = this->GetVertexIdByState(root);
    }

    template<class StateType>
    T_AOG<StateType>::T_AOG()
            :AOG_Graph<StateType>()
    { this->has_root_ = false; }
    
    //copy constructor for T_AOG
    template<class StateType>
    T_AOG<StateType>::T_AOG(const T_AOG &other)
            :AOG_Graph<StateType>(other), 
            all_rules_(other.all_rules_),
            all_leaf_states_(other.all_leaf_states_),
            state_to_vertex_(other.state_to_vertex_),
            root_id_(other.root_id_),
            has_root_(other.has_root_)
    { }

    template<class StateType>
    T_AOG<StateType>::T_AOG(const std::vector<Symbolic_Rule<StateType> > &rules)
            :AOG_Graph<StateType>(0, true)
    {
        has_root_ = false;
        for (auto rule: rules)
            this->AddRule(rule);
    }

    template<class StateType>
    T_AOG<StateType>::T_AOG(const std::vector<Symbolic_State<StateType> > &leaf_states)
            :AOG_Graph<StateType>(0, true)
    {
        //construct a root node
        has_root_ = false;
        std::shared_ptr<Symbolic_State<StateType> > root_state = std::make_shared<Symbolic_State<StateType> >();
        std::shared_ptr<AOG_Vertex<StateType> > root = std::make_shared<AOG_Vertex<StateType> >(*root_state, false, true);
        this->root_id_ = this->AddVertex(root);
        this->state_to_vertex_[*root_state] = this->root_id_;
        has_root_ = true;

        for (unsigned i = 0; i < leaf_states.size(); i++)
        {
            if (!leaf_states[i].GetIsBasic())
            {
                std::cerr << "Non basic states passed in\n";
                throw std::exception();
            }
            //duplicate states in the sequence
            if (this->all_leaf_states_.find(leaf_states[i]) != this->all_leaf_states_.end())
                continue;
            //take in leaf states as content, not an and-node, not a root node
            std::shared_ptr<AOG_Vertex<StateType> > ptr = std::make_shared<AOG_Vertex<StateType> >(leaf_states[i], false, false);
            VertexId new_v_id = this->AddVertex(ptr);
            this->all_leaf_states_.insert(leaf_states[i]);
            this->state_to_vertex_[leaf_states[i]] = new_v_id;
        }
    }

    template<class StateType>
    unsigned T_AOG<StateType>::NumOfRules() const
    { return unsigned(this->all_rules_.size()); }

    template<class StateType>
    unsigned T_AOG<StateType>::NumOfStates() const
    { return unsigned(this->state_to_vertex_.size()); }

    template<class StateType>
    unsigned T_AOG<StateType>::NumOfLeafStates() const
    { return unsigned(this->all_leaf_states_.size()); }



    template<class StateType>
    std::vector<Symbolic_Rule<StateType> > T_AOG<StateType>::GetRules() const
    {
        std::vector<Symbolic_Rule<StateType> > all_rules;
        for (auto rule: this->all_rules_)
            all_rules.push_back(rule);
        return all_rules;
    }


    template<class StateType>
    std::vector<Symbolic_Rule<StateType> > T_AOG<StateType>::GetRulesOfLevel(long level)
    {
        std::unordered_multimap<Symbolic_State<StateType>, Symbolic_Rule<StateType> > map;
        for(auto rule : all_rules_)
            map.insert(std::make_pair(rule.GetSource(), rule));

        int roll = 0;
        std::vector<Symbolic_State<StateType> > tmp;
        std::vector<std::vector<Symbolic_State<StateType> > > sources;
        sources.push_back(tmp), sources.push_back(tmp);

        //find the roots first
        VertexId root_id = this->GetRoot();
        sources[roll].push_back(this->GetStateByVertexId(root_id));

        //put
        std::vector<Symbolic_Rule<StateType> > rules;
        --level;
        while(level > 0){
            // @param src: iterator of vector<Symbolic_State<StateType> >
            for(auto src = sources[roll].begin(); src != sources[roll].end();++src){
                std::cerr << src->GetContent() << std::endl;
                // @param range: iterator of map<Symbolic_State<StateType>, Symbolic_Rule<StateType> >
                auto range = map.equal_range(*src);

                for(auto iter = range.first; iter != range.second;++iter)
                    for(auto s : iter->second.GetResults())
                        sources[1-roll].push_back(s);
                
            }
            --level;        
            sources[roll].clear();
            roll = 1-roll;
            
        }
        
        //find rules from the states
        for(auto src : sources[roll])
        {
            auto range = map.equal_range(src);
            for(auto iter = range.first; iter != range.second;++iter)
                rules.push_back(iter->second);
        }

        return rules;
    }

    template<class StateType>
    std::vector<Symbolic_State<StateType> > T_AOG<StateType>::GetLeafStates() const
    {
        std::vector<Symbolic_State<StateType> > all_states(this->all_leaf_states_.size());
        for (auto state: this->all_leaf_states_)
            all_states.push_back(state);
        return all_states;
    }

    template<class StateType>
    std::vector<Symbolic_State<StateType> > T_AOG<StateType>::GetStates() const
    {
        std::vector<Symbolic_State<StateType> > all_states(0);
        for (auto state: this->state_to_vertex_)
            all_states.push_back(state.first);
        return all_states;
    }

    template<class StateType>
    VertexId T_AOG<StateType>::GetRoot() const 
    {
        if (!this->has_root_){
            std::cerr << "The TAOG does not have a root\n";
            throw std::exception();
        }
        return this->root_id_;
    }

    template<class StateType>
    Symbolic_State<StateType> T_AOG<StateType>::GetStateByVertexId(const VertexId vid) const
    {
        if (this->IsValidVertex(vid))
            return this->GetVertexContent(vid)->GetState();
        std::cerr<<"Invalid State\n";
        throw std::exception();
    }

    template<class StateType>
    VertexId T_AOG<StateType>::GetVertexIdByState(const Symbolic_State<StateType> &state) const
    {
        auto iter = this->state_to_vertex_.find(state);
        if (iter == this->state_to_vertex_.end())
        {
            std::cerr << "Invalid State\n";
            throw std::exception();
        }
        return iter->second;
    }

    //T-AOG version add vertex, avoid adding duplicate vertices if they have same state
    template<class StateType>
      VertexId T_AOG<StateType>::AddVertex(const std::shared_ptr<AOG_Vertex<StateType> > aog_vertex)
    {
        Symbolic_State<StateType> state = aog_vertex->GetState();
        bool is_and = aog_vertex->IsAnd();
        bool is_root = aog_vertex->IsRoot();
        //check if user wants to add a second root
        if(is_root && has_root_)
        {
            std::cerr<<"Cannot add a root vertex.\n";
            throw std::exception();
        }
        //if the symbolic state in the given vertex does not exist, add a new vertex
        auto iter = this->state_to_vertex_.find(state);
        
        if (iter == this->state_to_vertex_.end())
        {
            std::shared_ptr<AOG_Vertex<StateType> > ptr(new AOG_Vertex<StateType>(state,is_and,is_root));
            VertexId new_id = AOG_Graph<StateType>::AddVertex(ptr);
            this->state_to_vertex_[state] = new_id;
            if (state.GetIsBasic())
                this->all_leaf_states_.insert(ptr->GetState());
            
            //if the vertex added want to be a root, add this vertex as root
            if(is_root && !(this->has_root_))
            {
                this->has_root_ = true;
                this->root_id_ = new_id;
            }
            return new_id;
        }

        //else return the existing vertex's id
        if(is_root && !(this->has_root_))
        {
            this->has_root_ = true;
            this->root_id_ = iter->second;
        }
        return iter->second;
    }

    template<class StateType>
    bool T_AOG<StateType>::DeleteVertex(const VertexId vid)
    {
        //update the data structures that keep track of all states and leaf states
        if (this->IsValidVertex(vid))
        {
            std::shared_ptr<AOG_Vertex<StateType> > aog_vertex = this->GetVertexContent(vid);
            Symbolic_State<StateType> state = aog_vertex->GetState();
            if (state.GetIsBasic())
                this->all_leaf_states_.erase(state);
            this->state_to_vertex_.erase(state);
        }

        //Invalid situation handled in base class
        return AOG_Graph<StateType>::DeleteVertex(vid);
    }


    template<class StateType>
    bool T_AOG<StateType>::AddEdge(const VertexId source, const VertexId target,
                                    const std::shared_ptr<AOG_Edge> aog_edge, bool multi_edge,int pos)
    {
        double weight = aog_edge->GetWeight();
        //bound checking
        if(weight < 0)
        {
            std::cerr << "Weight cannot be less than 0!\n";
            return false;
        }
        if(this->IsValidVertex(source))
        {
            //source of the edge to add should not be a leaf state
            Symbolic_State<StateType> src_state = GetStateByVertexId(source);
            if(src_state.GetIsBasic())
            {
                std::cerr<<"cannot add edge starting from leaf state\n";
                return false;
            }
        }

        //invalid situations handled in base class
        std::shared_ptr<AOG_Edge> ptr = std::make_shared<AOG_Edge>(weight);
        return AOG_Graph<StateType>::AddEdge(source, target, ptr, multi_edge, pos);
    }

    template<class StateType>
    bool T_AOG<StateType>::SetOutEdgeWeights(VertexId source, const std::unordered_map<VertexId, double> & weights)
    {

        // std::cout << "inside setoutedgeweights............." << std::endl;
        // for (auto i : weights){
        //     std::cout << "id :" << i.first << " weight: " << i.second << std::endl;
        // }
        if(this->GetVertexContent(source)->IsAnd())
        {
            std::cerr<<"Cannot get weights of edges from an And-node\n";
            return false;
        }
        //check if the number of weights inputs is the same as the number of out edges
        if(weights.size() != this->OutEdges(source).size())
        {
            std::cerr<<"The number of weights inputs is inconsistent with the number of out edges.\n";
            return false;
        }
        //update weights
        for(auto iter:weights){
            // if(iter.second == 0){
            //     std::cerr << "should not set weight to 0!" << std::endl;
            //     std::cerr << "(" << this->GetStateByVertexId(source).GetContent() << ", " << this->GetStateByVertexId(source).GetId() << ") ->" << std::endl;
            //     for(const auto& it : weights){
            //         std::cerr << "\t(" << this->GetStateByVertexId(it.first).GetContent() << ", " << this->GetStateByVertexId(it.first).GetId() << "), weight: " << it.second << std::endl;
            //     }
            //     throw std::exception();
            // }
            this->GetEdgeContent(source, iter.first)[0]->SetWeight(iter.second);
        }

        return true;
    }

    template<class StateType>
    std::unordered_map<VertexId, double> T_AOG<StateType>::GetOutEdgeWeights(VertexId source, bool is_normalized) const
    {

        //report error if trying to get weights of an and-node's edges
        if(this->GetVertexContent(source)->IsAnd())
        {
            std::cerr<<"Cannot get weights of edges from an And-node\n";
            throw std::exception();
        }

        std::unordered_map<VertexId, double> weights;
        double total_weights = 0.0;
        for(auto edge:this->OutEdges(source))
        {
            double weight = this->GetEdgeContent(source, edge.second)[0]->GetWeight();
            total_weights += weight;
            weights[edge.second] = weight;
        }
        // std::cout << "inside getoutedgeweights............." << std::endl;
        // for (auto i : weights){
        //     std::cout << "id :" << i.first << " weight: " << i.second << std::endl;
        // }

        //if user wants unnormalized weights
        if(!is_normalized)
            return weights;

        double coeff = 1.0/total_weights;
        for(auto &iter : weights)
            iter.second *= coeff;

        return weights;
    }


    /* template<class StateType> */
    /* VertexId T_AOG<StateType>::SetState(VertexId curr_vid, const Symbolic_State<StateType> &state) */
    /* { */
    /*     auto iter = state_to_vertex_.find(state); */
    /*     //if the new state does not exist in the T_AOG, add it */
    /*     if(iter == state_to_vertex_.end()) */
    /*     { */
    /*         this->GetVertexContent(curr_vid)->SetState(state); */
    /*         return curr_vid; */
    /*     } */

    /*     //does not allow to change the state of a dummy node */
    /*     if(!this->ParentsVertices(curr_vid)[0]->IsAnd()) */
    /*       { */
    /*         std::cerr<<"Cannot change the state of a dummy vertex!"<<std::endl; */
    /*         return curr_vid; */
    /*       } */
    /*     //if the state already exists, merge it with the original one */
    /*     VertexId ori_vid = iter->second; */

    /*     //find all vertices connected to the current vertex */
    /*     std::vector<VertexId> children = this->ChildrenVertices(curr_vid); */
    /*     std::vector<VertexId> parents = this->ParentsVertices(curr_vid); */

    /*     //used to keep track of duplicate children/parents */
    /*     std::unordered_set<VertexId> added_children; */
    /*     std::unordered_set<VertexId> added_parents; */
    /*     //store contents of all edges connected to the current vertex */
    /*     std::vector<std::shared_ptr<AOG_Edge> > children_contents; */
    /*     std::vector<std::shared_ptr<AOG_Edge> > parents_contents; */
    /*     for(auto child : children) */
    /*     { */
    /*         if(added_children.find(child) == added_children.end()) */
    /*         { */
    /*             added_children.insert(child); */
    /*             children_contents.push_back(std::make_shared<AOG_Edge>(this->GetEdgeContent(curr_vid, child)[0])); */
    /*         } */
    /*         //if one edge of the child has already added, add another edge */
    /*         else */
    /*             children_contents.push_back(std::make_shared<AOG_Edge>(this->GetEdgeContent(curr_vid, child)[1])); */
    /*     } */
            

    /*     for(auto parent : parents) */
    /*     { */
    /*         if(added_parents.find(parent) == added_parents.end()) */
    /*         { */
    /*             added_parents.insert(parent); */
    /*             parents_contents.push_back(std::make_shared<AOG_Edge>(this->GetEdgeContent(parent, curr_vid)[0])); */
    /*         } */
    /*         else */
    /*             parents_contents.push_back(std::make_shared<AOG_Edge>(this->GetEdgeContent(parent, curr_vid)[1])); */
    /*     } */


    /*     //delete the current vertex, along with its edges */
    /*     if(!this->DeleteVertex(curr_vid)) */
    /*         throw std::exception(); */

    
    /*     //TODO: and-> and node: if children are the same as the children of the state merged to, dont add, else chagne the state merged to to or-node */
    /*     //TODO: and -> or node: add a new dummy node, but what about the branching probability? what if the rule already exists?     */
    /*     //TODO: or -> and node: change the node merged to to Or-node, what about the branching probability of the original children? */
    /*     //TODO: or -> or node: normalize the probability? what if one rule already exist? */
    
    /*     //or -> or node */
    /*     // or -> and node */
    

    /*     //TODO: update rules */
    
    /*     //merge in/out edges to the original vertex */
    /*     for(int i = 0; i < children.size();i++) */
    /*       this->AddEdge(ori_vid, children[i], children_contents[i],true); */

    /*     for(int i = 0; i < parents.size();i++) */
    /*       this->AddEdge(parents[i], ori_vid, parents_contents[i],true); */

    /*     return ori_vid; */
    /* } */

    template<class StateType>
    void T_AOG<StateType>::SetIsAnd(VertexId source_id, const bool is_and)
    {
        if (this->IsValidVertex(source_id))
            this->GetVertexContent(source_id)->SetIsAnd(is_and);
        else
            std::cerr << "Vertex" << source_id << " doesn't exist.\n";

        return;
    }

    template<class StateType>
    std::unordered_map<VertexId,double> T_AOG<StateType>::Normalize(VertexId src_id)
    {
        //if it is an And-node, do nothing
        if(this->GetVertexContent(src_id)->IsAnd())
        {
            std::cerr << "Cannot normalize edges from And-node\n";
            throw std::exception();
        }

        //get normalized weights
        std::unordered_map<VertexId,double> weights = this->GetOutEdgeWeights(src_id,true);
        //set weights
        this->SetOutEdgeWeights(src_id,weights);
        return weights;
    }
    
    template<class StateType>
    bool T_AOG<StateType>::ExistRule(const Symbolic_Rule<StateType>& rule)
    {
        return this->all_rules_.find(rule) != this->all_rules_.end();
    }

    template<class StateType>
    bool T_AOG<StateType>::AddRule(const Symbolic_Rule<StateType>& rule)
    {
        // construct an edge with weight 1.0
        std::shared_ptr<AOG_Edge> edge = std::make_shared<AOG_Edge>();

        // check if the rule already exists
        if (this->ExistRule(rule))        
            return false;

        // get the source and result states from each rule and construct vertices for each of them.
        std::shared_ptr<Symbolic_State<StateType> > source_state = std::make_shared<Symbolic_State<StateType> >(rule.GetSource());
        std::shared_ptr<AOG_Vertex<StateType> > source = std::make_shared<AOG_Vertex<StateType> >(*source_state, true, false);
        
        // std::cerr<<"this rule is added from state id: "<<(*source_state).GetId()<<std::endl;
        VertexId src_id = this->AddVertex(source);               
        this->state_to_vertex_[*source_state] = src_id;
    
        std::vector<std::shared_ptr<Symbolic_State<StateType> > > result_states;
        std::vector<VertexId> rs_vtxs_id;
        auto results = rule.GetResults();

        for(auto result : results)
        {
            // std::cerr<<"this state is added to state id: "<<result.GetId();
            result_states.push_back(std::make_shared<Symbolic_State<StateType> >(result));
            std::shared_ptr<AOG_Vertex<StateType> > rs_vtx = std::make_shared<AOG_Vertex<StateType> >(result,true,false);
            VertexId rs_id = this->AddVertex(rs_vtx);
            rs_vtxs_id.push_back(rs_id);
            if(result.GetIsBasic())
                this->all_leaf_states_.insert(result);
            this->state_to_vertex_[result] = rs_id;
        }
                       
        std::vector<VertexId> children = AOG_Graph<StateType>::ChildrenVertices(src_id);


        //create a dummy And-node
        std::shared_ptr<Symbolic_State<StateType> > dummy_state = std::make_shared<Symbolic_State<StateType> >(*source_state);
        std::shared_ptr<AOG_Vertex<StateType> > dummy = std::make_shared<AOG_Vertex<StateType> >(*dummy_state, true, false);
        // if src is already a parent
        if (children.size()) {
            auto src_content = source_state->GetContent();

            // Add the dummy And-node for the new rule
            VertexId dummy_id = AOG_Graph<StateType>::AddVertex(dummy);

            // if src is an And-node, it has no dummy nodes yet
            if (this->GetVertexContent(src_id)->IsAnd())
            {
                // create an Or-node new parent between src and its parent
                // std::shared_ptr<Symbolic_State<StateType> > new_par_state
                //         = std::make_shared<Symbolic_State<StateType> >(src_content, false);
                std::shared_ptr<Symbolic_State<StateType> > new_par_state
                        = std::make_shared<Symbolic_State<StateType> >(*source_state);
                std::shared_ptr<AOG_Vertex<StateType> > new_par
                        = std::make_shared<AOG_Vertex<StateType> >(*new_par_state, false, false);
                //delete the old state and add the new state
                state_to_vertex_.erase(*source_state);
                VertexId new_par_id = this->AddVertex(new_par);
                this->state_to_vertex_[*new_par_state] = new_par_id;

                //change the root to the newly created node
                if(this->has_root_ && src_id == this->GetRoot())
                {
                    this->root_id_ = new_par_id;
                    this->GetVertexContent(src_id)->SetIsRoot(false);
                    this->GetVertexContent(new_par_id)->SetIsRoot(true);

                }

                // if src has parents, connect parents to new_par
                // delete parents to src
                for (auto parent_id : this->ParentsVertices(src_id)) 
                {                        
                    //find the source node's positions in all rules whose targets contains source
                    std::vector<VertexId> all_children = this->ChildrenVertices(parent_id);                
                    std::vector<int> result_pos;
                    auto result = std::find(all_children.begin(),all_children.end(),src_id);
                    while(result != all_children.end())
                    {
                        result_pos.push_back(result - all_children.begin());
                        result = std::find(result + 1,all_children.end(),src_id);
                    }
                    
                    // if(result_pos.empty())
                    // {
                    //     std::cerr<<"cannot find the source state!\n";
                    //     std::cerr<<"the rule: \n";
                    //     std::cerr<<"("<<rule.GetSource().GetContent()<<", "<<rule.GetSource().GetId()<<")";
                    //     for(auto state : rule.GetResults())
                    //     {
                    //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
                    //     }
                    //     std::cerr<<std::endl;
                    //     throw std::exception();
                    // }
                    
                    //delete all edges that have given source and target
                    this->DeleteEdge(parent_id, src_id);

                    //add all edges from parent to newly added state
                    for(int pos : result_pos)
                        this->AddEdge(parent_id, new_par_id, edge, true, pos);

                }
                // connect new_par to src and dummy
                this->AddEdge(new_par_id, src_id, edge);
                this->AddEdge(new_par_id, dummy_id, edge);
                // connect dummy to all results
                for(auto rs_id : rs_vtxs_id)
                    this->AddEdge(dummy_id, rs_id, edge);

            }
            // if src is an Or-node, simply add another dummy node
            else
            {
                // connect src to dummy
                this->AddEdge(src_id, dummy_id, edge);
                // connect dummy to all results
                for(auto rs_id : rs_vtxs_id)
                    this->AddEdge(dummy_id, rs_id, edge);
            }
        }

        // if src is not yet a parent, and src is root, add a dummy node below, set root to or_node, then add children under dummy  
        else if(!children.size() && this->GetVertexContent(src_id)->IsRoot()) 
        {
            VertexId dummy_id = AOG_Graph<StateType>::AddVertex(dummy);
            this->AddEdge(src_id,dummy_id,edge);
           
            for(auto rs_id : rs_vtxs_id)
                this->AddEdge(dummy_id, rs_id, edge);
            
            //change root to or-node
            this->GetVertexContent(src_id)->SetIsAnd(false);

        }

        // if src is not yet a parent, and src is not root, just add edges to its children
        else
        {
            // std::cerr<<"the edge is added from: "<<src_id<<std::endl;
            for(auto rs_id : rs_vtxs_id)
            {
                
                // std::cerr<<"the edge is added to : "<<rs_id<<std::endl;
                this->AddEdge(src_id, rs_id, edge);
                
            }
            
        }

        this->all_rules_.insert(rule);
        return true;
    }
    
    template<class StateType>
    bool T_AOG<StateType>::DeleteRule(const Symbolic_Rule<StateType>& rule)
    {
      //delete the rule in all_rules_
        if(!all_rules_.erase(rule))
        {
            //std::cerr<<"The rule to be deleted does not exist in the T_AOG!\n";
            return false;
        }
    
        //get vertices ids
        VertexId src_vtx = GetVertexIdByState(rule.GetSource());
        auto results = rule.GetResults();
        std:: vector<VertexId> end_vtxs;
        for (auto result : results)
            end_vtxs.push_back(GetVertexIdByState(result));
    
        size_t end_vtxs_size = end_vtxs.size();
    

        //if the source node is an and-node
        if(this->GetVertexContent(src_vtx)->IsAnd())
        {
            //delete the edges constructed from the rule
            for(auto end_vtx : end_vtxs)
                this->DeleteEdge(src_vtx,end_vtx);
        }

        //if source state is an or-node
        else
        {
            //find the dummy and delete it
            auto dummys = this->ChildrenVertices(src_vtx);
            bool found_rule;
            for(auto dummy : dummys)
            {
                found_rule = true;
                auto dummy_children = this->ChildrenVertices(dummy);
        
                if(dummy_children.size() != end_vtxs_size)
                    continue;
                //if the child is the dummy node that associated with the rule, delete the dummy
                for(size_t i = 0; i< end_vtxs_size;i++)
                    if(dummy_children[i] != end_vtxs[i])
                    {
                        found_rule = false;
                        break;
                    }
                if(found_rule)
                {
                    AOG_Graph<StateType>::DeleteVertex(dummy);
                    break;
                }
            }

            //if there is only one dummy left, and the source is not root, change the source vertex to an and-node
            dummys = this->ChildrenVertices(src_vtx);
            if(dummys.size() == 1 && (!this->GetVertexContent(src_vtx)->IsRoot()))
            {
                std::vector<VertexId> dummy_children = this -> ChildrenVertices(dummys[0]);
                for(auto dummy_child : dummy_children)
                    AddEdge(src_vtx,dummy_child,std::make_shared<AOG_Edge>());
                //delete the dummy node
                AOG_Graph<StateType>::DeleteVertex(dummys[0]);
        
                //change the property of the source vertex to and-node
                this->GetVertexContent(src_vtx)->SetIsAnd(true);
            }
        }
        return true;
    }

    template<class StateType>
    void T_AOG<StateType>::ReplaceFragment
            (const std::vector<Symbolic_Rule<StateType> >& fragment,
             const std::unordered_map<Symbolic_State<StateType>, std::vector<Symbolic_State<StateType> > >& match)
    {
        //add fragment into T_AOG
        for(auto f : fragment) {
            this->AddRule(f);
        }
        
        VertexId frag_root = this->GetVertexIdByState(fragment[0].GetSource());
        std::shared_ptr<AOG_Edge> edge = std::make_shared<AOG_Edge>();

        //each entry in match is the source of a rule, and the corresponding configuration
        //connect teh source to the fragment, and delete all edge from source to each state
        //of the configuration
        for(auto m : match){
            VertexId src_id = this->GetVertexIdByState(m.first);
            this->AddEdge(src_id, frag_root, edge);
            for(auto s : m.second)
                this->DeleteEdge(src_id, this->GetVertexIdByState(s));
        }

        return;
    }
  
    template<class StateType>
    std::shared_ptr<std::vector<VertexId>> T_AOG<StateType>::Sample(VertexId root, std::vector<VertexId>& res, double &prob) const
    {
        prob = 1;
        //if the node to sample from is a dummy node, start sampling from its parent
        auto root_parent = this->ParentsVertices(root);
        if(!root_parent.empty() && !this->GetVertexContent(root_parent[0])->IsAnd())
            root = root_parent[0];
        // the resultant graph that contains the key-value pairs
        // std::unordered_map<VertexId, std::vector<VertexId> > sample_graph;
        std::shared_ptr<std::vector<VertexId>> parse_tree = std::make_shared<std::vector<VertexId> >();
        // for level-order tree traversal
        std::stack<VertexId> stack;
        // will be used to obtain a seed for the random number engine
        std::random_device rd;
        // standard mersenne_twister_engine seeded with rd()
        std::mt19937 gen(rd());

        stack.push(root);
        while (!stack.empty())
        {
            VertexId parent_id = stack.top();
            stack.pop();
            parse_tree->push_back(parent_id);
            std::vector<VertexId> children = this->ChildrenVertices(parent_id);

            if (children.size() == 0){
                // file << this->GetStateByVertexId(parent_id).GetContent() << " ";
                res.push_back(parent_id);
            }
            
            // if the node is an And-node, keep all its children in sample_graph
            else if (AOG_Graph<StateType>::GetVertexContent(parent_id)->IsAnd())
            {
                // sample_graph.insert(std::make_pair(parent_id, children));
                for (int i = children.size()-1; i >= 0; i--)
                    stack.push(children[i]);
            }
            // if the node is an Or-node
            else
            {
                std::unordered_map<VertexId, double> edge_weights = this->GetOutEdgeWeights(parent_id, true);
                
                std::vector<VertexId> vertex;
                std::vector<double> weight;

                for (auto w:edge_weights)
                {
                    vertex.push_back(w.first);
                    weight.push_back(w.second);
                }

                // choose a child under the weight distribution and keep it in sample_graph
                std::discrete_distribution<> dis(weight.begin(), weight.end());
                unsigned index = dis(gen);
                VertexId sample = vertex[index];
                prob *= weight[index];
                
                // std::vector<VertexId> v = {sample};
                // sample_graph.insert(std::make_pair(parent_id, v));
                stack.push(sample);
            }
        }

        return parse_tree;
    }

 
    template<class StateType>
    void T_AOG<StateType>::Visualize(std::string filename, std::string dir)
    {
        std::ofstream file;
        if (dir.empty())
            dir = filename.append(".txt");
        else
            dir = dir.append("/").append(filename.append(".txt"));

        file.open(dir);
        if (file.is_open()) {
            std::unordered_map<VertexId, std::vector<VertexId> > map = this->PrintGraphJSON();

            file << "digraph g{" << std::endl;
            for (auto entry : map )
            {
                for (auto target : entry.second)
                {
                    std::vector<std::shared_ptr<AOG_Edge> > edge = this->GetEdgeContent(entry.first, target);
                    for (auto e : edge)
                        file << entry.first << " -> " << target << "[weight=\"" << e->GetWeight() << "\"];" << std::endl;
                }

                file << entry.first << "[label=\"" << entry.first << "\", name=\"" << entry.first << "\", ";
                if (this->GetVertexContent(entry.first)->IsAnd())
                    file << "AO=\"A\"];" << std::endl;
                else
                    file << "AO=\"O\"];" << std::endl;
            }

            std::vector<VertexId> root;
            for(VertexId v = 0;v < this->NumberOfVertices();v++)
                if(this->IsValidVertex(v))
                    if(this->ParentsVertices(v).size() == 0)
                        root.push_back(v);

            if(root.size()>1){
                VertexId rid = this->NumberOfVertices()+1;
                file << rid << "[label=\"" << rid << "\", name=\"" << rid << "\", AO=\"X\"];";
                for(auto r : root)
                    file << rid << " -> " << r << "[weight=\"1\"];" << std::endl;
            }

            file << "}";
        }
        else
            std::cout << "Unable to open " << dir << std::endl;

        return;
    }

    template<class StateType>    
    std::vector<Symbolic_State<StateType> > T_AOG<StateType>::GetTopLevelStates()
    {
        std::vector<Symbolic_State<StateType> > top_level_states;
        std::unordered_set<Symbolic_State<StateType> > sources;
        std::unordered_set<Symbolic_State<StateType> > results;
        for (Symbolic_Rule<StateType> rule : this->GetRules())
        {
            sources.insert(rule.GetSource());
            results.insert(rule.GetResults().begin(), rule.GetResults().end());
        }

        // check which source is not other sources' result
        for (Symbolic_State<StateType> source : sources)
        {
            if (results.find(source) == results.end())
                top_level_states.push_back(source);
        }
        return top_level_states;
    }

    template<class StateType>    
    std::vector<Symbolic_Rule<StateType> > T_AOG<StateType>::GetTopLevelRules()
    {
        std::vector<Symbolic_Rule<StateType> > top_level_rules;
        Symbolic_State<StateType> root_state = this->GetStateByVertexId(this->GetRoot());
        for (Symbolic_Rule<StateType> rule : this->GetRules())
            if (rule.GetSource() == root_state)
                top_level_rules.push_back(rule);

        return top_level_rules;
    }

    // template<class StateType>
    // void T_AOG<StateType>::CleanRules()
    // {
    //     // // create a map from source to all possible results (each results are vector of Symbolic_State)
    //     // std::unordered_map<Symbolic_State<StateType>, std::vector<std::vector<Symbolic_State<StateType> > > > rules_map;
    //     // for (auto it = this->all_rules_.begin(); it != this->all_rules_.end(); it++)
    //     // {   
    //     //     // if the source does not exist in map, create a vector
    //     //     if (rules_map.find() == rules_map.end())
    //     //         rules_map.insert(std::make_pair(it->GetSource(), std::vector<std::vector<Symbolic_State<StateType> > >(0,0)));
    //     //     rules_map[it->GetSource()].push_back(it->GetResults());
    //     // }

    //     auto it = this->GetRules().begin();
    //     while (it != this->GetRules().end())
    //     {
    //         std::vector<Symbolic_Rule<StateType> > rules = this->GetRules();
            
    //         // if the rule only has one children (i.e. it is a redundant rule)
    //         // connect its children to its parent
    //         // i.e. rule A->B->C
    //         if (it != this->GetRules().end() && it->GetResults().size() == 1)
    //         {
    //             bool has_parent = false;
    //             // the state to be connected to its parent's parent
    //             Symbolic_State<StateType> state_to_be_connected = it->GetResults()[0];
    //             Symbolic_State<StateType> state_to_be_deleted = it->GetSource();

    //             // find parent's parent
    //             for (Symbolic_Rule<StateType> rule : rules)
    //             {
    //                 std::vector<Symbolic_State<StateType> > result_vector = rule.GetResults();
    //                 auto result_it = std::find(result_vector.begin(), result_vector.end(), state_to_be_deleted);
    //                 // if this dummy 
    //                 if (result_it != result_vector.end())
    //                 {
    //                     // delete B->C
    //                     this->DeleteRule(*it);
    //                     // add A->C
    //                     std::vector<Symbolic_State<StateType> > result(result_vector.begin(), result_it);
    //                     result.push_back(state_to_be_connected);
    //                     result.insert(result.end(), result_it+1, result_vector.end());
    //                     this->AddRule(Symbolic_Rule<StateType>(rule.GetSource(), result));
    //                     // delete A->B
    //                     this->Delete(rule);
    //             }
    //             if (has_parent)Rule(rule);
    //                     has_parent = true;
    //                     // update iterator
    //                 }
    //                 it = this->GetRules().begin();
    //             else 
    //                 it++;
    //         }
    //         else if (it == this->GetRules().end())
    //             break;
    //         else
    //             it++;
    //     }
    // }

    template <class StateType>
    double T_AOG<StateType>:: GrammarSize(double leaf_bias, bool simplified_graph)
    {
        double grammar_size = 0;
        if(!simplified_graph){
            for(const auto& rule : this->GetRules())
            {
                ++grammar_size;
                for(const auto& state : rule.GetResults())
                {
                    if(state.GetId() == -1)
                        grammar_size += leaf_bias;
                    else    
                        ++grammar_size;
                }
            }
        }
        else{
            std::queue<VertexId> q;
            q.push(this->GetRoot());
            while(!q.empty()){
                VertexId curr = q.front();
                q.pop();
                std::vector<VertexId> children = this->ChildrenVertices(curr);
                if(this->GetVertexContent(curr)->IsAnd()){
                    grammar_size++;
                }
                else{
                    grammar_size += children.size();
                }
                for(VertexId child : children){
                    if(this->GetStateByVertexId(child).GetIsBasic())
                        grammar_size += leaf_bias;
                    else{
                        q.push(child);
                        grammar_size += 1;
                    }
                }
            }
        }

        return grammar_size;
    }


    template<class StateType>
    void T_AOG<StateType>::SimplifyGraph(){
        std::cout << "Simplifying..." << std::endl;

        //get max state id
        int stateId = 0;
        for(auto it : this->state_to_vertex_){
            if(it.first.GetId() > stateId){
                stateId = it.first.GetId();
            }
        }
        stateId++;
        
        auto delSameParAndChildState = [this](){
            std::queue<VertexId> q;
            q.push(this->GetRoot());

            bool OuterChanged = false;

            while(!q.empty()){
                VertexId curr = q.front();
                bool changed = false;
                unsigned childPos = 0;
                std::vector<VertexId> children = this->ChildrenVertices(curr);

                for(VertexId child : children){
                    
                    bool childIsAnd = this->GetVertexContent(child)->IsAnd();

                    if(!this->GetStateByVertexId(child).GetIsBasic() && 
                        childIsAnd == this->GetVertexContent(curr)->IsAnd() &&
                        this->ParentsVertices(child).size() == 1){
                        
                        OuterChanged = true;
                        changed = true;
                        this->DeleteEdge(curr, child);

                        if(childIsAnd){
                            for(VertexId grandChild : this->ChildrenVertices(child)){
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(curr, grandChild, tmpEdge, true, childPos);
                                childPos++;
                                if(!this->GetStateByVertexId(grandChild).GetIsBasic()){
                                    q.push(grandChild);
                                }
                            }
                            childPos--;
                        }
                        else {
                            for(VertexId grandChild : this->ChildrenVertices(child)){
                                auto it = std::find(children.begin(), children.end(), grandChild);
                                double weight = this->GetEdgeContent(child, grandChild)[0]->GetWeight();
                                if(it == children.end()){
                                    std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>(weight);
                                    this->AddEdge(curr, grandChild, tmpEdge, false); //or does not have multiedge
                                    if(!this->GetStateByVertexId(grandChild).GetIsBasic()){
                                        q.push(grandChild);
                                    }
                                }
                                else{ //rare case
                                    std::unordered_map<VertexId, double> m = this->GetOutEdgeWeights(curr, false);
                                    m[*it] += weight;
                                    this->SetOutEdgeWeights(curr, m);
                                }
                            }
                        }

                        this->DeleteVertex(child);
                    }

                    else if(!this->GetStateByVertexId(child).GetIsBasic()){
                        q.push(child);
                    }

                    childPos++;
                }
                
                if(!changed){
                    q.pop();
                }
            }

            return OuterChanged;
        };

        auto moveCommonChild = [this, &stateId](){

            std::queue<VertexId> q;
            q.push(this->GetRoot());

            bool OuterChanged = false;

            while(!q.empty()){
                VertexId curr = q.front();
                q.pop();

                const std::vector<VertexId>& children = this->ChildrenVertices(curr);

                
                if(!this->GetVertexContent(curr)->IsAnd() && children.size() > 1){
                    // //record the common grandchild at the first and last position of all children of all children of current OR node
                    // VertexId fstGrandChild = 0; //cannot be root node, so used for initial null state
                    // VertexId lstGrandChild = 0;

                    // //check common grandChild and its position in children
                    // //check children are all AND nodes
                    
                    // for(VertexId child : children){
                    //     if(!this->GetVertexContent(child)->IsAnd()){
                    //         fstGrandChild = 0;
                    //         lstGrandChild = 0;
                    //         break;
                    //     }

                    //     std::vector<VertexId> grandChildren = this->ChildrenVertices(child);
                    //     if(grandChildren.empty()){
                    //         fstGrandChild = 0;
                    //         lstGrandChild = 0;
                    //         break;
                    //     }

                    //     if(fstGrandChild == 0 && lstGrandChild == 0){
                    //         fstGrandChild = grandChildren.front();
                    //         lstGrandChild = grandChildren.back();
                    //     }
                    //     else{
                    //         bool fstGrandChildPassed = true;
                    //         if(fstGrandChild != grandChildren.front()){
                    //             fstGrandChildPassed = false;
                    //             fstGrandChild = 0;
                    //         }
                    //         if(lstGrandChild != grandChildren.back()){
                    //             lstGrandChild = 0;
                    //             if(!fstGrandChildPassed){
                    //                 break;
                    //             }
                    //         }
                    //     }
                    // }

                    // if(fstGrandChild == 0 && lstGrandChild == 0){
                    //     for (VertexId child : children){
                    //         if(!this->GetStateByVertexId(child).GetIsBasic()){
                    //             q.push(child);
                    //         }
                    //     }
                    //     continue;
                    // }

                    // if(fstGrandChild == lstGrandChild){
                    //     lstGrandChild = 0;
                    // }

                    //find map from grandchildren in 1st position and last position to its parents (at child level wrt. curr)
                    std::unordered_map<VertexId, std::vector<VertexId>> fstGrandChildren, sndGrandChildren;
                    for (VertexId child : children){
                        const std::vector<VertexId>& grandChildren = this->ChildrenVertices(child);
                        if(grandChildren.empty()){
                            continue;
                        }
                        VertexId fstGrandChild = grandChildren.front();
                        VertexId sndGrandChild = grandChildren.back();
                        auto fst_it = fstGrandChildren.find(fstGrandChild);
                        if(fst_it == fstGrandChildren.end()){
                            fstGrandChildren[fstGrandChild] = {child};
                        }
                        else{
                            fst_it->second.push_back(child);
                        }
                        auto snd_it = sndGrandChildren.find(sndGrandChild);
                        if(snd_it == sndGrandChildren.end()){
                            sndGrandChildren[sndGrandChild] = {child};
                        }
                        else{
                            snd_it->second.push_back(child);
                        }
                    }
                    
                    //find map from child level to their common grandchildren
                    std::map<std::vector<VertexId>, std::pair<VertexId, VertexId>> Children2Grand;
                    for(const auto& x : fstGrandChildren){
                        if(x.second.size() >= 2)
                            Children2Grand[x.second] = {x.first, 0};
                    }
                    for(const auto& x : sndGrandChildren){
                        if(x.second.size() >= 2){
                            auto it = Children2Grand.find(x.second);
                            if(it == Children2Grand.end()){
                                Children2Grand[x.second] = {0, x.first};
                            }
                            else if(it->second.first != x.first){
                                it->second.second = x.first;
                            }
                        }
                    }
                    
                    if(Children2Grand.empty()){
                        for (VertexId child : children){
                            if(!this->GetStateByVertexId(child).GetIsBasic()){
                                q.push(child);
                            }
                        }
                        continue;
                    }

                    //if conflict in child level occurs, keep the child level with larger size
                    std::vector<std::pair<std::vector<VertexId>, std::pair<VertexId, VertexId>>> ch2g; 
                    std::copy(Children2Grand.begin(), Children2Grand.end(), back_inserter(ch2g));
                    std::sort(ch2g.begin(), ch2g.end(), 
                            [](const std::pair<std::vector<VertexId>, std::pair<VertexId, VertexId>> &a, const std::pair<std::vector<VertexId>, std::pair<VertexId, VertexId>> & b)
                            {return a.first.size() < b.first.size();});
                    
                    std::unordered_set<VertexId> existVertex;
                    for(auto it = ch2g.begin(); it != ch2g.end();){
                        bool found = false;
                        for(VertexId x : it->first){
                            if(existVertex.find(x) != existVertex.end()){
                                it = ch2g.erase(it);
                                found = true;
                                break;
                            }
                        }
                        if(found){
                            continue;
                        }
                        for(VertexId x : it->first){
                            existVertex.insert(x);
                        }
                        ++it;
                    }

                    //push no common grandchildren child to the queue
                    for(VertexId child : children){
                        if(existVertex.find(child) == existVertex.end() && !this->GetStateByVertexId(child).GetIsBasic()){
                            q.push(child);
                        }
                    }
                    
                    for(const auto& x : ch2g){
                        //test for children, should have no parents other than curr
                        const std::vector<VertexId>& andVertices = x.first;
                        VertexId fstGrandChild = x.second.first;
                        VertexId sndGrandChild = x.second.second;

                        bool noOtherParent = true;
                        for(VertexId child : andVertices){
                            std::vector<VertexId> otherParents = this->ParentsVertices(child);
                            if(otherParents.size() != 1){
                                noOtherParent = false;
                                break;
                            }
                        }

                        if(!noOtherParent){
                            for (VertexId child : andVertices){
                                if(!this->GetStateByVertexId(child).GetIsBasic()){
                                    q.push(child);
                                }
                            }
                            continue;
                        }

                        std::cout << "[simplify] pass all tests" << std::endl;

                        //separate cases for:
                        //1. all and nodes share common children, or
                        //2. some of the and nodes do
                        //the separation of cases is for the purpose of keep state id unchanged after truncate() for case 1
                        if(andVertices.size() == children.size()){ //case 1
                            std::shared_ptr<Symbolic_State<StateType> > dummy_state = std::make_shared<Symbolic_State<StateType>>(stateId++);
                            std::shared_ptr<AOG_Vertex<StateType> > dummy = std::make_shared<AOG_Vertex<StateType>>(*dummy_state, true, false);
                            unsigned dummy_id = this->AddVertex(dummy);

                            //connect dummy to parents
                            for(VertexId parent : this->ParentsVertices(curr)){
                                if(!this->GetVertexContent(parent)->IsAnd()){ //need to consider weight, no multiedge case, and no position requirement
                                    std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>(this->GetEdgeContent(parent, curr)[0]->GetWeight());
                                    this->AddEdge(parent, dummy_id, tmpEdge, false);
                                }
                                else{
                                    std::vector<VertexId> childrenOfPar = this->ChildrenVertices(parent);
                                    std::vector<unsigned> pos;
                                    auto iter = childrenOfPar.begin();
                                    auto iter_begin = childrenOfPar.begin();
                                    while ((iter = std::find(iter, childrenOfPar.end(), curr)) != childrenOfPar.end()){
                                        pos.push_back(iter - iter_begin);
                                        iter++;
                                    }
                                    this->DeleteEdge(parent, curr);
                                    for(unsigned innerPos : pos){
                                        std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                        this->AddEdge(parent, dummy_id, tmpEdge, true, innerPos);
                                    }
                                }
                            }

                            if(curr == this->GetRoot()){
                                this->root_id_ = dummy_id;
                                this->GetVertexContent(curr)->SetIsRoot(false);
                                this->GetVertexContent(dummy_id)->SetIsRoot(true);
                            }

                            //connect dummy to fstGrandChild (if shared), curr, and lstGrandChild (if shared)
                            if(fstGrandChild != 0){
                                for(VertexId child : children){
                                    this->DeleteEdge(child, fstGrandChild);
                                }
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(dummy_id, fstGrandChild, tmpEdge, true);
                                if(!this->GetStateByVertexId(fstGrandChild).GetIsBasic()){
                                    q.push(fstGrandChild);
                                }
                            }

                            std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                            this->AddEdge(dummy_id, curr, tmpEdge, true);
                            q.push(curr);

                            if(sndGrandChild != 0){
                                for(VertexId child : children){
                                    this->DeleteEdge(child, sndGrandChild);
                                }
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(dummy_id, sndGrandChild, tmpEdge, true);
                                if(!this->GetStateByVertexId(sndGrandChild).GetIsBasic()){
                                    q.push(sndGrandChild);
                                }
                            }

                            OuterChanged = true;
                        }
                        
                        else{ //case 2
                            std::shared_ptr<Symbolic_State<StateType> > dummy_and_state = std::make_shared<Symbolic_State<StateType>>(stateId++);
                            std::shared_ptr<AOG_Vertex<StateType> > dummy_and = std::make_shared<AOG_Vertex<StateType>>(*dummy_and_state, true, false);
                            unsigned dummy_and_id = this->AddVertex(dummy_and);

                            std::shared_ptr<Symbolic_State<StateType> > dummy_or_state = std::make_shared<Symbolic_State<StateType>>(stateId++);
                            std::shared_ptr<AOG_Vertex<StateType> > dummy_or = std::make_shared<AOG_Vertex<StateType>>(*dummy_or_state, false, false);
                            unsigned dummy_or_id = this->AddVertex(dummy_or);
                            
                            // if(this->GetVertexContent(curr)->IsRoot()){
                            //     this->GetVertexContent(this->root_id_)->SetIsRoot(false);

                            //     std::shared_ptr<Symbolic_State<StateType> > dummy_root_state = std::make_shared<Symbolic_State<StateType>>(stateId++);
                            //     std::shared_ptr<AOG_Vertex<StateType> > dummy_root = std::make_shared<AOG_Vertex<StateType>>(*dummy_root_state, false, false);
                            //     unsigned dummy_root_id = this->AddVertex(dummy_root);
                            //     this->GetVertexContent(dummy_root_id)->SetIsRoot(true);
                            //     this->root_id_ = dummy_root_id;

                            //     unsigned weight = 0;
                            //     for(VertexId tmpChild : this->ChildrenVertices(curr)){
                            //         weight += this->GetEdgeContent(curr, tmpChild)[0]->GetWeight();
                            //     }
                            //     std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>(weight);
                            //     this->AddEdge(dummy_root_id, dummy_id, tmpEdge, false);
                            // }

                            //connect dummy_and to curr, disconnect selected children to curr, connect children to dummy_or
                            double dummy_and_weight = 0;
                            for(VertexId andVtx : andVertices){
                                double weight = this->GetEdgeContent(curr, andVtx)[0]->GetWeight();
                                this->DeleteEdge(curr, andVtx);
                                std::shared_ptr<AOG_Edge> dummy_or_edge = std::make_shared<AOG_Edge>(weight);
                                this->AddEdge(dummy_or_id, andVtx, dummy_or_edge, false);
                                dummy_and_weight += weight;
                            }
                            std::shared_ptr<AOG_Edge> dummy_and_edge = std::make_shared<AOG_Edge>(dummy_and_weight);
                            this->AddEdge(curr, dummy_and_id, dummy_and_edge, false);

                            // if(curr == this->GetRoot()){
                            //     this->root_id_ = dummy_id;
                            //     this->GetVertexContent(curr)->SetIsRoot(false);
                            //     this->GetVertexContent(dummy_id)->SetIsRoot(true);
                            // }

                            //connect dummy to fstGrandChild (if shared), curr, and lstGrandChild (if shared)
                            if(fstGrandChild != 0){
                                for(VertexId andVtx : andVertices){
                                    this->DeleteEdge(andVtx, fstGrandChild);
                                }
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(dummy_and_id, fstGrandChild, tmpEdge, true);
                                if(!this->GetStateByVertexId(fstGrandChild).GetIsBasic()){
                                    q.push(fstGrandChild);
                                }
                            }

                            std::shared_ptr<AOG_Edge> dummy_and_to_or_edge = std::make_shared<AOG_Edge>();
                            this->AddEdge(dummy_and_id, dummy_or_id, dummy_and_to_or_edge, true);
                            q.push(dummy_or_id);

                            if(sndGrandChild != 0){
                                for(VertexId andVtx : andVertices){
                                    this->DeleteEdge(andVtx, sndGrandChild);
                                }
                                std::shared_ptr<AOG_Edge> tmpEdge = std::make_shared<AOG_Edge>();
                                this->AddEdge(dummy_and_id, sndGrandChild, tmpEdge, true);
                                if(!this->GetStateByVertexId(sndGrandChild).GetIsBasic()){
                                    q.push(sndGrandChild);
                                }
                            }

                            OuterChanged = true;
                        }
                    }
                }
                else{
                    for (VertexId child : children){
                        if(!this->GetStateByVertexId(child).GetIsBasic()){
                            q.push(child);
                        }
                    }
                }
            }

            return OuterChanged;
        };

        while(true){
            bool changed = delSameParAndChildState();
            bool changed2 = moveCommonChild();
            changed = changed || changed2;
            if(!changed){
                break;
            }
            this->TruncateGraph();
        }
    }

    template<class StateType>
    void T_AOG<StateType>::TruncateGraph(){
        VertexId root_id = this->GetRoot();
        std::queue<VertexId> traverse_q;
        traverse_q.push(root_id);
        int count = 1;
        while (!traverse_q.empty())
        {
            VertexId cur_vertex_id = traverse_q.front();
            traverse_q.pop();

            std::vector<VertexId> children_id = this->ChildrenVertices(cur_vertex_id);
            for (VertexId child_id : children_id)
            {
                traverse_q.push(child_id);
            }

            if (children_id.size() == 1 && cur_vertex_id != this->GetRoot())
            {
                // find its parent, and connect child to parent
                std::vector<VertexId> parents_id = this->ParentsVertices(cur_vertex_id);

                // the reconnection is needed for all parents of cur_vertex_id
                for (VertexId p_id : parents_id)
                {
                    // delete old edge parent to cur
                    std::vector<VertexId> p_children = this->ChildrenVertices(p_id);
                    int pos = -1;
                    std::vector<std::vector<std::shared_ptr<AOG_Edge> > > parent_to_p;
                    for (int i = 0; i < p_children.size(); i++)
                    {
                        if (p_children[i] == cur_vertex_id)
                            pos = i;
                        parent_to_p.push_back(this->GetEdgeContent(p_id, p_children[i]));
                        this->DeleteEdge(p_id, p_children[i]);
                    }

                    for (int i = 0; i < pos; i++)
                    {
                        std::shared_ptr<AOG_Edge> new_edge = std::make_shared<AOG_Edge>();                        
                        // Or node will not have multi-edge
                        if (!this->GetVertexContent(p_id)->IsAnd())
                        {
                            if (parent_to_p[i].size() == 0)
                                throw std::exception();
                            if (parent_to_p[i].size() > 1)
                            {
                                std::cerr << "Multiedge\n";
                                throw std::exception();
                            }
                            new_edge->SetWeight(parent_to_p[i][0]->GetWeight());
                        }
                        // add new edge with weight if parent is Or
                        this->AddEdge(p_id, p_children[i], new_edge);
                    }

                    // repeat same operation for cur_vertex_id
                    // if parent is Or-node, copy its weight when creating the edge
                    // Or node will not have multi-edge
                    std::shared_ptr<AOG_Edge> new_edge = std::make_shared<AOG_Edge>();                        
                    if (!this->GetVertexContent(p_id)->IsAnd())
                    {
                        if (parent_to_p[pos].size() == 0)
                            throw std::exception();
                        if (parent_to_p[pos].size() > 1)
                        {
                            std::cerr << "Multiedge\n";
                            throw std::exception();
                        }
                        new_edge->SetWeight(parent_to_p[pos][0]->GetWeight());
                    }
                    // add new edge with weight if parent is Or
                    this->AddEdge(p_id, children_id[0], new_edge);

                    for (int i = pos+1; i < p_children.size(); i++)
                    {
                        std::shared_ptr<AOG_Edge> new_edge = std::make_shared<AOG_Edge>();                                                
                        // Or node will not have multi-edge
                        if (!this->GetVertexContent(p_id)->IsAnd())
                        {
                            if (parent_to_p[i].size() == 0)
                                throw std::exception();
                            if (parent_to_p[i].size() > 1)
                            {
                                std::cerr << "Multiedge\n";
                                throw std::exception();
                            }
                            new_edge->SetWeight(parent_to_p[i][0]->GetWeight());
                        }
                        // add new edge with weight if parent is Or
                        this->AddEdge(p_id, p_children[i], new_edge);
                    }
                }
                // delete old edge cur to child
                this->DeleteEdge(cur_vertex_id, children_id[0]);
                // delete dummy node

                this->DeleteVertex(cur_vertex_id);
            }
        }
    }

    template<class StateType>
    void T_AOG<StateType> :: DeleteNoParentRules(const Symbolic_State<StateType>& src_state)
    {
        VertexId vid = this -> GetVertexIdByState(src_state);
        
        
        //if the passed in vertex has parent, or the state is a leaf, directly return
        if(this->ParentsVertices(vid).size() || src_state.GetIsBasic() )
            return;
        
        // if(src_state.GetIsBasic())
        // {
        //     std::cerr<<"about to delete a leaf node in AOG!\n";
        //     throw std::exception();
        // }

        std::vector<VertexId> children_vtx = this->ChildrenVertices(vid);
        std::cerr<<"all children id: \n";
        for(auto child_vtx : children_vtx)
            std::cerr<<child_vtx<<" ";
        std::cerr<<std::endl;
        // if it is an and-node
        if(this->GetVertexContent(vid)->IsAnd())
        {
            //delete curent rule
            std::vector<Symbolic_State<StateType> > children_states;
            for(auto vertex : children_vtx )
                children_states.push_back(this->GetStateByVertexId(vertex));
            
            if(!this->DeleteRule(Symbolic_Rule<StateType>(src_state ,children_states)))
            {
                std::cerr<<"fail to delete rules in and-node:\n";
                // std:: cerr<<"("<<src_state.GetContent()<<", "<<src_state.GetId()<<")->\n";
                // for(auto child : children_states)
                //     std::cerr<<"("<<child.GetContent()<<", "<<child.GetId()<<")";
                // std::cerr<<std::endl;
                // int a;
                // std::cin >>a;
            }
            
            //recursively delete child states if they have no parents
            for(auto child_state : children_states)
                this->DeleteNoParentRules(child_state);
        }

        //if it is an or-node
        else
        {
            std::vector< std::vector< Symbolic_State<StateType> > > grand_children_states;

            //first create all rules
            for(auto child : children_vtx)            
            {

                std::vector<VertexId> grand_children_vtx = this->ChildrenVertices(child);
                std::vector<Symbolic_State<StateType> > seq;
                for(auto grand_child : grand_children_vtx )
                    seq.push_back(this->GetStateByVertexId(grand_child));
                grand_children_states.push_back(seq);
            }
                // std::cerr<< "current child vertex id: "<<child<<std::endl;

                
                // std::cerr<<"the child's grand children id: \n";
                // for(auto child_vtx : grand_children_vtx)
                //     std::cerr<<child_vtx<<" ";
                // std::cerr<<std::endl;

                // auto parents = this->ParentsVertices(406);
                // auto children = this->ChildrenVertices(651);
                // std::cerr<<"406's parent vertices: \n";
                // for(auto parent : parents)
                //     std::cerr << parent<<" ";
                // std::cerr <<std::endl;

                // std::cerr<<"651's children vertices: \n";                
                // for(auto child : children)
                //     std::cerr<<child<<" ";
                // std::cerr<<std::endl;

              
                //delete all rules
                for(int i = 0; i < grand_children_states.size(); ++i)
                {
                    if(!this->DeleteRule(Symbolic_Rule<StateType>(src_state ,grand_children_states[i])))
                    {
                        std::cerr<<"fail to delete rules in or-node\n";
                        // std:: cerr<<"("<<src_state.GetContent()<<", "<<src_state.GetId()<<")->\n";
                        // for(auto grand_child : grand_children_states)
                        //     std::cerr<<"("<<grand_child.GetContent()<<", "<<grand_child.GetId()<<")";
                        // std::cerr<<std::endl;

                        // for(auto rule : this->GetRules())
                        // {
                        //     std::cerr<<"("<<rule.GetSource().GetContent()<<", "<<rule.GetSource().GetId()<<")->\n";
                        //     for(auto state : rule.GetResults())
                        //     {
                        //         std::cerr<<"("<<state.GetContent()<<", "<<state.GetId()<<")";
                        //     }
                        //     std::cerr<<std::endl;
                        // }
                        // int a;
                        // std::cin >>a;
                    }

                    for(auto state : grand_children_states[i])
                        this->DeleteNoParentRules(state);
                }
               
                
            
        }
    }



    template<class StateType>
    void T_AOG<StateType>::OutputGraph(std::string filename, std::string dir, bool truncate, bool simplify)
    {
    
        std::cerr << "PRINTED GRAPH\n";
        std::shared_ptr<T_AOG<StateType> > printed_graph = std::make_shared<T_AOG<StateType> >(*this);

        //truncate flag -----------------------------
        if (truncate)
        {
            printed_graph->TruncateGraph();
        }

        std::string dir_copy = dir;
        std::string filename_copy = filename;
        std::ofstream file;
        if (dir_copy.empty())
            dir_copy = filename_copy.append(".txt");
        else
            dir_copy = dir_copy.append("/").append(filename_copy.append(".txt"));

        file.open(dir_copy, std::ofstream::out|std::ofstream::trunc);
        if (file.is_open()) {

            VertexId root_id = printed_graph->GetRoot();
            std::queue<VertexId> q;
            q.push(root_id);
            std::queue<unsigned> level_q;
            level_q.push(1);

            std::unordered_set<VertexId> plot_visited;

            while (!q.empty())
            {
                VertexId cur_vertex = q.front();
                unsigned level = level_q.front();
                q.pop();
                level_q.pop();

                std::vector<VertexId> children_vertices_id = printed_graph->ChildrenVertices(cur_vertex);

                if(!printed_graph->GetVertexContent(cur_vertex))
                {
                    int a;
                }
                bool isAnd = printed_graph->GetVertexContent(cur_vertex)->IsAnd();
                std::unordered_map<VertexId, double> weights;
                Symbolic_State<StateType> cur_vertex_state = printed_graph->GetStateByVertexId(cur_vertex);
                
                if (!isAnd)
                    weights = printed_graph->GetOutEdgeWeights(cur_vertex, true);
                    
                for (int i = 0; i < children_vertices_id.size(); i++)
                {
                    Symbolic_State<StateType> child_state = printed_graph->GetStateByVertexId(children_vertices_id[i]);
                    if (isAnd)
                    {    
                        file << cur_vertex << "," << children_vertices_id[i] << "," << cur_vertex_state.GetContent() << "_" << cur_vertex_state.GetId() 
                            << "," << child_state.GetContent() << "_" << child_state.GetId() << ",," << i+1 << "," << level << "\n";
                    }
                    else
                    {
                        file << cur_vertex << "," << children_vertices_id[i] << "," << cur_vertex_state.GetContent() << "_" << cur_vertex_state.GetId() 
                            << "," << child_state.GetContent() << "_" << child_state.GetId() << "," << weights[children_vertices_id[i]] << "," << i+1 <<"," << level << "\n";
                    }
                    q.push(children_vertices_id[i]);
                    level_q.push(level+1);
                }
            }
            file.close();
        }
        else
            std::cout << "Unable to open " << dir_copy << std::endl;



        //simplify flag -----------------------------
        if (simplify)
        {
            printed_graph->SimplifyGraph();
        }
        else
        {
            return;
        }

        std::string dir_copy2 = dir;
        std::string filename_copy2 = filename;
        std::ofstream file_simplify;
        if (dir_copy2.empty())
            dir_copy2 = filename_copy2.append("_simplify.txt");
        else
            dir_copy2 = dir_copy2.append("/").append(filename_copy2.append("_simplify.txt"));

        file_simplify.open(dir_copy2, std::ofstream::out|std::ofstream::trunc);
        if (file_simplify.is_open()) {

            VertexId root_id = printed_graph->GetRoot();
            std::queue<VertexId> q;
            q.push(root_id);
            std::queue<unsigned> level_q;
            level_q.push(1);

            std::unordered_set<VertexId> plot_visited;

            while (!q.empty())
            {
                VertexId cur_vertex = q.front();
                unsigned level = level_q.front();
                q.pop();
                level_q.pop();

                std::vector<VertexId> children_vertices_id = printed_graph->ChildrenVertices(cur_vertex);

                if(!printed_graph->GetVertexContent(cur_vertex))
                {
                    int a;
                }
                bool isAnd = printed_graph->GetVertexContent(cur_vertex)->IsAnd();
                std::unordered_map<VertexId, double> weights;
                Symbolic_State<StateType> cur_vertex_state = printed_graph->GetStateByVertexId(cur_vertex);
                
                if (!isAnd)
                    weights = printed_graph->GetOutEdgeWeights(cur_vertex, true);
                    
                for (int i = 0; i < children_vertices_id.size(); i++)
                {
                    Symbolic_State<StateType> child_state = printed_graph->GetStateByVertexId(children_vertices_id[i]);
                    if (isAnd)
                    {    
                        file_simplify << cur_vertex << "," << children_vertices_id[i] << "," << cur_vertex_state.GetContent() << "_" << cur_vertex_state.GetId() 
                            << "," << child_state.GetContent() << "_" << child_state.GetId() << ",," << i+1 << "," << level << "\n";
                    }
                    else
                    {
                        file_simplify << cur_vertex << "," << children_vertices_id[i] << "," << cur_vertex_state.GetContent() << "_" << cur_vertex_state.GetId() 
                            << "," << child_state.GetContent() << "_" << child_state.GetId() << "," << weights[children_vertices_id[i]] << "," << i+1 <<"," << level << "\n";
                    }
                    q.push(children_vertices_id[i]);
                    level_q.push(level+1);
                }
            }
            file_simplify.close();
        }
        else
            std::cout << "Unable to open " << dir_copy2 << std::endl;

        return;
    }

    template<class StateType>
    void T_AOG<StateType>::PrintOrNodeAndChildrenCount(std::ofstream &file)
    {
        int or_node_count = 0;
        int or_children_count = 0;
        VertexId root_id = this->GetRoot();
        std::queue<VertexId> q;
        q.push(root_id);
        while(!q.empty())
        {
            VertexId cur_ver = q.front();
            q.pop();
            if(!this->GetVertexContent(cur_ver)->IsAnd())
            {
                ++or_node_count;
                or_children_count += this->ChildrenVertices(cur_ver).size();
            }
            auto children_vtx = this->ChildrenVertices(cur_ver);
            for(auto vtx : children_vtx)
                q.push(vtx);
        }
        file << "total or nodes: "<<or_node_count<<" total_or_children: "<< or_children_count<<std::endl;
    }

    template<class StateType>
    void T_AOG<StateType>::OutputLearnedTree(std::string path, std::string appendix)
    {
        std::vector<Symbolic_Rule<StateType>> rules = this->GetRules();
        std::string file_name;

        file_name = appendix == "" ? path + "learned_tree.txt" : path + "learned_tree_" + appendix + ".txt";
        
        std::ofstream file;
        file.open(file_name, std::ofstream::out|std::ofstream::trunc);

        if (file.is_open())
        {
            std::cerr << "File opened. \n";
            for (const Symbolic_Rule<std::string>& rule : rules)
            {
                // size of result, souce id, source content, [child id, child content]*
                file << rule.GetResults().size();
                std::vector<Symbolic_State<std::string> > states = rule.GetResults();
                VertexId source_vtx_id = this->GetVertexIdByState(rule.GetSource());
                bool isAnd = this->GetVertexContent(source_vtx_id)->IsAnd();
                // std::cerr << "state_id: " << rule.GetSource().GetId() << ", vertex_id: " << source_vtx_id << ", isAnd: " << isAnd << std::endl;
                if(!isAnd){
                    std::vector<VertexId> children_vector;
                    for(const Symbolic_State<std::string>& state : states){
                        children_vector.push_back(this->GetVertexIdByState(state));
                    }
                    VertexId tmpDummy;
                    bool found = false;
                    for (VertexId dummy : this ->ChildrenVertices(source_vtx_id))
                    {
                        if (this->ChildrenVertices(dummy) == children_vector)
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
                    double weight = this->GetOutEdgeWeights(source_vtx_id, false)[tmpDummy];
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

            file.close();
        }
        else
            std::cout << "Unable to open " << "learned_tree.txt" << std::endl;


    }

}

#endif //AOG_LIB_T_AOG_H
