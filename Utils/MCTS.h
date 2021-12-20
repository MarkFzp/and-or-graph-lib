//
// Created by Luyao Yuan on 18/2/21.
//

#ifndef AOG_LIB_MCTS_H
#define AOG_LIB_MCTS_H

#include <memory>
#include <vector>
#include "../Learner/Online_Learner.h"
#include <chrono>
#define HAS_RESOURCES 1
extern double Cp;

namespace AOG_LIB
{
    class Search_Node; // forward declaration
    typedef std::shared_ptr<Search_Node> SN_Ptr;

    class Search_Node
    {
        SN_Ptr parent_;
        std::vector<SN_Ptr> children_;
        unsigned visit_count_;
        double accumulated_value_;
		int level_;

    protected:
	    //used in constructors of derived classes
	    void Clean()
	    {
			this->level_ = 0;
		    this->parent_ = 0;
		    this->children_.resize(0);
		    this->visit_count_ = 0;
		    this->accumulated_value_ = 0;
	    }

    public:

		// Default constructor
		Search_Node();

		// Constructor to be called by Checkpoint class
		Search_Node(SN_Ptr, std::vector<SN_Ptr>, unsigned, double, int);
		//TODO basic usage functions
		virtual void SetParent(SN_Ptr);
		virtual void AddChildren(const SN_Ptr&);
		//void ClearChildren(){this->children_.resize(0);};
		SN_Ptr GetParent() const;
		std::vector<SN_Ptr> GetChildren() const; 
		unsigned GetVisitCount() const;
		double GetValue() const;
		int GetLevel() const;
		void SetLevel(int);

		//Override this function to expand
        virtual SN_Ptr Expand() = 0;

	    //Override this function to roll out
	    virtual SN_Ptr Simulate() = 0;

	    //Override this function to select children
	    //The argument is the parameter used for UCT
	    virtual SN_Ptr Select(double,bool) = 0;

		virtual void Backtrack(SN_Ptr) = 0;

	    //update mean reward as value of this search node
        void Update(double reward)
        {
            this->accumulated_value_ += (reward - this->accumulated_value_) / (++this->visit_count_);
        }
    };

	Search_Node::Search_Node()
	{ this->Clean(); }

	Search_Node::Search_Node(SN_Ptr parent, std::vector<SN_Ptr> children, 
		unsigned visit_count, double accumulated_value, int level)
	{  
		this->parent_ = parent;
		this->children_ = children;
		this->visit_count_ = visit_count;
		this->accumulated_value_ = accumulated_value;
		this->level_ = level;
	}

	void Search_Node::SetParent(SN_Ptr parent)
	{ this->parent_ = parent; }

	void Search_Node::AddChildren(const SN_Ptr& child)
	{ this->children_.push_back(child); }

	SN_Ptr Search_Node::GetParent() const
	{ return this->parent_; }

	std::vector<SN_Ptr> Search_Node::GetChildren() const
	{ return this->children_; }

	unsigned Search_Node::GetVisitCount() const 
	{ return this->visit_count_; }

	double Search_Node::GetValue() const 
	{ return this->accumulated_value_; }

	int Search_Node::GetLevel() const
	{ return this->level_; }

	void Search_Node::SetLevel(int level)
	{ this->level_ = level; }

    class MCTS
    {
        SN_Ptr root_;
		unsigned resource_;
    public:
	    //TODO MCTS core function
	    /*
		 * This implementation of MCTS is different from vanilla MCTS
		 * Since in every cycle of search, we want to choose among the best leaf search nodes,
		 * which are nodes unexpandable. Whereas vanilla MCTS is used to decide the next state
		 * to go after current search root. Hence, we need to return all leaf search nodes
		 * so that we can choose the best among them.
	     *
	     * Potential implementation:
	     *      call iteratively select() until current end node
	     *      call expand() to generate a new children
	     *      call simulate() [P.S. no new search nodes should be created,
	     *                       since we only care about the final reward
	     *                       and final result. Consider this as a tail recursion flavor thing]
	     *      simulate() should take care of the update part
	     *      memorize the final result return from simulate()
	     *      return all final results
		*/
	    std::vector<SN_Ptr> Search();

	    //TODO MCTS update root with the input of the function
		SN_Ptr GetRoot() { return this->root_; }
	    void UpdateRoot(SN_Ptr root){ this->root_ = root; }
		void DeleteRoot() { this->root_.reset(); }
		void UpdateResource(const unsigned resource) { this-> resource_ = resource; }
    };

	std::vector<SN_Ptr> MCTS::Search() 
	{
		/** std::cerr << "MCTS Search:" << std::endl;	 */	
		// vector to memorize the result of rewards
		std::vector<SN_Ptr> search_res(0, SN_Ptr());
		// variable used for UCT

		// while still have resources to run mcts
		while (this->resource_--){
			/** std::cerr<<"------------------------Start Next Round---------------------------\n"; 
			// find the most urgent expandable node
			std::cerr<<"------------------------Before Select---------------------------\n"; */
			std::cerr << "RESOURCE LEFT: " << this->resource_ << "\n";
			
			auto begin = std::chrono::system_clock::now();

			SN_Ptr node_to_expand = this->root_->Select(Cp,false);
			/** std::cerr << "MCTS Search: resource left:" << this->resource_ << std::endl; */

			// expand the node, and now node_to_expand becomes the newly expanded node
			/** std::cerr << "MCTS Search: before expand" << std::endl;		 */
			SN_Ptr expanded_node = node_to_expand->Expand();
			
			SN_Ptr termination_node;
			bool fail_to_expand = false;
			while (expanded_node == node_to_expand)
			{
				// if its children size is 0, this is the end state
				if (expanded_node->GetChildren().size() == 0)
				{
					std::cout << "[Cannot Expand, Already terminal]\n";
					termination_node = expanded_node;
					fail_to_expand = true;
					break;
				}
				else
				{
					std::cerr << "[Reselecting]\n";
					node_to_expand = expanded_node->Select(Cp,true);
					expanded_node = node_to_expand->Expand();
				}
			}

			// if the node expanded successfully
			if (!fail_to_expand)
			{
				termination_node = expanded_node->Simulate();
			}
			// if the node cannot expand
			// expanded_node itself is termination node
			expanded_node->Backtrack(termination_node);
			search_res.push_back(termination_node);
			auto end = std::chrono::system_clock::now();

			std::chrono::duration<double> learning_time = end-begin;
			std::cout << "The learning time for resource "<<this->resource_<< ": "<<learning_time.count()<<std::endl;

		}


		return search_res;
	}

}

#endif //AOG_LIB_MCTS_H
