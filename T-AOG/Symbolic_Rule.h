//
// Created by yuanluyao on 10/23/17.
//

#ifndef AOG_LIB_SYMBOLIC_RULE_H
#define AOG_LIB_SYMBOLIC_RULE_H

#include <utility>
#include "Symbolic_State.h"
#include<boost/functional/hash.hpp>
namespace AOG_LIB
{
    /* This class defines symbolic rules used to construct an AOG Graph
     * rules are all in chomsky normal form; each rule looks like
     * source_ -> result1_ result2_
     */
    template<class StateType>
    class Symbolic_Rule
    {
        Symbolic_State<StateType> source_;
        /* Symbolic_State<StateType> result1_; */
        /* Symbolic_State<StateType> result2_; */
	    std::vector<Symbolic_State<StateType> > results_;

    public:
        // constructs a rule with given source and resulting symbolic states
        Symbolic_Rule(const Symbolic_State<StateType> &,
                      const std::vector<Symbolic_State<StateType> >&);
        //return the source symbolic state of the rule
        const Symbolic_State<StateType>& GetSource() const;
        //return the rule's results
        const std::vector<Symbolic_State<StateType> >& GetResults() const;    
      

        /* This function compares two symbolic rules
        * @params: rhs: the other symbolic rule to be compared with
        * @return: true if the source and resulting symbolic states of the two rules are all the same
        */
        bool operator==(const Symbolic_Rule<StateType> &) const;

        // friend std::ostream& operator <<(std::ostream& os, const Symbolic_Rule<StateType> &);
    };
}

namespace std
{
    template<class StateType>
    struct hash<AOG_LIB::Symbolic_Rule<StateType> >
    {
        size_t operator()(const AOG_LIB::Symbolic_Rule<StateType> &rule) const noexcept
        {
            
            using std::string;
	        size_t vec_seed = boost::hash_range(rule.GetResults().begin(), rule.GetResults().end());
	    
            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:

            return (hash<AOG_LIB::Symbolic_State<StateType> >{}(rule.GetSource())
                     ^ vec_seed);
        }
    };


}


template <class StateType>
AOG_LIB::Symbolic_Rule<StateType>::Symbolic_Rule(const Symbolic_State<StateType> & source,
						 const std::vector<Symbolic_State<StateType> > & results)
        :source_(source), results_(results){}


template <class StateType>
const AOG_LIB::Symbolic_State<StateType>& AOG_LIB::Symbolic_Rule<StateType>::GetSource() const
{return this->source_;}

template <class StateType>
const std::vector<AOG_LIB::Symbolic_State<StateType> >&
AOG_LIB::Symbolic_Rule<StateType>::GetResults() const
{return this->results_;}

template <class StateType>
bool AOG_LIB::Symbolic_Rule<StateType>::operator== (const Symbolic_Rule<StateType>& rhs) const
{
  //compare # of results and the source
  if(this->results_.size() != rhs.results_.size() || this->source_ != rhs.source_ )
      return false;

  size_t length = this->results_.size();
  for(int i = 0; i < length; i++)
  {
      if(this->results_[i] != rhs.results_[i])
          return false;
  }
  return true;
}

// template <class StateType>
// std::ostream& operator <<(std::ostream& os, const AOG_LIB::Symbolic_Rule<StateType>& rule) 
// {
//     // os << "Rule is:\n";
//     // AOG_LIB::Symbolic_State<StateType> source = rule.GetSource();
//     // std::vector<AOG_LIB::Symbolic_State<StateType> > results = rule.GetResults();
//     // for (AOG_LIB::Symbolic_State<StateType> result : results)
//     // {
//     //     os << "\t" << "(" << source.GetContent() << "," << source.GetId() << ")" 
//     //         << " -> " << "(" << result.GetContent() << "," << result.GetId() << ")" << "\n";
//     // }
//     return os;
// }

#endif //AOG_LIB_SYMBOLIC_RULE_H

