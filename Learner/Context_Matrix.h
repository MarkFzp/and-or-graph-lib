//
// Created by Luyao Yuan on 18/2/22.
//

#ifndef AOG_LIB_CONTEXT_MATRIX_H
#define AOG_LIB_CONTEXT_MATRIX_H

#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/member.hpp>

#include "Learner.h"

template <class T>
using SequenceType = std::vector<AOG_LIB::Symbolic_State<T> >;

template<class T>
bool operator==(const SequenceType<T>& rhs, const SequenceType<T>& lhs)
{
    if(rhs.size() != lhs.size()) return false;
    unsigned i = 0;
    while(i < rhs.size())
        if(!(rhs[i] == lhs[i++]))
            return false;
    return true;
}

template <class T>
using ContextType = std::vector<SequenceType<T> >;

template <class T>
bool operator==(const ContextType<T>& rhs, const ContextType<T>& lhs)
{
    if(rhs.size() != lhs.size())
        return false;
    
    if(!rhs.size() || !lhs.size() )
    {
        std::cerr<<"one of the contexts is empty!\n";
        throw std::exception();
    }

    unsigned i = 0;
    while(i < rhs.size())
        if(!(rhs[i] == lhs[i++]))
            return false;
    return true;

}

template <class StateType>
size_t hash_value(const SequenceType<StateType> &seq)
{
    // using std::string;
    size_t vec_seed = boost::hash_range(seq.begin(),seq.end());
    return vec_seed;
}

namespace std
{
    template<class StateType>
    struct hash<SequenceType<StateType> >
    {
        size_t operator()(const SequenceType<StateType> &seq) const noexcept
        {

            // using std::string;
            boost::hash<SequenceType<StateType> > seq_hasher;
            return seq_hasher(seq);
        }
    };

    template<class StateType>
    struct hash<ContextType<StateType> >
    {
        size_t operator()(const ContextType<StateType> &context) const noexcept
        {

            using std::string;
            size_t vec_seed = boost::hash_range(context.begin(), context.end());
            // size_t vec_seed_2 = boost::hash_range(context.second.begin(), context.second.end());

            // Compute individual hash values for first,
            // second and third and combine them using XOR
            // and bit shifting:

            return vec_seed;
        }
    };
}

namespace AOG_LIB
{
    //An entry in the context matrix, the data can be recovered by
    //[context_.first, configuration, context_.second]
    template <class StateType>
    struct CM_Entry
    {
        ContextType<StateType> context_;
        std::vector<SequenceType<StateType> >  configuration_;
        unsigned count_;
    };
}

template <class T>
using Context_Matrix = boost::multi_index::multi_index_container
<
    AOG_LIB::CM_Entry<T>,
    boost::multi_index::indexed_by
    <
        boost::multi_index::hashed_non_unique //query by context
        <
            boost::multi_index::member<AOG_LIB::CM_Entry<T>, ContextType<T>, &AOG_LIB::CM_Entry<T>::context_>
        >,
        boost::multi_index::hashed_non_unique //query by configuration
        <
            boost::multi_index::member<AOG_LIB::CM_Entry<T>, std::vector<SequenceType <T> >, &AOG_LIB::CM_Entry<T>::configuration_>
        >,
	    boost::multi_index::hashed_non_unique //query by count, may not be necessary, remove if not
	    <
		    boost::multi_index::member<AOG_LIB::CM_Entry<T>, unsigned , &AOG_LIB::CM_Entry<T>::count_>
        >,
        boost::multi_index::random_access<> //enable vector like access, .size() can get memory usage
    >
>;

/*
 * naive implementation of a table using standard library
 * remove after prove the boost multi_index works
 *
 *
template <class StateType>
using Config2Context = std::unordered_map
<
    SequenceType<StateType>,
    std::unordered_map
    <
        ContextType<StateType>, std::shared_ptr<CM_Entry<StateType> >
    >
>;
template <class StateType>
using Context2Config = std::unordered_map
<
    ContextType<StateType>,
    std::unordered_map
    <
        SequenceType<StateType>, std::shared_ptr<CM_Entry<StateType> >
    >
>;

//Define a multi index container so that context matrix can be
//accessed by both configuration and context
template <class T>
struct Context_Matrix
{
    Config2Context<T> config_context_;
    Context2Config<T> context_config_;
};
*/
#endif //AOG_LIB_CONTEXT_MATRIX_H
