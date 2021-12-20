//
// Created by Luyao Yuan on 17/10/23.
//

#ifndef AOG_LIB_SYMBOLIC_STATE_H
#define AOG_LIB_SYMBOLIC_STATE_H

#include <functional>
#include <boost/functional/hash.hpp>
#include <iostream>
#include <fstream>
#include <string>

namespace AOG_LIB
{
/*
     * This class defines symbolic states, which are the semantic meanings of each AOG vertex in AOG graph.
     */
template <class StateType>
class Symbolic_State
{
    StateType state_content_; // the content of the symbolic state
    bool is_basic_;           // true if the symbolic state is used in a leaf node
    int id_;                  //if the state has semantic meaning, then -1, else other integers used to distinguish states with no semantic meaning

  public:
    Symbolic_State();
    //constructor for state with no actual semantic meaning
    Symbolic_State(int);
    Symbolic_State(int id, const StateType &);
    /* constructs a symbolic state with a given state and nodes' properties
         * @param:
         *      content: the semantic meaning of the state
         *      is_leaf: true if the symbolic state is used in a leaf node
         */
    Symbolic_State(const StateType &, bool);

    //true if the symbolic state is used in a leaf node.
    bool GetIsBasic() const;

    //return the content(the semantic meaning) of the symbolic state.
    const StateType &GetContent() const;

    int GetId() const { return this->id_; };

    /* This function compares two symbolic states
         * @params: rhs: the other symbolic state to be compared with
         * @return: true if the content, is_basic of the two symbolic states are all the same
         */
    bool operator==(const Symbolic_State<StateType> &) const;
    bool operator!=(const Symbolic_State<StateType> &) const;
};
template <class StateType>
size_t hash_value(const AOG_LIB::Symbolic_State<StateType> &state)
{

    using std::string;

    // Compute individual hash values for first,
    // second  and combine them using XOR
    // and bit shifting:
    boost::hash<StateType> hash_state;
    boost::hash<bool> hash_bool;
    return hash_state(state.GetContent()) ^ (hash_bool(state.GetIsBasic()) << 1);
    // return ((boost::hash_value(state.GetContent())
    // 	     ^ (boost::hash_value(state.GetIsBasic()) << 1)));
}
}

namespace std
{
template <class StateType>
struct hash<AOG_LIB::Symbolic_State<StateType> >
{
    size_t operator()(const AOG_LIB::Symbolic_State<StateType> &state) const noexcept
    {
        boost::hash<AOG_LIB::Symbolic_State<StateType>> hasher;
        return hasher(state);
    }
};
}

template <class StateType>
AOG_LIB::Symbolic_State<StateType>::Symbolic_State()
    : id_(-1), is_basic_(true)
{
}

template <class StateType>
AOG_LIB::Symbolic_State<StateType>::Symbolic_State(int id)
    : id_(id), is_basic_(false)
{
    /** std::cerr << "[CSCTR] " << id << std::endl; */
}

template <class StateType>
AOG_LIB::Symbolic_State<StateType>::Symbolic_State(int id, const StateType &content)
    : id_(id), is_basic_(false), state_content_(content)
{
    /** std::cerr << "[CSCTR] " << id << std::endl; */
}

template <class StateType>
AOG_LIB::Symbolic_State<StateType>::Symbolic_State(const StateType &content, bool is_leaf)
    : state_content_(content), is_basic_(is_leaf), id_(-1)
{
}

template <class StateType>
bool AOG_LIB::Symbolic_State<StateType>::GetIsBasic() const
{
    return this->is_basic_;
}

template <class StateType>
const StateType &AOG_LIB::Symbolic_State<StateType>::GetContent() const
{
    return this->state_content_;
}

template <class StateType>
bool AOG_LIB::Symbolic_State<StateType>::operator==(const Symbolic_State<StateType> &rhs) const
{
    return (this->is_basic_ == rhs.GetIsBasic() &&
            this->state_content_ == rhs.GetContent() &&
            this->id_ == rhs.GetId());
}

template <class StateType>
bool AOG_LIB::Symbolic_State<StateType>::operator!=(const Symbolic_State<StateType> &rhs) const
{
    return !(*this == rhs);
}

template <class StateType>
std::vector<std::vector<AOG_LIB::Symbolic_State<StateType> > > FileParser(std::string filename)
{
    std::ifstream file(filename);
	if(!file.good())
	{
		std::cout << "Cannot open file " << filename << std::endl;
		throw std::exception();
	}
    std::string line;
    std::vector<std::vector<AOG_LIB::Symbolic_State<StateType> > > ParsedVector;

    while (std::getline(file, line))
    {
        std::istringstream iss(line);

        std::vector<AOG_LIB::Symbolic_State<StateType> > subVector;
        do
        {
            std::string subs;
            iss >> subs;        
            AOG_LIB::Symbolic_State<StateType> subState(subs, true);
            subVector.push_back(subState);
        } while (iss);
        subVector.pop_back();
        ParsedVector.push_back(subVector);
    }
    return ParsedVector;
}

#endif //AOG_LIB_SYMBOLIC_STATE_H
