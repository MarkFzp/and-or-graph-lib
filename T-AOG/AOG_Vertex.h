//
// Created by yuanluyao on 10/23/17.
//

#ifndef AOG_LIB_AOG_VERTEX_H
#define AOG_LIB_AOG_VERTEX_H

#include <string>
#include "Symbolic_State.h"
#include "../Core/Graph.hpp"
namespace AOG_LIB
{
    //forward declaration
    template <class StateType> class T_AOG;
    class AOG_Edge;

    /* This class defines vertices in an AOG graph.
     * Mutators in this class are set as private to prevent unwanted changes of vertices' properties from users
     * This class declares T_AOG as friend class to enable T_AOG to change its vertices' properties.
     */
    template<class StateType>
    class AOG_Vertex
    {
        bool is_and_; // true if the AOG Vertex is an And-node
        bool is_root_;// true if the AOG Vertex is a root node
        Symbolic_State<StateType> state_; // The Symbolic state contained in the AOG Vertex


        //friend class Graph<AOG_Vertex<StateType>,AOG_Edge>;
        friend class T_AOG<StateType>;

        /* This function sets the symbolic state of a AOG Vertex.
         * @param:
         *      state: The symbolic state to be set
         * @return: void
         */
        void SetState(const Symbolic_State<StateType> &state){this->state_ = state;};
        /* This function sets a vertex to be an and-node or an or-node in AOG graph
         * @param:
         *      is_and: a boolean value that specifies whether the vertex is an and-node
         * @return: void
         */
        void SetIsAnd(const bool is_and){this->is_and_ = is_and;};
        void SetIsRoot(bool is_root){this->is_root_ = is_root;};
    public:
        /* Construct an AOG Vertex from a symbolic state along with the node's properties.
         * @param:
         *      state: the symbolic state used to construct the vertex
         *      bool is_and: specify whether the vertex is an and-node or an or-node
         *      bool is_root: specify whether the vertex is a root node
         */
        AOG_Vertex(const Symbolic_State<StateType> &, bool, bool);

        // return the current symbolic state of the AOG vertex
        const Symbolic_State<StateType>& GetState();

        // return a boolean value indicating whether the vertex is an and-node.
        bool IsAnd() const;

        // return a boolean value indicating whether the vertex is a root node.
        bool IsRoot() const;

    };


    template<class StateType>
    AOG_Vertex<StateType>::AOG_Vertex(const Symbolic_State<StateType> &state, bool is_and, bool is_root)
            :state_(state), is_and_(is_and), is_root_(is_root)
    {}

    template<class StateType>
    const Symbolic_State<StateType>& AOG_Vertex<StateType>::GetState() 
    { return this->state_; }

    template<class StateType>
    bool AOG_Vertex<StateType>::IsAnd() const
    { return this->is_and_; }

    template<class StateType>
    bool AOG_Vertex<StateType>::IsRoot() const
    { return this->is_root_; }
    
}
#endif //AOG_LIB_AOG_VERTEX_H
