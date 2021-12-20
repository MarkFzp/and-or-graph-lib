//
// Created by yuanluyao on 10/25/17.
//

#ifndef AOG_LIB_LEARNER_H
#define AOG_LIB_LEARNER_H

#include "../T-AOG/T_AOG.h"

namespace AOG_LIB
{
    template<class StateType>
    class Learner
    {
    protected:
        std::shared_ptr<AOG_LIB::T_AOG<StateType> > related_graph_ptr_;
    public:
        Learner(std::shared_ptr<AOG_LIB::T_AOG<StateType> >);
        std::shared_ptr<AOG_LIB::T_AOG<StateType>> GetGraphPtr();
    };

    template <class StateType>
    Learner<StateType>::Learner(std::shared_ptr<AOG_LIB::T_AOG<StateType> > graph_ptr)
    {
        this->related_graph_ptr_ = graph_ptr;
    }

    template <class StateType>
     std::shared_ptr<AOG_LIB::T_AOG<StateType>> Learner<StateType> :: GetGraphPtr()
    {
        return this->related_graph_ptr_;
    }

}
#endif //AOG_LIB_LEARNER_H
