//
// Created by yuanluyao on 10/23/17.
//

#ifndef AOG_LIB_AOG_EDGE_H
#define AOG_LIB_AOG_EDGE_H

namespace AOG_LIB
{
    //forward declaration
    template<class StateType> class AOG_Vertex;

    /* This class defines the edges in an AOG graph.
     * SetWeight(double) in this class is set as private to prevent unwanted changes of edge's weight from users.
     * This class declares T_AOG class as a friend class to enable T_AOG to change its edges' weights.
     */
    class AOG_Edge
    {
        double weight_; // the weight of an edge in AOG Graph

        //not sure if it is correct
        template <class StateType>
        friend class T_AOG;

        /* This function changes the weight of an AOG edge.
         * @param:
         *      weight: the weight to be set to.
         * @throw:
         *      if the given weight is less than zero, throw an exception.
         * @return: void
         */
        void SetWeight(double);

    public:
        // The default constructor that constructs an edge with weight 1.
        AOG_Edge():weight_(1.0){}
        // Construct an AOG edge with the given weight.
        explicit AOG_Edge(double);
        // return the current weight of the edge.
        double GetWeight() const;

    };


    AOG_Edge::AOG_Edge(double weight)
            : weight_(weight)
    {}

    double AOG_Edge::GetWeight() const
    { return this->weight_; }

    void AOG_Edge::SetWeight(double weight)
    {
        if(weight < 0)
        {
            std::cerr << "Weight cannot be less than zero\n";
            throw std::exception();
        }
        this->weight_ = weight;
    }
}
#endif //AOG_LIB_AOG_EDGE_H
