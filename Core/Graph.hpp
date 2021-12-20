//
// Creator:          Feng Shi
//
// Date of Creation: 8/16/17.
//
// Description:      This class implements a lightweight graph library;
//                   the graph can be in directed or undirected structure.

#ifndef GRAPH_GRAPH_HPP
#define GRAPH_GRAPH_HPP

#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>
#include <utility>
#include <vector>
#include <queue>
#include <list>
#include <map>
#include<unordered_map>

namespace AOG_LIB
{

    typedef unsigned int VertexId;

    template<class TV, class TE>
    class Graph
    {
    private:
        class Vertex
        {
        public:
            VertexId id_;

            std::shared_ptr<TV> attrs_;

            explicit Vertex(std::shared_ptr<TV>);

            Vertex(VertexId, std::shared_ptr<TV>);
            Vertex(const Vertex &other);
            ~Vertex();

            std::shared_ptr<TV> GetAttributes() const;

            void SetAttributes(std::shared_ptr<TV>);

            VertexId GetId() const;

            void SetId(VertexId);
        };

        class Edge
        {
        public:
            VertexId source_;  // the index or id of source vertex

            VertexId target_;  // the index or id of target vertex

            std::shared_ptr<TE> attrs_;

            // constructor: we suppose that T can be null
            Edge(VertexId, VertexId, std::shared_ptr<TE>);
            Edge(const Edge &other);
            // destructor
            ~Edge();

            // get the value attribute
            std::shared_ptr<TE> GetAttributes() const;

            // set the value attribute
            void SetAttributes(std::shared_ptr<TE>);

            // the source of the edge, cannot be modified
            VertexId source() const;

            // the target of the edge, cannot be modified
            VertexId target() const;

            // toString function
            std::ostream &operator<<(std::ostream &os);
        };

        /*Graph member variables*/
        unsigned number_vertices_;       // number of vertices

        unsigned number_edges_;          // number of edges

        bool is_directed_;           // is directed graph?

        bool reuse_memory_;          //whether to use a more memory efficient representation?

        std::vector<std::shared_ptr<Vertex> > vertices_;

        std::vector<std::vector<std::shared_ptr<Edge> > > adjacency_lists_;

        std::vector<std::vector<VertexId> > reverse_adjacency_lists_;

        std::queue<VertexId> empty_spots_;//a queue that keeps track of empty spots
        //Vertex class inside

    public:
        //default constructor
        Graph();
        

        // @constructor:
        //    unsigned: number of vertices
        //    bool: is directed-graph?
        //    bool: use a more memory efficient representation for vertices and adjacency list ?
        Graph(unsigned, bool, bool=false);

        //    std::istream: construct a graph from input, e.g. file
        Graph(std::istream &, bool);

        // copy constructor for Graph
        Graph(const Graph& other);
        // @destructor
        ~Graph();

        // number of vertices
        virtual unsigned NumberOfVertices() const;

        // number of edges
        unsigned NumberOfEdges() const;

        // test if the graph is directed
        bool IsDirected() const;

        //check if the vertex id is valid
        bool IsValidVertex(VertexId) const;

        // add a vertex to graph
        virtual VertexId AddVertex(std::shared_ptr<TV>);

        // add an edge to graph
        // return success of the addition
        virtual bool AddEdge(VertexId, VertexId, std::shared_ptr<TE>, bool=false,  int = -1);

        //remove a vertex given its vertex id
        virtual bool DeleteVertex(VertexId);

        //remove all edges that have given source vertex id and target vertex id
        bool DeleteEdge(VertexId, VertexId);

        //get children nodes of a vertex given its id
        //based on the edge's source and target
        std::vector<VertexId> ChildrenVertices(VertexId) const;
      
        //get parent nodes of a vertex given its id
        //based on the edge's source and target
        std::vector<VertexId> ParentsVertices(VertexId) const;

        //get edges of a vertex as target given its id
        std::vector<std::pair<VertexId, VertexId> > InEdges(VertexId) const;

        //get edges of a vertex as source given its id
        std::vector<std::pair<VertexId, VertexId> > OutEdges(VertexId) const;

        //get attribute of a vertex given its vertex id
        std::shared_ptr<TV> GetVertexContent(VertexId) const;

        //get attribute of an edge given its source and target
        std::vector<std::shared_ptr<TE> > GetEdgeContent(VertexId, VertexId) const;

        // return all edges
        std::vector<std::pair<VertexId, VertexId> > AllEdges() const;

        // return an unordered map containing parent VertexId (key) and a vector of child VertexId (value),
        // and output the graph in JSON format in file graph.json given an destination directory
        std::unordered_map<VertexId, std::vector<VertexId> > PrintGraphJSON(std::string="");
    };


/*templated functions must be in header file*/


/*Vertex Functions*/
//
// Created by Feng Shi on 8/16/17.
//

    template<class TV, class TE>
    Graph<TV, TE>::Vertex::Vertex(const std::shared_ptr<TV> attr)
    {
        this->attrs_ = attr;
    }

    template<class TV, class TE>
    Graph<TV, TE>::Vertex::Vertex(const VertexId id,
                                  const std::shared_ptr<TV> attrs)
            :id_(id), attrs_(attrs)
    {}

    template<class TV, class TE>
    Graph<TV, TE>::Vertex::Vertex(const Vertex &other)
            
    {
        
        this->id_ = other.id_;
        this->attrs_ = std::make_shared<TV>(*(other.attrs_));
    }
    

    template<class TV, class TE>
    Graph<TV, TE>::Vertex::~Vertex()
    {
        // std::cout << "deleting Vertex(" << this->id_ << ")" << std::endl;
    }

    template<class TV, class TE>
    std::shared_ptr<TV> Graph<TV, TE>::Vertex::GetAttributes() const
    {
        return this->attrs_;
    }

    template<class TV, class TE>
    void Graph<TV, TE>::Vertex::SetAttributes(const std::shared_ptr<TV> newVal)
    {
        this->attrs_ = newVal;
    }

    template<class TV, class TE>
    void Graph<TV, TE>::Vertex::SetId(VertexId id)
    {
        this->id_ = id;
    }

    template<class TV, class TE>
    VertexId Graph<TV, TE>::Vertex::GetId() const
    {
        return this->id_;
    }

/*Edge Functions*/
//
// Created by Feng Shi on 8/16/17.
//
    template<class TV, class TE>
    Graph<TV, TE>::Edge::Edge(const VertexId source,
                              const VertexId target,
                              const std::shared_ptr<TE> attrs)
            : source_(source), target_(target), attrs_(attrs)
    {}

    template<class TV, class TE>
    Graph<TV, TE>::Edge::Edge(const Edge &other)
            : source_(other.source_), target_(other.target_)
    {
        this->attrs_ = std::make_shared<TE>(*(other.attrs_));
    }

    template<class TV, class TE>
    Graph<TV, TE>::Edge::~Edge()
    {
        // std::cout << "Edge between Vertex(" << this->source_
        //           << ") and Vertex(" << this->target_
        //           << ") has been deleted" << std::endl;
    }

    template<class TV, class TE>
    std::shared_ptr<TE> Graph<TV, TE>::Edge::GetAttributes() const
    {
        return this->attrs_;
    }

    template<class TV, class TE>
    void Graph<TV, TE>::Edge::SetAttributes(const std::shared_ptr<TE> newVal)
    {
        this->attrs_ = newVal;
    }

    template<class TV, class TE>
    VertexId Graph<TV, TE>::Edge::source() const
    { return this->source_; }

    template<class TV, class TE>
    VertexId Graph<TV, TE>::Edge::target() const
    { return this->target_; }

    template<class TV, class TE>
    std::ostream &Graph<TV, TE>::Edge::operator<<(std::ostream &os)
    {
        return os << "(" << this->source_ << " -> " << this->target_ << ")";
    }

/*Graph Functions*/

    template<class TV, class TE>
    Graph<TV, TE>::Graph()
            :number_edges_(0), number_vertices_(0), is_directed_(true),reuse_memory_(false)
    {}

    template<class TV, class TE>
    Graph<TV, TE>::Graph(unsigned nv, bool directed,bool reuse):
            number_edges_(0), adjacency_lists_(nv),
            vertices_(nv), reverse_adjacency_lists_(nv),reuse_memory_(reuse)
    {
        this->number_vertices_ = nv;
        this->is_directed_ = directed;

        for (unsigned i = 0; i < nv; i++)
        {
            this->adjacency_lists_[i] = std::vector<std::shared_ptr<Edge> >();
            this->reverse_adjacency_lists_[i] = std::vector<VertexId>();
            std::shared_ptr<TV> temp_ptr = 0;
            this->vertices_[i] = std::make_shared<Vertex>(i, temp_ptr);
        }
    }

    template<class TV, class TE>
    Graph<TV,TE>::Graph(const Graph& other):number_vertices_(other.number_vertices_),
                                            number_edges_(other.number_edges_),
                                            is_directed_(other.is_directed_),
                                            reuse_memory_(other.reuse_memory_),
                                            reverse_adjacency_lists_(other.reverse_adjacency_lists_),
                                            empty_spots_(other.empty_spots_)

    {
         
        for(const auto& vtx_ptr : other.vertices_)
        {
            if (vtx_ptr)       
                this->vertices_.push_back(std::make_shared<Vertex>(*vtx_ptr));
            else
                this->vertices_.push_back(0);
        }

        for(int i = 0; i < other.adjacency_lists_.size(); ++i)
        {
            this->adjacency_lists_.push_back(std::vector<std::shared_ptr<Edge> >());
            for(int j = 0; j < other.adjacency_lists_[i].size();++j)
            {
                std::shared_ptr<Edge> edge_ptr = other.adjacency_lists_[i][j];
                this->adjacency_lists_[i].push_back(std::make_shared<Edge>(*edge_ptr));
            }
        }


    }


    template<class TV, class TE>
    Graph<TV, TE>::Graph(std::istream &is, bool directed)
    {
        //TODO read in data from file "is"
        this->is_directed_ = directed;
    }

    template<class TV, class TE>
    Graph<TV, TE>::~Graph()
    {
        this->vertices_.clear();
        this->adjacency_lists_.clear();
        this->reverse_adjacency_lists_.clear();
    }

    template<class TV, class TE>
    unsigned Graph<TV, TE>::NumberOfVertices() const
    {
        return this->number_vertices_;
    }

    template<class TV, class TE>
    unsigned Graph<TV, TE>::NumberOfEdges() const
    {
        return this->number_edges_;
    }

    template<class TV, class TE>
    bool Graph<TV, TE>::IsDirected() const
    {
        return this->is_directed_;
    }

    template<class TV, class TE>
    bool Graph<TV, TE>::IsValidVertex(const VertexId vid) const
    {
        //in range and not been deleted
        return (vid < this->vertices_.size() && this->vertices_[vid]);
    }

    template<class TV, class TE>
    std::vector<VertexId> Graph<TV, TE>::ChildrenVertices(VertexId vid) const
    {
        std::vector<VertexId> children(0);
        if(!this->IsValidVertex(vid))
            return children;
        for (const auto & outEdge: this->adjacency_lists_[vid])
            children.push_back(outEdge->target());
        return children;
    }

    template<class TV, class TE>
    std::vector<VertexId> Graph<TV, TE>::ParentsVertices(VertexId vid) const
    {
        if(!this->IsValidVertex(vid))
        {
            std::vector<VertexId> empty_v;
            return empty_v;
        }
        return this->reverse_adjacency_lists_[vid];
    }

    template<class TV, class TE>
    std::vector<std::pair<VertexId, VertexId> > Graph<TV, TE>::OutEdges(VertexId vid) const
    {
        std::vector<std::pair<VertexId, VertexId> > edges(0);
        if(!this->IsValidVertex(vid))
            return edges;
        for (const auto & outEdge: this->adjacency_lists_[vid])
            edges.push_back(std::make_pair(vid, outEdge->target()));
        return edges;
    }

    template<class TV, class TE>
    std::vector<std::pair<VertexId, VertexId> > Graph<TV, TE>::InEdges(VertexId vid) const
    {
        std::vector<std::pair<VertexId, VertexId> > edges(0);
        if(!this->IsValidVertex(vid))
            return edges;
        for (const auto & parent_vid: this->reverse_adjacency_lists_[vid])
            edges.push_back(std::make_pair(parent_vid, vid));
        return edges;
    }

    template<class TV, class TE>
    VertexId Graph<TV, TE>::AddVertex(const std::shared_ptr<TV> attr)
    {
        VertexId new_id;
        if (reuse_memory_ && !empty_spots_.empty())
        {
            new_id = empty_spots_.front();
            empty_spots_.pop();
            this->vertices_[new_id] = std::make_shared<Vertex>(new_id,attr);
        }
        else
        {
            new_id = unsigned(this->vertices_.size());
            this->vertices_.push_back(std::make_shared<Vertex>(new_id, attr));
            this->adjacency_lists_.push_back(std::vector<std::shared_ptr<Edge> >());
            this->reverse_adjacency_lists_.push_back(std::vector<VertexId>());
        }
        this->number_vertices_++;
        return new_id;
    }

    template<class TV, class TE>
    bool Graph<TV, TE>::AddEdge(const VertexId source,
                                const VertexId target,
                                const std::shared_ptr<TE> attr,bool multi_edge, int pos)
    {
        //check source and target validity
        if ((source >= this->vertices_.size() || target >= this->vertices_.size()) ||
            !(this->vertices_[source] && this->vertices_[target]))
        {
            std::cerr << "Invalid source or target vertex id\n";
            return false;
        }
        if(source == target)
        {
            std::cerr << "Source shouldn't be same with target\n";
	        std::cerr<<"The source is "<<source<<" and the target is "<<target<<std::endl;
            // std::throw exception();
            return false;
        }
        
        if(!multi_edge && std::find_if(this->adjacency_lists_[source].begin(),
                        this->adjacency_lists_[source].end(),
                        [target](const std::shared_ptr<Edge> edge) {
                            return edge->target() == target;
                        })
           != this->adjacency_lists_[source].end())
        {
	        std::cerr << "Edge already exists between Vertex " << source << " and Vertex " << target << std::endl;
            std::cerr<<"the edge: "<< source <<"-> "<<target<<std::endl;
            throw std :: exception();
            return false;
        }

        std::shared_ptr<Edge> newEdge(new Edge(source, target, attr));
       
        //if add edge at the end of adjacency list(default behavior)
        if(pos == -1)
        {
            this->adjacency_lists_[source].push_back(newEdge);
        }
        else
        {
            auto iter = adjacency_lists_[source].begin() + pos;
            this->adjacency_lists_[source].insert(iter,newEdge);
        }
       
        this->reverse_adjacency_lists_[target].push_back(source);
        this->number_edges_++;
        return true;
    }


    template<class TV, class TE>
    bool Graph<TV, TE>::DeleteVertex(const VertexId vid)
    {
        if (vid >= this->vertices_.size() || !this->vertices_[vid])
        {
            std::cerr << "Invalid vertex to delete\n";
            return false;
        }
        this->number_vertices_--;
        this->vertices_[vid] = 0;

        //remember the empty spot for later use
        if(reuse_memory_)
            empty_spots_.push(vid);

        //remove children's edges
        for (const auto & edge: this->adjacency_lists_[vid])
        {
            //this is the erase-remove idiom in C++
            this->reverse_adjacency_lists_[edge->target()].erase(
                    std::remove(this->reverse_adjacency_lists_[edge->target()].begin(),
                                this->reverse_adjacency_lists_[edge->target()].end(), vid),
                    this->reverse_adjacency_lists_[edge->target()].end());
        }
        this->number_edges_ -= unsigned(this->adjacency_lists_[vid].size());
        this->adjacency_lists_[vid].clear();

        //remove parent's edges
        for (const auto & parent: this->reverse_adjacency_lists_[vid])
        {
            //erase-remove idiom in C++ using remove if
            this->adjacency_lists_[parent].erase(std::remove_if(this->adjacency_lists_[parent].begin(),
                                                                this->adjacency_lists_[parent].end(),
                                                                [vid](const std::shared_ptr<Edge> edge) {
                                                                    return edge->target() == vid;
                                                                }),
                                                 this->adjacency_lists_[parent].end());
        }
        this->number_edges_ -= unsigned(this->reverse_adjacency_lists_[vid].size());
        this->reverse_adjacency_lists_[vid].clear();
        return true;
    }

    template<class TV, class TE>
    bool Graph<TV, TE>::DeleteEdge(VertexId source, VertexId target)
    {
        if ((source >= this->vertices_.size() || target >= this->vertices_.size()) ||
            !(this->vertices_[source] && this->vertices_[target]))
        {
            std::cerr << "Invalid source or target vertex id\n";
            return false;
        }
        else if (std::find_if(this->adjacency_lists_[source].begin(),
                              this->adjacency_lists_[source].end(),
                              [target](const std::shared_ptr<Edge> edge) {
                                  return edge->target() == target;
                              })
                 == this->adjacency_lists_[source].end())
        {
            std::cerr << "Edge from " << source << " to " << target << " doesn't exist (not an error in AddRule() / ReconstructTree())" << std::endl;
            return false;
        }
        else
        {
            //remove from children's side
            this->reverse_adjacency_lists_[target].erase(
                    std::remove(this->reverse_adjacency_lists_[target].begin(),
                                this->reverse_adjacency_lists_[target].end(), source),
                    this->reverse_adjacency_lists_[target].end());

            //remove from parent's side
            this->adjacency_lists_[source].erase(std::remove_if(this->adjacency_lists_[source].begin(),
                                                                this->adjacency_lists_[source].end(),
                                                                [target](const std::shared_ptr<Edge> edge) {
                                                                    return edge->target() == target;
                                                                }),
                                                 this->adjacency_lists_[source].end());
            this->number_edges_--;
            return true;
        }
    }

    template<class TV, class TE>
    std::shared_ptr<TV> Graph<TV, TE>::GetVertexContent(const VertexId vid) const
    {
        if (vid >= this->vertices_.size() || !this->vertices_[vid])
            return 0;
        return this->vertices_[vid]->GetAttributes();
    }

    template<class TV, class TE>
    std::vector<std::shared_ptr<TE> > Graph<TV, TE>::GetEdgeContent(const VertexId source, const VertexId target) const
    {
	    std:: vector<std::shared_ptr<TE> > edges;
        if ((source >= this->vertices_.size() || target >= this->vertices_.size()) ||
                !(this->vertices_[source] && this->vertices_[target]))
        {
            std::cerr << "Invalid source or target vertex id\n";
            return edges;
        }

        for (const auto & edge: this->adjacency_lists_[source])
        {
            if (edge->target() == target)
                edges.push_back(edge->GetAttributes());
        }

        if(edges.empty())
        {
            std::cerr << "Vertex " << source << "doesn't connect with vertex " << target << std::endl;
            return edges;
        }
	  
	    return edges;
    }

    template<class TV, class TE>
    std::vector<std::pair<VertexId, VertexId> > Graph<TV, TE>::AllEdges() const
    {
        std::vector<std::pair<VertexId, VertexId> > edges;
        for (const auto & v_prt: this->vertices_)
        {
            if (v_prt)
            {
                for (const auto & edge: this->adjacency_lists_[v_prt->GetId()])
                {
                    edges.push_back(std::make_pair(v_prt->GetId(), edge->target()));
                }
            }
        }
        return edges;
    }


    /*
     * @param dir : the destination of the JSON file
     * @return : an unordered_map storing VertexId as keys and vectors of child VertexId as values
     */
    template<class TV, class TE>
    std::unordered_map<VertexId, std::vector<VertexId> > Graph<TV, TE>::PrintGraphJSON(std::string dir)
    {
        // graph_map is an auxiliary map for printing the JSON file
        std::map<VertexId, std::vector<VertexId> > graph_map;
        // graph_umap is the actual map returned
        std::unordered_map<VertexId, std::vector<VertexId> > graph_umap;

        // empty graph
        if (this->number_vertices_ == 0)
        {
            std::cout << "This graph is empty" << std::endl;
            return graph_umap;
        }

        // construct maps with the parent VertexId as key and a vector of child VertexId as value
        for (const auto & v_prt: this->vertices_)
        {
            if (v_prt)
            {
                VertexId source = v_prt->GetId();
                std::vector<VertexId> targets;
                for (const auto & edge: this->adjacency_lists_[source])
                    targets.push_back(edge->target());
                graph_map.insert(std::make_pair(source, targets));
                graph_umap.insert(std::make_pair(source, targets));
            }
        }

        std::ofstream file;
        if (dir.empty())
            dir = std::string("graph.json");
        else
            dir = dir.append("/graph.json");

        // output the graph in the following format:
        // {
        //      "parent1" : [child1, child2],
        //      "parent2" : [child1, child2]
        // }
        file.open(dir);
        if (file.is_open()) {
            file << "{\n\t";
            for (auto it: graph_map) {
                VertexId source = it.first;
                std::vector<VertexId> target = it.second;
                file << "\"" << source << "\":[";
                for (auto t: target) {
                    file << t;
                    if (t != target.back())
                        file << ", ";
                }
                file << "]";
                if (it.first != graph_map.rbegin()->first)
                    file << ",\n\t";
                else
                    file << "\n";
            }
            file << "}\n";
            file.close();
        }
        else
            std::cout << "Unable to open " << dir << std::endl;

        return graph_umap;
    }
}
#endif //GRAPH_GRAPH_HPP
