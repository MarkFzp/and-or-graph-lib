# Temporal And-Or Graph(T-AOG) Library
This is a C++ implementation of a T-AOG library that can be used to model decomposition of objects and events.

## Requirements
* C++ Boost Libraries
* CMake is useful for building

## Instructions
To use this library, first install C++ boost libraries. You can download it from <http://www.boost.org/>, or by your system's package manager. For example, if you are using macOS, run the following command in your terminal:  
`brew install boost`

For CMake, you can download it from <https://cmake.org/download/>, or by your system's package manager, for example:  
`brew install cmake`

Then git clone this repository:  
`git clone https://github.com/vcla/aog_lib.git`

And change to the directory:  
`cd aog_lib/`

Now create a directory for running CMake:  
`mkdir cmake-build-debug`

And change to the newly created folder:  
`cd cmake-build-debug/`

After that, in the "cmake-build-debug" folder, run:  
`cmake ../`

And run:  
`make`

to compile the library. Finally run:  
`./graph`

to start the main sample program.

## AOG_LIB :: T_AOG
### Constructors
Constructors | Description
------------|--------------
T_AOG() | Default Constructor
T_AOG(Symbolic_Rule\<StateType\> rules) | Construct a T_AOG from a set of symbolic rules
T_AOG(Symbolic_State\<StateType\> leaf_states) | Construct a T_AOG from a set of leaf symbolic states. No edge will be added to the T_AOG
  
### Methods
Methods | Description
--------|--------------
unsigned NumOfRules() | Return the number of rules contained in the T_AOG.
unsigned NumOfStates() | Return the number of symbolic states in the T_AOG.
unsigned NumofLeafStates() | Return the number of leaf symbolic states in T_AOG.
bool ExistRule(Symbolic_Rule\<StateType\> rule) | Return true if the given rule exists in the T_AOG. Note that A->BC and A->CB are considered different rules
vector\<Symbolic_Rule \<StateType\> \> GetRules() | Return all the rules contained in the T_AOG.
vector\<Symbolic_State \<StateType\> \> GetStates() | Return all symbolic states in the T_AOG.
vector\<Symbolic_State \<StateType\> \> GetLeafStates() | Return all leaf symbolic states in the T_AOG.
VertexId GetVertexIdByState(Symbolic_State\<StateType\> state) | Return the vertex id of the given state.
Symbolic_State\<StateType\> GetStateByVertexId(VertexId source) | Return the symbolic state that has the given vertex id
VertexId AddVertex(shared_ptr\<AOG_Vertex\<StateType\>\> aog_vertex) | Add the given vertex into T_AOG and return the vertex id assigned to this vertex.
bool AddEdge(VertexId source, VertexId target, shared_ptr\<AOG_Edge\> aog_edge,bool multi_edge = true)| parameters:<br/> &nbsp;&nbsp;&nbsp;&nbsp;source: the vertex id of the source of the edge<br/>&nbsp;&nbsp;&nbsp;&nbsp;target: the vertex id of the target of the edge<br/>&nbsp;&nbsp;&nbsp;&nbsp; aog_edge: an AOG_Edge object that contains the weight of the edge to be added
bool DeleteVertex(VertexId vid) | Return true if the vertex indicated by the given vertex id is successfully deleted
void SetIsAnd(VertexId source_id, bool is_and) | set a vertex to an And-node or an Or-node
bool SetOutEdgeWeights(vertexId source, unordered_map\<VertexId ,double\> weights) | Reset the weights of all the outedges of a given Or-node. Return true if weights are set successfully
unordered_map\<VertexId,double\> GetOutEdgeWeights(VertexId source, bool is_normalized) | return an unordered_map that maps all the ids of outedges' target vertices to their corrensponding weights. The weights are normalized if is_normalized = true
unordered_map\<VertexId,double\> Normalize(VertexId src_id) | Normalize the branching probabilities of the given Or-node, and return the normalized weights.
unordered_map\<VertexId, vector\<VertexId\>\> Sample(VertexId root) | Sample a parse graph of the T_AOG starting from the given root vertex. 
void AddRule(Symbolic_Rule\<StateType\> rule) | Add a rule to the T_AOG and update it.
bool DeleteRule(Symbolic_Rule\<StateType\> rule) | Delete a rule in T_AOG and update it. Return false if the given rule does not exist.
void Visualize(string="aog", string="") | Visualize the T_AOG. See aogvisualizer folder for detailed usage.
