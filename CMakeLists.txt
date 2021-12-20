cmake_minimum_required(VERSION 3.8)
project(aog_lib)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost 1.58.0 COMPONENTS filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lboost_system -lboost_filesystem")

set(GRAPH_SOURCE_FILES
        Core/Graph.hpp
        T-AOG/AOG_Edge.h
        T-AOG/AOG_Vertex.h
        T-AOG/Symbolic_Rule.h
        T-AOG/T_AOG.h
        T-AOG/Symbolic_State.h)

find_package(Boost 1.58.0 COMPONENTS filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lboost_system -lboost_filesystem")

set(GRAPH_TEST_FILES
        Tests/test_graph.cpp)

set(LEANER_SOURCE_FILES
        Learner/Learner.h
        Learner/Context_Matrix.h
        Learner/Online_Learner.h
        Learner/Partial_Parser_Utils.h
        Learner/Checkpoint_Utils.h
        Learner/Checkpoint_Functions.h)

set(MCTS_SOURCE_FILES
        Utils/MCTS.h
        Utils/Earley_Parser.h)

set(ALGORITHM_TEST_FILES
        Tests/test_run.cpp)

add_executable(graph ${GRAPH_SOURCE_FILES} ${GRAPH_TEST_FILES})

add_executable(run ${LEARNER_SOURCE_FILES}
                   ${MCTS_SOURCE_FILES}
                   ${GRAPH_SOURCE_FILES}
                   ${ALGORITHM_TEST_FILES})

target_link_libraries(run
        ${Boost_LIBRARIES})
