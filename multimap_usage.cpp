#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/random_access_index.hpp>

#include "T-AOG/Symbolic_State.h"
#include "Learner/Online_Learner.h"

/*
 * The reason this file stays in our noble folder is to be an example
 * of using boost::multi_index_container. As long as we have usage
 * in our code, please delete this file immediately.
 */

using namespace boost::multi_index;

typedef std::vector<AOG_LIB::Symbolic_State<std::string> > CodeType;

struct animal
{
    std::string name;
    CodeType code;
    int legs;
};

typedef multi_index_container<
        animal,
        indexed_by<
                hashed_non_unique<
                        member<
                                animal, std::string, &animal::name
                        >
                >,
                hashed_non_unique<
                        member<
                                animal, int, &animal::legs
                        >
                >,
                hashed_non_unique<
                        member<
                                animal, CodeType, &animal::code
                        >
                >,
                random_access<>
        >
> animal_multi;

int main()
{
    AOG_LIB::MCTS seacher;
    auto tg = std::make_shared<AOG_LIB::T_AOG<std::string> >();
    AOG_LIB::Online_Learner<std::string> ol(tg);
    std::string square_str = "square";
    std::string square_adjcnt_edges_str = "square adjacent edges";
    std::string square_sprt_edges_str = "square separate edges";
    std::string atom_edge_str = "edge";

    AOG_LIB::Symbolic_State<std::string> square(square_str,false);
    AOG_LIB::Symbolic_State<std::string> sqr_adjcnt_edges(square_adjcnt_edges_str,false);
    AOG_LIB::Symbolic_State<std::string> sqr_sprt_edges(square_sprt_edges_str,false);
    AOG_LIB::Symbolic_State<std::string> edge(atom_edge_str,true);

    CodeType code1 = {square, sqr_adjcnt_edges};
    CodeType code2 = {edge, edge};
    CodeType code3 = {sqr_adjcnt_edges, square};

    animal_multi animals;

    animals.insert({"cat", code1, 4});
    animals.insert({"shark", code2, 0});
    animals.insert({"pig", code3, 4});
    animals.insert({"spider", code1, 8});
    animals.insert({"cat", code2, 9});

    auto &legs_index = animals.get<1>();
    auto its = legs_index.equal_range(4);
    for(; its.first != its.second; its.first++)
    {
        std::cout << its.first->name << std::endl;
    }

    auto &name_index = animals.get<0>();
    auto itts = name_index.equal_range("cat");
    for(; itts.first != itts.second; itts.first++)
    {
        std::cout << itts.first->legs << std::endl;
    }

    auto &code_index = animals.get<2>();
    auto ittts = code_index.equal_range(code1);
    for(; ittts.first != ittts.second; ittts.first++)
    {
        std::cout << ittts.first->name << std::endl;
    }

    auto & random_index = animals.get<3>();
    unsigned len = random_index.size();
    std::cout << "length is: " << len << std::endl;
}
