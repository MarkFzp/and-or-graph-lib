#include <vector>
#include <string>
#include <unordered_map>

#include "./T-AOG/T_AOG.h"
#include "./Learner/Online_Learner.h"

using namespace std;
using namespace AOG_LIB;
//using ContextType = std::pair<SequenceType<T>, SequenceType<T>>;

int main()
{
    string a11_str = "c1";
    string a12_str = "c4";
    string a21_str = "c3";
    string a22_str = "c2";
    Symbolic_State<string> a11(a11_str, false);
    Symbolic_State<string> a12(a12_str, false);
    Symbolic_State<string> a21(a21_str, false);
    Symbolic_State<string> a22(a22_str, false);

    vector<Symbolic_State<string> > firstChild = {a11, a12};
    vector<Symbolic_State<string> > secondChild = {a21, a22};

    ConfigBuffer<string> or_children = {firstChild, secondChild};

    // struct AOFStruct
    // {
    //     ConfigBuffer<StateType> or_children_;
    //     std::vector<std::vector<int>> weights_;
    //     AOG_LIB::Symbolic_State<StateType> and_node_;
    //     static int id_;

    vector<int> weight1_col(4,10);
    vector<vector<int> > weight;
    weight.push_back(weight1_col);
    Symbolic_State<string> temp(a11_str, false);
    AOFStruct<string> AOF = {
        or_children,
        weight,
        temp,
    };


    //this->CM_.insert({context, configuration, 1});
    //configurations in the CM
    SequenceType<string> first_config = {a11, a21};
    SequenceType<string> second_config = {a11, a22};
    SequenceType<string> third_config = {a12, a21};
    SequenceType<string> fourth_config = {a21, a11};

    //contexts in the CM
    SequenceType<string> context1_1 = {a11};
    SequenceType<string> context1_2 = {a22, a21};
    SequenceType<string> context2_1 = {a22};
    SequenceType<string> context2_2 = {a21, a12};
    SequenceType<string> context3_1 = {a12};
    SequenceType<string> context3_2 = {a12, a22};

    ContextType<string> context1 = make_pair(context1_1, context1_2);
    ContextType<string> context2 = make_pair(context2_1, context2_2);
    ContextType<string> context3 = make_pair(context3_1, context3_2);

    Context_Matrix<string> CM;
    CM.insert({context1, first_config, 1});
    CM.insert({context1, second_config, 0});
    CM.insert({context1, third_config, 0});
    CM.insert({context1, fourth_config, 3});

    CM.insert({context2, first_config, 4});
    CM.insert({context2, second_config, 1});
    CM.insert({context2, third_config, 5});
    CM.insert({context2, fourth_config, 0});

    CM.insert({context3, first_config, 0});
    CM.insert({context3, second_config, 2});
    CM.insert({context3, third_config, 3});
    CM.insert({context3, fourth_config, 1});

    //Reduction_full
    //using RdBuffer = std::unordered_map<SequenceType<T>, std::vector<std::pair<unsigned, unsigned> > >;
    pair<unsigned, unsigned> dummy(0,1);
    vector<pair<unsigned, unsigned> > rd_first_pos;
    vector<pair<unsigned, unsigned>> rd_second_pos;
    vector<pair<unsigned, unsigned>> rd_third_pos;
    vector<pair<unsigned, unsigned>> rd_fourth_pos;
    
    rd_first_pos.push_back(dummy);
    rd_first_pos.push_back(dummy);
    rd_first_pos.push_back(dummy);
    rd_first_pos.push_back(dummy);
    rd_first_pos.push_back(dummy);

    rd_second_pos.push_back(dummy);
    rd_second_pos.push_back(dummy);
    rd_second_pos.push_back(dummy);

    rd_third_pos.push_back(dummy);
    rd_third_pos.push_back(dummy);
    rd_third_pos.push_back(dummy);
    rd_third_pos.push_back(dummy);
    rd_third_pos.push_back(dummy);
    rd_third_pos.push_back(dummy);
    rd_third_pos.push_back(dummy);
    rd_third_pos.push_back(dummy);

    rd_fourth_pos.push_back(dummy);
    rd_fourth_pos.push_back(dummy);
    rd_fourth_pos.push_back(dummy);
    rd_fourth_pos.push_back(dummy);

    RdBuffer<string> reduction_full;
    reduction_full.insert({first_config, rd_first_pos});
    reduction_full.insert({second_config, rd_second_pos});
    reduction_full.insert({third_config, rd_third_pos});
    reduction_full.insert({fourth_config, rd_fourth_pos});

    double likelihood_gain = CalculateLikelihoodGain(AOF, reduction_full, CM);
    cout << likelihood_gain;

    return 0;
}