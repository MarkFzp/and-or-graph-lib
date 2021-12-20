#ifndef AOG_LIB_EARLEY_H
#define AOG_LIB_EARLEY_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <sstream>
#include <set>
#include <algorithm>
#include <utility>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <boost/functional/hash.hpp>

#include "../T-AOG/Symbolic_Rule.h"
#include "../T-AOG/Symbolic_State.h"


#ifndef CONTAINER_HASH
#define CONTAINER_HASH
template <typename Container> // (hasher for containers) we make this generic for any container
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};
#endif

namespace AOG_LIB
{
    // class EarleyParser(const grammar &); // Forward declaration

    // // Form the cartesian product of two sequences
    // template <class InputIt1, class InputIt2, class OutputIt>
    // OutputIt cartesian_product(InputIt1 begin1, InputIt1 end1, InputIt2 begin2, InputIt2 end2, OutputIt it)
    // {
    //     for (InputIt1 it1 = begin1; it1 != end1; ++it1) {
    //         for (InputIt2 it2 = begin2; it2 != end2; ++it2) {
    //             *it++ = std::make_pair(*it1, *it2);
    //         }
    //     }
    //     return it;
    // }
    
    // // Find out if container contains all of range of items
    // template <class AssociativeContainer, class InputIt>
    // bool find_all(const AssociativeContainer& container, InputIt begin, InputIt end)
    // {
    //     bool all = true;
    //     for (InputIt it = begin; it != end && all; ++it) {
    //         all = container.find(*it) != container.end();
    //     }
    //     return all;
    // }

    template<class StateType>
    class rule
    {
    public:
        typedef std::shared_ptr<rule> ptr;
        typedef std::shared_ptr<const rule> const_ptr;
        rule(const rule&) = delete;
        rule(){}
        rule& operator=(const rule&) = delete;
        // operator bool()
        // {
        //     return succeeded_;
        // }
        size_t size(unsigned int right) const
        {
            return right_[right].size();
        }
    private:
        Symbolic_State<StateType> left_;
        std::vector<std::vector<Symbolic_State<StateType> > > right_;
        // bool succeeded_;

        //  TODO: check whether "Symbolic_State<StateType>" or vector"Symbolic_State<StateType>"
        rule(const Symbolic_State<StateType>& left, const Symbolic_State<StateType>& right)
            : left_(left)
        {
            right_.emplace_back(1, right);
        }

        // TODO: check constructor
        rule(const Symbolic_Rule<StateType> &str)
        {
            left_ = str.GetSource();

            right_.push_back(str.GetResults());
        }
    public:
        // make right a single state vector
        static ptr create(const Symbolic_State<StateType> & left, const Symbolic_State<StateType>& right)
        {
            return ptr(new rule(left, right));
        }
        static ptr create(const Symbolic_Rule<StateType> & line)
        {
            return ptr(new rule(line));
        }
        const Symbolic_State<StateType>& left() const
        {
            return left_;
        }
        const std::vector<std::vector<Symbolic_State<StateType> > >& right() const
        {
            return right_;
        }
        friend std::ostream& operator <<(std::ostream& os, const rule& r)
        {
            os << r.left_.GetContent() << " -> ";
            unsigned int i = 0;
            for (const std::vector<Symbolic_State<StateType> >& alternative : r.right_) {
                for (const Symbolic_State<StateType> & symbol: alternative) {
                    os << symbol.GetContent() << " ";
                }
                if (i++ < r.right_.size() - 1) {
                    os << "| ";
                }
            }
            return os;
        }
    };
    
    template<class StateType>
    class grammar
    {
    public:
        // Create from a stream
        // grammar(std::istream& is)
        grammar(){}
        grammar(const std::vector<Symbolic_Rule<StateType> >& rules, const std::vector<Symbolic_State<StateType> >& top_level_rules)
        {
            tops_ = top_level_rules;
            // std::string line;
            // while (std::getline(is, line)) {
            //     line.erase(line.find_last_not_of("\r\n") + 1);
            //     if (line.size()) {
            //         rule::ptr r = rule::create(line);
            //         if (*r) {
            //             rules_.push_back(r);
            //         }
            //     }
            // }
            // std::cerr << "--------rules size inside grammar is: " << rules.size() << std::endl;
            for (int i = 0; i < rules.size(); i++){
                typename rule<StateType>::ptr r = rule<StateType>::create(rules[i]);
                // std::cerr << "Grammar initialization: rules.left is: " << r->left().GetContent() << std::endl;
                rules_.push_back(r);
            }
            // Get the terminals
            // TODO: debug here, changed set to unordered_set
            std::unordered_set<Symbolic_State<StateType> > nonterminals;
            std::unordered_set<Symbolic_State<StateType> > symbols;
            for (typename rule<StateType>::const_ptr r : rules_) {
                nonterminals.insert(r->left());
                for (const std::vector<Symbolic_State<StateType> >& alternative : r->right()) {
                    for (const Symbolic_State<StateType> & symbol : alternative) {
                        symbols.insert(symbol);
                    }
                }
            }
            for (const Symbolic_State<StateType> & symbol : symbols) {
                if (nonterminals.find(symbol) == nonterminals.end()) {
                    terminals_.push_back(symbol);
                }
            }
            // std::cerr << "private rules size inside grammar is: " << this->rules().size() << std::endl;
        }

        // grammar constructor used when doing delete;
        // The sequence could have and node or or node as terminals
        grammar(const std::vector<Symbolic_Rule<StateType>> &rules, const std::vector<Symbolic_State<StateType>> &top_level_rules, 
                const std::vector<Symbolic_State<StateType>> &added_terminals)
        {
            tops_ = top_level_rules;

            for (int i = 0; i < rules.size(); i++)
            {
                typename rule<StateType>::ptr r = rule<StateType>::create(rules[i]);
                // std::cerr << "Grammar initialization: rules.left is: " << r->left().GetContent() << std::endl;
                rules_.push_back(r);
            }
            // Get the terminals
            // TODO: debug here, changed set to unordered_set
            std::unordered_set<Symbolic_State<StateType>> nonterminals;
            std::unordered_set<Symbolic_State<StateType>> symbols;
            for (typename rule<StateType>::const_ptr r : rules_)
            {
                nonterminals.insert(r->left());
                for (const std::vector<Symbolic_State<StateType>> &alternative : r->right())
                {
                    for (const Symbolic_State<StateType> &symbol : alternative)
                    {
                        symbols.insert(symbol);
                    }
                }
            }
            for (const Symbolic_State<StateType> &symbol : symbols)
            {
                if (nonterminals.find(symbol) == nonterminals.end())
                {
                    terminals_.push_back(symbol);
                }
            }

            for (const Symbolic_State<StateType> added_terminal : added_terminals)
            {
                if (std::find(terminals_.begin(), terminals_.end(), added_terminal) == terminals_.end() && !added_terminal.GetIsBasic())
                    terminals_.push_back(added_terminal);
            }

            // std::cerr << "Print terminals at construction of grammar:\n";
            // print_terminal();
        }



        void update(const std::vector<Symbolic_Rule<StateType> >& rules, const std::vector<Symbolic_State<StateType> >& top_level_rules, const SequenceType<StateType> &input_seq)
        {
            this->rules_.clear();
            this->start_rules_.clear();
            this->terminals_.clear();
            this->tops_.clear();
            /** std::cerr << "Update Grammar" << std::endl; */
            tops_ = top_level_rules;

            for (int i = 0; i < rules.size(); i++){
                typename rule<StateType>::ptr r = rule<StateType>::create(rules[i]);
                rules_.push_back(r);
            }
            // Get the terminals
            // TODO: debug here, changed set to unordered_set
            std::unordered_set<Symbolic_State<StateType> > nonterminals;
            std::unordered_set<Symbolic_State<StateType> > symbols;
            for (typename rule<StateType>::const_ptr r : rules_) {
                nonterminals.insert(r->left());
                for (const std::vector<Symbolic_State<StateType> >& alternative : r->right()) {
                    for (const Symbolic_State<StateType> & symbol : alternative) {
                        symbols.insert(symbol);
                    }
                }
            }
            for (const Symbolic_State<StateType> & symbol : symbols) {
                if (nonterminals.find(symbol) == nonterminals.end()) {
                    terminals_.push_back(symbol);
                }
            }

            //add and/or nodes in the input_seq to the terminals_
            for(auto iter = input_seq.begin(); iter != input_seq.end(); iter++)
            {
                if(!(*iter).GetIsBasic())
                {
                    if(std::find(terminals_.begin(), terminals_.end(), (*iter)) == terminals_.end())
                        terminals_.push_back((*iter));
                }
            }
        }

        const std::vector<typename rule<StateType>::ptr>& rules() const
        {
            return rules_;
        }
    
        // Get all of the rules where symbol is the subject
        template <class OutputIt>
        OutputIt get_rules_for_left(const Symbolic_State<StateType> & symbol, OutputIt it) const
        {
            for (typename rule<StateType>::const_ptr r : rules_) {
                if (r->left() == symbol) {
                    *it++ = r;
                }
            }
            return it;
        }
    
        // Is this symbol a terminal (doesn't occur as the subject of a rule)?
        bool symbol_is_terminal(const Symbolic_State<StateType> & symbol) const
        {
            return std::find(terminals_.begin(), terminals_.end(), symbol) != terminals_.end();
        }

        void update_terminal(const Symbolic_State<StateType>& symbol)
        {
            terminals_.push_back(symbol);
        }

            // Get the rule(s) whose left-hand side is the start symbol
            template <class OutputIt>
            OutputIt get_start_rules(OutputIt it) const
        {
            Symbolic_State<StateType> start_symbol;
            bool started = false;
            for (typename rule<StateType>::const_ptr r : rules_) {
                if (!started || r->left() == start_symbol) {
                    *it++ = r;
                }
                if (!started) {
                    started = true;
                    start_symbol = r->left();
                }
            }
            return it;
        }
    
        // TODO: what this function does
        // Get the rules where any of these pairs of symbols occur together as a right hand side alternative
        // template <class InputIt, class OutputIt>
        // OutputIt get_rules_for_symbols(InputIt begin, InputIt end, OutputIt it) const
        // {
        //     for (InputIt init = begin; init != end; ++init) {
        //         for (rule::const_ptr r : rules_) {
        //             for (const std::vector<Symbolic_State<StateType> >& alternative : r->right()) {
        //                 if (alternative.size() == 2
        //                         && alternative[0] == init->first
        //                         && alternative[1] == init->second) {
        //                     *it++ = r;
        //                 }
        //             }
        //         }
        //     }
        //     return it;
        // }
    
        // Get the rules where this symbol occurs in a right hand side alternative
        template <class OutputIt>
        OutputIt get_rules_for_symbol(const Symbolic_State<StateType> & symbol, OutputIt it) const
        {
            for (typename rule<StateType>::const_ptr r : rules_) {
                for (const std::vector<Symbolic_State<StateType> >& alternative : r->right()) {
                    if (std::find(alternative.begin(), alternative.end(), symbol) != alternative.end()) {
                        *it++ = r;
                    }
                }
            }
            return it;
        }
    
        // Get the subject of the start symbol
        const Symbolic_State<StateType> & get_start_left() const
        {
            return rules_.front()->left();
        }
    
        // Is this symbol the left (subject) of a terminal?
        bool symbol_is_left_of_terminal(const Symbolic_State<StateType> & symbol) const
        {
            for (typename rule<StateType>::const_ptr r : rules_) {
                if (r->left() == symbol) {
                    for (const std::vector<Symbolic_State<StateType> >& alternative : r->right()) {
                        if (alternative.size() == 1 && symbol_is_terminal(alternative[0])) {
                            return true;
                        }
                    }
                }
            }
            return false;
        }
    
        const std::vector<Symbolic_State<StateType> > get_top_level_rules() const
        { return tops_; }

        // Pretty-print
        friend std::ostream& operator <<(std::ostream& os, const grammar& gram)
        {
            for (typename rule<StateType>::const_ptr r : gram.rules_) {
                os << *r;
                os << '\n';
            }
            return os;
        }

        void print_terminal() const
        {
            for(int i = 0; i < terminals_.size(); i++)
            {
                std::cerr << "(" << terminals_[i].GetContent() << ", " << terminals_[i].GetId() << ") ";
            }
        }

      private:
        std::vector<Symbolic_State<StateType>> terminals_;
        std::vector<typename rule<StateType>::ptr> rules_;
        std::vector<typename rule<StateType>::ptr> start_rules_;
        std::vector<Symbolic_State<StateType> > tops_;
    };

    template<class StateType>
    struct state
    {
        typename rule<StateType>::const_ptr rule_; // The grammar rule
        unsigned int right_; // Index of right hand side alternative
        unsigned int dot_; // Position of dot within symbols on right hand side
        unsigned int i_, j_; // Positions within the input
        char added_by_; // Which function added this state
        std::vector<std::pair<int, int> > back_pointer_; // back pointer to next level parsed rules: statelist #, position in statelist
        state(typename rule<StateType>::const_ptr rule, unsigned int right, unsigned int i, unsigned int j)
            : rule_(rule), right_(right), dot_(0), i_(i), j_(j), added_by_(0)
        {        }
    
        // Is dot all the way to the right?
        bool completed() const
        {
            return dot_ == rule_->right()[right_].size();
        }
    
        // Get symbol to the right of dot
        Symbolic_State<StateType> next_symbol() const
        {
            return rule_->right()[right_][dot_];
        }

        // void print_state_rules()
        // {
        //     std::vector<typename rule<StateType>::const_ptr> rules;
        //     grammar_.get_rules_for_left(st.next_symbol(), std::back_inserter(rules));
        //     for (typename rule<StateType>::const_ptr r : rules)
        //     {
        //         std::cout << "Left: ";
        //         std::cout << "(" << r->left().GetContent() << ", " << r->left().GetId() << ")\n";
        //         std::cout << "Right: \n";
        //         for (auto right : r->right())
        //         {
        //             for (auto iter = right.begin(); iter != right.end(); iter++)
        //             {
        //                 std::cout << "(" << (*iter).GetContent() << ", " << (*iter).GetId() << ") ";
        //             }
        //             std::cout << "\n";
        //         }
        //     }
        // }
    };

    // Pretty-print state
    template <class StateType>
    std::ostream& operator <<(std::ostream& os, const state<StateType>& st)
    {
        const std::vector<Symbolic_State<StateType> >& right = st.rule_->right()[st.right_];
        size_t rlen = right.size();
        os << '(';
        os << '('<<st.rule_->left().GetContent()<<", "<<st.rule_->left().GetId()<<")" << " -> ";
        unsigned int s;
        for (s = 0; s < rlen; ++s) {
            if (s == st.dot_) {
                os << "@ ";
            }
            os << "("<<right[s].GetContent()<<", "<<right[s].GetId()<<")";
            if (s < rlen - 1) {
                os << ' ';
            }
        }
        if (s == st.dot_) {
            os << " @";
        }
        os << ", [" << st.i_ << " , " << st.j_ << "]) ";
        switch (st.added_by_) {
            case 'P':
                os << "predictor";
                break;
            case 'S':
                os << "scanner";
                break;
            case 'C':
                os << "completer";
                break;
            default:
                os << "start state";
        }
        return os;
    }
        
    template <class StateType>
    bool operator == (const rule<StateType>& rule1, const rule<StateType>& rule2){
        if(rule1.left() == rule2.left() && rule1.right() == rule2.right()){
            return true;
        }
        else{
            return false;
        }
    }
        
    // Needed to check for duplicate states
    template <class StateType>
    bool operator ==(const state<StateType>& state1, const state<StateType>& state2)
    {
        if (*state1.rule_ == *state2.rule_
            && state1.right_ == state2.right_
            && state1.dot_ == state2.dot_
            && state1.i_ == state2.i_
            && state1.j_ == state2.j_
            && state1.back_pointer_== state2.back_pointer_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    
    // A statelist is a list of states
    template <class StateType>
    struct StateList
    {
        typedef std::vector<state<StateType> > statelist;
    };
    // A chart is a vector of statelists
    template <class StateType>
    struct chart
    {
        const grammar<StateType>& grammar_;
        std::vector<typename StateList<StateType>::statelist> chart_;
        std::vector<typename std::unordered_set<state<StateType> > > chart_set_;
        std::vector<unsigned int> succeeded_chart_idx;
    
        chart(const grammar<StateType> & grammar)
            : grammar_(grammar)
            // Chart begins with 1 statelist containing 1 dummy state used to track
            // successful parse
            // chart_(1, typename StateList<StateType>::statelist(1, state<StateType>(rule<StateType>::create(Symbolic_State<StateType>("$", false), grammar.get_start_left()),
            //                 0, 0, 0)))
        {
            typename StateList<StateType>::statelist top_list;
            std::unordered_set<state<StateType> > top_list_set;
            // std::cerr << "Initialize grammar initialization:before get_top_level_rules" << std::endl;            
            std::vector<Symbolic_State<StateType> > top_level_rules = grammar.get_top_level_rules();
            // std::cerr << "Initialize grammar initialization:after get_top_level_rules" << std::endl;                        
            for (Symbolic_State<StateType> source : top_level_rules)
            {
                state<StateType> source_state(rule<StateType>::create(Symbolic_State<StateType>("$", false), source), 0, 0, 0);
                top_list.push_back(source_state);
                top_list_set.insert(source_state);
            }
            // std::cerr << "Initialize grammar initialization:before push_back" << std::endl;
            chart_.push_back(top_list);
            chart_set_.push_back(top_list_set);
        }
    
        // Add state st to statelist s
        void add_state(const state<StateType> & st, unsigned int s)
        {
            // std::cerr << "add_state: chart_ size is: " << chart_.size() << std::endl;
            // std::cerr << "add_state: s is: " << s << std::endl;
            if (st.rule_->left().GetContent() == "$")
                succeeded_chart_idx.push_back(s);
            if (s < chart_.size()) {
                // Adding to the last statelist
                if(chart_set_[s].find(st) == chart_set_[s].end()){
                    chart_[s].push_back(st);
                    chart_set_[s].insert(st);
                }
            }
            else {
                // std::cerr << "add_state: inside else" << std::endl;                
                // Adding to a new statelist
                chart_.emplace_back(1, st);
                chart_set_.push_back(std::unordered_set<state<StateType> >({st}));
            }
            // std::cerr << "before return from add_state" << std::endl;
        }
    
        // Add predictions for the next symbol in this state
        void predictor(const state<StateType>& st)
        {   
            std::vector<typename rule<StateType>::const_ptr> rules;
            grammar_.get_rules_for_left(st.next_symbol(), std::back_inserter(rules));
            // std::cerr << "chart: for loop: predictor: outside for loop" << std::endl; 
            // std::cout << "For the next symbol: (" << st.next_symbol().GetContent() << ", " << st.next_symbol().GetId() << ")\n";
            // std::cout << "Print out all the rules: \n";
            // for (typename rule<StateType>::const_ptr r : rules)
            // {
            //     std::cout << "Left: ";
            //     std::cout << "(" << r->left().GetContent() << ", " << r->left().GetId() << ")\n";
            //     std::cout << "Right: \n";
            //     for(auto right : r->right())
            //     {
            //         for(auto iter = right.begin(); iter != right.end(); iter++)
            //         {
            //             std::cout << "(" << (*iter).GetContent() << ", " << (*iter).GetId() << ") ";
            //         }
            //         std::cout << "\n";
            //     }
            // }
            for (typename rule<StateType>::const_ptr r : rules) {
                // std::cerr << "chart: for loop: predictor: inside for loop" << std::endl; 
                // std::cerr << "rule left is: " << r->left().GetContent() << std::endl;                 
                for (unsigned int a = 0; a < r->right().size(); ++a) {
                    // std::cerr << "chart: for loop: predictor: inside for loop: inside for loop" << std::endl;
                    // std::cerr << "rule right size is: " << r->right()[0].size() << std::endl;
                    // std::cerr << "rule right is: ";
                    // for (int i = 0; i < r->right()[0].size(); i++)
                    // {
                    //     std::cerr << r->right()[0][i].GetContent() << " ";     
                    // }
                    // std::cerr << std::endl;  
                    // std::cerr << "chart: for loop: predictor: inside for loop: inside for loop: before prediction" << std::endl;                       
                    state<StateType> prediction = state<StateType>(r, a, st.j_, st.j_);
                    prediction.added_by_ = 'P';
                    // std::cerr << "chart: for loop: predictor: inside for loop: inside for loop: before add state" << std::endl;
                    add_state(prediction, st.j_);
                    // std::cerr << "chart: returned from add_state" << std::endl;
                    // std::cerr << "size of r->right().size(): " << r->right().size() << std::endl;
                }
            }
            // std::cerr << "before return from predictor" << std::endl;
        }
    
        // Scan input for next symbol
        void scanner(const state<StateType>& st, const std::vector<Symbolic_State<StateType> >& input, bool back_track = true)
        {
            // position of scanned word
            const Symbolic_State<StateType> & word = input[st.j_];
            // std::cout << "Print the scanned word: \n";
            // std::cout << "ID: " << word.GetId() << "Content: " << word.GetContent();
            // std::cout << "\nPrint out the pos of its statelist: " << st.j_ << "\n";
            // std::cerr << "chart: for loop: scanner: outside if" << std::endl;
            // std::cout << "The value of the variable right_: " << st.right_ << "\n";
            // std::cout << "The value of the variable dot_: " << st.dot_ << "\n";
            // std::cout << "Content in the right hand side: \n";
            // std::cout << "ID: \n";
            // for (auto it = st.rule_->right()[st.right_].begin(); it != st.rule_->right()[st.right_].end(); it++)
            // {
            //     std::cout << it->GetId() << " ";
            // }
            // std::cout << "\nContent: \n";
            // for (auto it = st.rule_->right()[st.right_].begin(); it != st.rule_->right()[st.right_].end(); it++)
            // {
            //     std::cout << it->GetContent() << " ";
            // }
            // std::cout << "\n";
            if (word == st.rule_->right()[st.right_][st.dot_]) {
                // std::cout << "Entered the if!\n";
                state<StateType> scanned = state<StateType>(st.rule_, st.right_, st.i_, st.j_ + 1);
                scanned.dot_ = st.dot_ + 1;
                scanned.added_by_ = 'S';
                if(back_track){
                    scanned.back_pointer_ = st.back_pointer_;
                }
                add_state(scanned, st.j_ + 1);
            }
        }
    
        // Complete states
        void completer(const state<StateType>& param_st, const int chart_i, const int chart_j, bool back_track = true)
        {
            state<StateType> st = param_st;
            Symbolic_State<StateType> left = st.rule_->left();
            std::vector<Symbolic_State<StateType> > right = st.rule_->right()[st.right_];
            const unsigned int i = st.i_;
            const unsigned int j = st.j_;
            
            unsigned kk = chart_[0][0].rule_->right().size();
            
            // std::cerr << "chart: for loop: completer: outside for loop, j is: " << j << ", i is: " << i <<  std::endl; 
            int k = 0;
            for (const typename StateList<StateType>::statelist& list = chart_[st.i_]; k < list.size(); k++) {
                // std::cerr << st.rule_->left().GetContent() << std::endl;
                // std::cerr << "chart: for loop: completer: for loop" << std::endl;
                if (list[k].j_ == i
                        && !list[k].completed()
                        && list[k].next_symbol() == st.rule_->left())
                {
                    Symbolic_State<StateType> left = list[k].rule_->left();
                    std::vector<Symbolic_State<StateType> > right = list[k].rule_->right()[st.right_];
                    // std::cerr << "chart: for loop: completer: for loop: before list" << std::endl;                                         
                    state<StateType> completed = state<StateType>(list[k].rule_, list[k].right_, list[k].i_, j);
                    // std::cerr << "chart: for loop: completer: for loop: after list" << std::endl;                                         
                    completed.dot_ = list[k].dot_ + 1;
                    completed.added_by_ = 'C';
                    if (back_track){
                        if(list[k].back_pointer_.size() != 0){
                            completed.back_pointer_ = list[k].back_pointer_; 
                        }               
                        completed.back_pointer_.push_back(std::pair<int, int>(chart_i, chart_j));
                    }
                    add_state(completed, j);

                    // std::cerr << "chart: for loop: completer: for loop: returned from add_state" << std::endl; 
                }
            }
            // std::cerr << "chart: for loop: returned from completer" << std::endl;             
        }
    
        // // We have succeeded if the completed dummy state is in the final statelist
        // bool succeeded() const
        // {
        //     const typename StateList<StateType>::statelist& list = chart_[chart_.size() - 1];
        //     return std::find_if(list.begin(), list.end(),
        //             [](const state<StateType> &st)->bool {
        //             return st.rule_->left() == "$" && st.completed(); })
        //         != list.end();
        // }
    
        int parsed_position() const
        {
            // chart 0 is dummy
            if(succeeded_chart_idx.empty())
                return 0;
            
            return succeeded_chart_idx.back();
        }

        // The main algorithm
        int parse(const std::vector<Symbolic_State<StateType> >& input, std::ostream *os,
                  std::vector<std::vector<typename StateList<StateType>::statelist> >& parsed_results, bool &parsing_success, bool back_track=true)
        {
            // std::cerr << "\n------ Inside Parse ------\n";
            parsing_success = false;
            // std::cerr << "The sequence to be parsed: \n";
            // for(auto it = input.begin(); it != input.end(); it++)
            // {
            //     std::cerr << "(" << (*it).GetContent() << ", " << (*it).GetId() << ") ";
            // }
            // std::cerr << "\n";

            // std::cerr << "Print terminal before parse: \n";
            
            // grammar_.print_terminal();
            
            for (unsigned int i = 0; i <= input.size(); ++i) {
                // std::cerr << "chart: for loop: chart_.size(): " << chart_.size() << ", i is: " << i << std::endl;
                if (chart_.size() > i) { // Check for running out of statelists when parse fails
                    // for (typename StateList<StateType>::statelist::iterator it = chart_[i].begin(); it != chart_[i].end(); ++it) {
                    typename StateList<StateType>::statelist::iterator it = chart_[i].begin();
                    int pos = 0;
                    while(it != chart_[i].end())
                    {
                        // std::cerr << "******inside iterator******" << std::endl;
                        const state<StateType> st = *it;
                        // if (!st.completed() && (std::find(and_or_node_ids.begin(), and_or_node_ids.end(), st.next_symbol().GetId()) != and_or_node_ids.end()))
                        // {
                        //     std::cout << "Print the content of the next symbol: " << "(" << st.next_symbol().GetContent() << ", " << st.next_symbol().GetId() << ")";
                            
                        //     std::cout << "Tracking the target and node:\n";
                        //     if (grammar_.symbol_is_terminal(st.next_symbol()))
                        //     {
                        //         std::cout << "Successfully added the and node to terminals_\n";                                
                        //     }
                        //     else
                        //         std::cout << "Probably failed!\n";
                        //     if (i < input.size())
                        //         scanner(st, input);
                        // }

                        if (!st.completed()
                                && !grammar_.symbol_is_terminal(st.next_symbol())) {
                            // std::cerr << "chart: for loop: called predictor" << std::endl;
                            // std::cout << "------ Inside Predictor -------" << std::endl;
                            predictor(st);
                            // std::cout << "------ Out of Predictor -------" << std::endl;
                            // std::cerr << "returned from predictor" << std::endl;
                        }
                        else if (!st.completed()) {
                            if (i < input.size()) {
                                // std::cout << "------ Inside Scanner -------" << std::endl;
                                // std::cerr << "Scanned sequence: " << std::endl;
                                // for(auto it = input.begin(); it != input.end(); it++)
                                // {
                                //     std::cerr << it->GetId() << " ";
                                // }
                                // std::cerr << "\n";
                                scanner(st, input, back_track);
                                // std::cout << "------ Out of Scanner -------" << std::endl;
                            }
                        }
                        else {
                            // std::cerr << "chart: for loop: called completer" << std::endl;                             
                            completer(st, i, pos, back_track);
                        }
                        pos++;
                        it = chart_[i].begin()+pos;
                    }
                    // if (os) {
                    //     *os << *this;
                    //     *os << '\n';
                    // }
                }
            }
            // return succeeded();
            // std::cerr << "Is parsing a success? " << succeeded << "\n";
            // std::cerr << "The parsed position is: " << parsed_position() << "\n";

            int parsed_pos = parsed_position();
            if (parsed_pos != 0 && !succeeded_chart_idx.empty())
            {
                parsing_success = true;
                std::vector<typename StateList<StateType>::statelist> parsed_chart(chart_.begin(),
                                                                    chart_.begin() + parsed_pos + 1);
                parsed_results.push_back(parsed_chart);
                return parsed_pos;
            }
            // std::cerr << "------ Outside Parse ------\n";

            // std::cerr << "chart:parsed position: " << parsed_position() << std::endl;
            return 0;
        }

        int parse2(const std::vector<Symbolic_State<StateType> >& input, std::ostream *os,
                  std::vector<std::vector<typename StateList<StateType>::statelist> >& parsed_results, bool &parsing_success, double &prob, std::shared_ptr<T_AOG<std::string> > graph_ptr)
        {
            // std::cerr << "\n------ Inside Parse ------\n";
            parsing_success = false;
            // std::cerr << "The sequence to be parsed: \n";
            // for(auto it = input.begin(); it != input.end(); it++)
            // {
            //     std::cerr << "(" << (*it).GetContent() << ", " << (*it).GetId() << ") ";
            // }
            // std::cerr << "\n";

            // std::cerr << "Print terminal before parse: \n";
            
            // grammar_.print_terminal();
            
            for (unsigned int i = 0; i <= input.size(); ++i) {
                // std::cerr << "chart: for loop: chart_.size(): " << chart_.size() << ", i is: " << i << std::endl;
                if (chart_.size() > i) { // Check for running out of statelists when parse fails
                    // for (typename StateList<StateType>::statelist::iterator it = chart_[i].begin(); it != chart_[i].end(); ++it) {
                    typename StateList<StateType>::statelist::iterator it = chart_[i].begin();
                    int pos = 0;
                    while(it != chart_[i].end())
                    {
                        // std::cerr << "******inside iterator******" << std::endl;
                        state<StateType>& st = *it;
                        // if (!st.completed() && (std::find(and_or_node_ids.begin(), and_or_node_ids.end(), st.next_symbol().GetId()) != and_or_node_ids.end()))
                        // {
                        //     std::cout << "Print the content of the next symbol: " << "(" << st.next_symbol().GetContent() << ", " << st.next_symbol().GetId() << ")";
                            
                        //     std::cout << "Tracking the target and node:\n";
                        //     if (grammar_.symbol_is_terminal(st.next_symbol()))
                        //     {
                        //         std::cout << "Successfully added the and node to terminals_\n";                                
                        //     }
                        //     else
                        //         std::cout << "Probably failed!\n";
                        //     if (i < input.size())
                        //         scanner(st, input);
                        // }

                        if (!st.completed()
                                && !grammar_.symbol_is_terminal(st.next_symbol())) {
                            // std::cerr << "chart: for loop: called predictor" << std::endl;
                            // std::cout << "------ Inside Predictor -------" << std::endl;
                            predictor(st);
                            // std::cout << "------ Out of Predictor -------" << std::endl;
                            // std::cerr << "returned from predictor" << std::endl;
                        }
                        else if (!st.completed()) {
                            if (i < input.size()) {
                                // std::cout << "------ Inside Scanner -------" << std::endl;
                                // std::cerr << "Scanned sequence: " << std::endl;
                                // for(auto it = input.begin(); it != input.end(); it++)
                                // {
                                //     std::cerr << it->GetId() << " ";
                                // }
                                // std::cerr << "\n";
                                scanner(st, input);
                                // std::cout << "------ Out of Scanner -------" << std::endl;
                            }
                        }
                        else {
                            // std::cerr << "chart: for loop: called completer" << std::endl;                             
                            completer(st, i, pos);
                        }
                        pos++;
                        it = chart_[i].begin()+pos;
                    }
                    // if (os) {
                    //     *os << *this;
                    //     *os << '\n';
                    // }
                }
            }
            // return succeeded();
            // std::cerr << "Is parsing a success? " << succeeded << "\n";
            // std::cerr << "The parsed position is: " << parsed_position() << "\n";

            int parsed_pos = parsed_position();
            if (parsed_pos != 0 && !succeeded_chart_idx.empty())
            {
                parsing_success = true;
                std::vector<typename StateList<StateType>::statelist> parsed_chart(chart_.begin(),
                                                                    chart_.begin() + parsed_pos + 1);
                parsed_results.push_back(parsed_chart);                
                
                double total_prob = 0;
                for (const std::vector<typename StateList<StateType>::statelist> &chart : parsed_results)
                {
                    // std::cout << "~" << std::endl;
                    // best_parsing_in_chart: best parsing in each run of the parser
                    // Symbolic_State<StateType> best_parsing_start;
                    int last_statelist_index = chart.size() - 1;
                    // find the top most root in the last statelist (i.e. all possible parsing for this seq)
                    for (state<StateType> st : chart[last_statelist_index])
                    {
                        
                        // find top most root
                        if (st.i_ == 0 && st.j_ != 0 && st.rule_->left().GetContent() == "$" && st.completed())
                        {
                            // std::cout << "haha" << std::endl;

                            /*std::cerr << "The state is: [" << st.i_ << "," << st.j_ << "], FROM " << st.rule_->left().GetContent() << " TO ";
                            auto results =  st.rule_->right()[0];
                            for (auto result : results)
                            {
                                std::cerr << result.GetContent() << " ";
                            }
                            std::cerr << std::endl; */
                            // @param prob: probability of this particular possible parsing
                            double prob_ = 1;
                            
                            // @param q: queue used to backtrack all subparsing in one possible parsing
                            std::queue<state<StateType> > q;
                            // @param possible_parsing: possible parsing in one chart
                            Symbolic_State<StateType> possible_parsing_start;
                            q.push(st);
                            while (!q.empty())
                            {
                                st = q.front();
                                q.pop();
                                // this rule must be fully parsed to be considered as a valid parse
                                /** std::cerr << "back pointer size is: " << st.back_pointer_.size() << std::endl; */
                                for (std::pair<int, int> back_pointer : st.back_pointer_)
                                {
                                    auto pointed_state = chart[back_pointer.first][back_pointer.second];
                                    /** std::cerr << "Back Pointer points to: [" << pointed_state.i_ << "," << pointed_state.j_ << "], FROM " << pointed_state.rule_->left().GetContent() << " TO ";
                                    auto results =  pointed_state.rule_->right()[0];
                                    for (auto result : results)
                                    {
                                        std::cerr << result.GetContent() << " ";
                                    }
                                    std::cerr << std::endl; */
                                    q.push(chart[back_pointer.first][back_pointer.second]);
                                }

                                // skip the dummy rule from $ -> rest
                                if (st.rule_->left().GetContent() == "$")
                                {
                                    assert(st.rule_->right()[st.right_].size() == 1);
                                    possible_parsing_start = st.rule_->right()[st.right_][0];
                                    continue;
                                }
                                /** std::cerr << "parent state is:" << st.rule_->left().GetContent() << std::endl; */
                                VertexId parentId = graph_ptr->GetVertexIdByState(st.rule_->left());
                                // try
                                // {
                                //     // auto all_states = graph_ptr->GetStates();
                                //     /** std::cerr<<"size of states in AOG: "<<all_states.size()<<std::endl; */
                                //     parentId = graph_ptr->GetVertexIdByState(st.rule_->left());
                                    
                                // }
                                // catch(std::exception e)
                                // {
                                    
                                //     // for(auto state : all_states)
                                //     //     std::cerr<<state.GetId()<<std::endl;
                                    
                                //     // AOFStruct<StateType> peeker;
                                //     /** std::cerr << "AOF id_: " << peeker.id_; */
                                //     // auto all_states = graph_ptr->GetStates();
                                //     /** std::cerr<<"size of states in AOG: "<<all_states.size()<<std::endl;
                                //     for(auto state : all_states)
                                //     {
                                //         std::cerr<<"State Content: "<<state.GetContent()<<"_\n"<<"State ID: "<<state.GetId()<<std::endl;
                                //     }

                                //     std::cout << e.what() << std::endl; */
                                // }
                                // std::cout << "xixi" << std::endl;

                                // if it is Or-node, find the corresponding weight and update likelihood of this possible parsing
                                if (!graph_ptr->GetVertexContent(parentId)->IsAnd())
                                {
                                    /** std::cerr << "The rule's source is an Or-node" << std::endl; */
                                    // find all out edge weights of the source
                                    std::unordered_map<VertexId, double> outEdgeWeights =
                                            graph_ptr->GetOutEdgeWeights(parentId, true);

                                    // locate the outedge we are looking for
                                    // find all right-hand-side vertexids
                                    // @param right_hand_state_ids: the ids of right hand side states (e.g. S->NP VP. ids of NP, VP)
                                    std::vector<VertexId> right_hand_state_ids;
                                    for (Symbolic_State<StateType> right_hand_state : st.rule_->right()[st.right_])
                                        right_hand_state_ids.push_back(
                                                graph_ptr->GetVertexIdByState(right_hand_state));
                                    // find the children that has all the right-hand-side vertex,
                                    // i.e. corresponds to the right-hand-rule
                                    // @param dummy_vertices: the ids of the dummy vertices under the Or-node of parent (e.g. S->NP VP. dummy nodes under S)
                                    std::vector<VertexId> dummy_vertices =
                                            graph_ptr->ChildrenVertices(parentId);

                                    for (VertexId id : dummy_vertices)
                                    {
                                        // @param found_destination: flag indicating the (e.g. S->NP VP) branch is found
                                        // @param children_vertices: children vertices of the dummy node we are looking at
                                        bool found_destination = true;
                                        std::vector<VertexId> children_vertices =
                                                graph_ptr->ChildrenVertices(id);

                                        // if the children of the dummy node is of different size than our target, then this is not the branch we are looking for
                                        if (children_vertices.size() != right_hand_state_ids.size())
                                            continue ;
                                        // if size is equal, check all children are same
                                        for (int i = 0; i < children_vertices.size(); i++)
                                        {
                                            if (children_vertices[i] != right_hand_state_ids[i])
                                            {
                                                found_destination = false;
                                                break;
                                            }
                                        }
                                        if (found_destination)
                                        {
                                            double weight = outEdgeWeights[id];
                                            prob_ *= weight;
                                            // std::cout << "weight is: " << weight << ", probability is: " << prob_ << std::endl;
                                            break;
                                        }
                                    }
                                } 
                            }


                            // else find the largest probability parsing
                            // if (prob_ > best_prob)
                            // {
                            //     best_prob = prob_;
                            //     best_parsing_start = possible_parsing_start;
                            // }
                            total_prob += prob_;
                        }
                        // else{
                        //     std::cerr << "The state is: [" << st.i_ << "," << st.j_ << "], FROM " << st.rule_->left().GetContent() << " TO ";
                        // }
                    }
                }
                prob = total_prob;
                return parsed_pos;
            }
            // std::cerr << "------ Outside Parse ------\n";

            // std::cerr << "chart:parsed position: " << parsed_position() << std::endl;
            return 0;
        }
    
        // Pretty-print chart
        friend std::ostream& operator <<(std::ostream& os, const chart &ch)
        {
            for (unsigned int i = 0; i < ch.chart_.size(); ++i) {
                os << "S" << i << ": ";
                os << '[';
                unsigned int j = 0;
                for (const state<StateType>& st : ch.chart_[i]) {
                    os << st;
                    if (j < ch.chart_[i].size() - 1) {
                        os << ",\n";
                    }
                    ++j;
                }
                os << ']';
                os << "\n\n";
            }
            return os;
        }
    };

    template <class StateType>
    class EarleyParser
    {
    public:
        EarleyParser(){}
        EarleyParser(const grammar<StateType>& grammar)
            : grammar_(grammar)
        {
            // for (auto it = grammar_.rules().begin(); it != grammar_.rules().end(); it++){
            //     std::cerr << "From " << (*it)->left().GetContent() << std::endl;
            //     for (int i = 0; i < (*it)->right()[0].size(); i++)
            //         std::cerr << "TO " << (*it)->right()[0][i].GetContent() << std::endl;
            // }      
            // std::cerr << "EarleyParser initialization: private grammar size is: " << grammar_.rules().size() << std::endl;
        }
        template <class InputIt>
        int parse(InputIt begin, InputIt end)
        {
            std::vector<Symbolic_State<StateType> > input;
            std::copy(begin, end, std::back_inserter(input));
            return chart<StateType>(grammar_).parse(input, nullptr, parsed_results);
        }
        template <class InputIt>
        int parse(InputIt begin, InputIt end, std::ostream& os, bool & parsing_success, bool back_track = true)
        {
            // std::cerr << "EarleyParser: grammar size is: " << grammar_.rules().size() << std::endl;
            // for (auto it = grammar_.rules().begin(); it != grammar_.rules().end(); it++){
            //     std::cerr << "From " << (*it)->left().GetContent() << std::endl;
            //     for (int i = 0; i < (*it)->right()[0].size(); i++)
            //         std::cerr << "TO " << (*it)->right()[0][i].GetContent() << std::endl;
            // }
            std::vector<Symbolic_State<StateType> > input;
            std::copy(begin, end, std::back_inserter(input));
            // std::cerr << "EarleyParser: before chart construction " << std::endl;
            int pos = chart<StateType>(grammar_).parse(input, &os, parsed_results, parsing_success, back_track);
            // std::cerr << "EarleyParser: position is: " << pos << std::endl;
            return pos;
        }
        template <class InputIt>
        int parse2(InputIt begin, InputIt end, std::ostream& os, bool & parsing_success, double &prob, std::shared_ptr<T_AOG<std::string> > graph_ptr)
        {
            // std::cerr << "EarleyParser: grammar size is: " << grammar_.rules().size() << std::endl;
            // for (auto it = grammar_.rules().begin(); it != grammar_.rules().end(); it++){
            //     std::cerr << "From " << (*it)->left().GetContent() << std::endl;
            //     for (int i = 0; i < (*it)->right()[0].size(); i++)
            //         std::cerr << "TO " << (*it)->right()[0][i].GetContent() << std::endl;
            // }
            std::vector<Symbolic_State<StateType> > input;
            std::copy(begin, end, std::back_inserter(input));
            // std::cerr << "EarleyParser: before chart construction " << std::endl;
            int pos = chart<StateType>(grammar_).parse2(input, &os, parsed_results, parsing_success, prob, graph_ptr);
            // std::cerr << "EarleyParser: position is: " << pos << std::endl;
            return pos;
        }

        double prob(std::shared_ptr<T_AOG<std::string> > graph_ptr){
            double total_prob = 0;
            for (const std::vector<typename StateList<StateType>::statelist> &chart : this->parsed_results)
            {
                int last_statelist_index = chart.size() - 1;
                for (const state<StateType>& st_last : chart[last_statelist_index])
                {
                    // find top most root
                    if (st_last.i_ == 0 && st_last.j_ != 0 && st_last.rule_->left().GetContent() == "$" && st_last.completed())
                    {
                        double prob_ = 1;
                        
                        // @param q: queue used to backtrack all subparsing in one possible parsing
                        std::queue<state<StateType> > q;
                        // @param possible_parsing: possible parsing in one chart
                        Symbolic_State<StateType> possible_parsing_start;
                        q.push(st_last);

                        while (!q.empty())
                        {
                            state<StateType> st = q.front();
                            q.pop();
                            for (const std::pair<int, int>& back_pointer : st.back_pointer_)
                            {
                                auto pointed_state = chart[back_pointer.first][back_pointer.second];
                                q.push(chart[back_pointer.first][back_pointer.second]);
                            }

                            if (st.rule_->left().GetContent() == "$")
                            {
                                assert(st.rule_->right()[st.right_].size() == 1);
                                possible_parsing_start = st.rule_->right()[st.right_][0];
                                continue;
                            }

                            VertexId parentId = graph_ptr->GetVertexIdByState(st.rule_->left());

                            // if it is Or-node, find the corresponding weight and update likelihood of this possible parsing
                            if (!graph_ptr->GetVertexContent(parentId)->IsAnd())
                            {
                                std::unordered_map<VertexId, double> outEdgeWeights = graph_ptr->GetOutEdgeWeights(parentId, true);

                                // locate the outedge we are looking for
                                // find all right-hand-side vertexids
                                // @param right_hand_state_ids: the ids of right hand side states (e.g. S->NP VP. ids of NP, VP)
                                std::vector<VertexId> right_hand_state_ids;
                                for (const Symbolic_State<StateType>& right_hand_state : st.rule_->right()[st.right_])
                                    right_hand_state_ids.push_back(graph_ptr->GetVertexIdByState(right_hand_state));

                                // find the children that has all the right-hand-side vertex,
                                // i.e. corresponds to the right-hand-rule
                                // @param dummy_vertices: the ids of the dummy vertices under the Or-node of parent (e.g. S->NP VP. dummy nodes under S)
                                std::vector<VertexId> dummy_vertices = graph_ptr->ChildrenVertices(parentId);

                                for (VertexId id : dummy_vertices)
                                {
                                    // @param found_destination: flag indicating the (e.g. S->NP VP) branch is found
                                    // @param children_vertices: children vertices of the dummy node we are looking at
                                    std::vector<VertexId> children_vertices =
                                            graph_ptr->ChildrenVertices(id);

                                    // if the children of the dummy node is of different size than our target, then this is not the branch we are looking for
                                    if (children_vertices != right_hand_state_ids){
                                        prob_ *= outEdgeWeights[id];
                                        break;
                                    }
                                }
                            } 
                        }

                        total_prob += prob_;
                    }
                }
            }
            return total_prob;
        }

        std::vector<std::vector<typename StateList<StateType>::statelist > > GetPartialParse(){
            return parsed_results;
        }
        grammar<StateType> GetGrammar(){
            return grammar_;
        }

        void update(const std::vector<Symbolic_Rule<StateType> > & rules, const SequenceType<StateType> &input_seq)
        {
            SequenceType<StateType> top_level_rules = this->get_top_level_rules(rules);
            this->grammar_.update(rules, top_level_rules, input_seq);
        }

        std::vector<Symbolic_State<StateType> > get_top_level_rules(std::vector<Symbolic_Rule<StateType> > rules)
        {
            std::vector<Symbolic_State<StateType> > top_level_rules;
            std::vector<Symbolic_State<StateType> > sources;
            std::unordered_set<Symbolic_State<StateType> > results;
            for (Symbolic_Rule<StateType> rule : rules)
            {
                sources.push_back(rule.GetSource());
                results.insert(rule.GetResults().begin(), rule.GetResults().end());
            }
            
            // check which source is not other sources' result
            for (Symbolic_State<StateType> source : sources)
            {
                if (results.find(source) == results.end())
                    top_level_rules.push_back(source);
            }
            
            return top_level_rules;
        }
    private:
        grammar<StateType> grammar_;
        std::vector<std::vector<typename StateList<StateType>::statelist > > parsed_results;
    };

}

namespace std
{
    template <typename StateType>
    size_t hash_value(const AOG_LIB::rule<StateType> &rule){
        std::size_t seed = 0;
        boost::hash_combine(seed, rule.left());
        boost::hash_combine(seed, rule.right());
        return seed;
    }

    template <typename StateType>
    struct hash<AOG_LIB::rule<StateType> >
    {
        size_t operator()(const AOG_LIB::rule<StateType> &rule) const noexcept
        {
            boost::hash<AOG_LIB::rule<StateType> > hasher;
            return hasher(rule);
        }
    };

    template <typename StateType>
    struct hash<AOG_LIB::state<StateType> >
    {
        size_t operator()(const AOG_LIB::state<StateType> &state) const noexcept
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, *state.rule_);
            boost::hash_combine(seed, state.right_);
            boost::hash_combine(seed, state.dot_);
            boost::hash_combine(seed, state.i_);
            boost::hash_combine(seed, state.j_);
            boost::hash_combine(seed, state.back_pointer_);
            return seed;
        }
    };
}


#endif //AOG_LIB_EARLEY_H
