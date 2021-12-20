#include <iostream>
#include<fstream>
#include<cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <ctime>
#include <algorithm>
#include <iterator>
#include <boost/functional/hash.hpp>

using namespace std;


template <typename Container> // (hasher for containers) we make this generic for any container
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};

double Bleu (string reference_file, string hypothesis_file){
    istream_iterator<string> eos;
    istreambuf_iterator<char> eoc;
    
    // clock_t start;

    // ifstream reference_stream_tmp(reference_file);
    ifstream reference_stream(reference_file);
    // size_t ref_line_count = count(istreambuf_iterator<char>(reference_stream_tmp), eoc, '\n') + 1;
    
    vector<vector<string>> reference;
    // reference.reserve(ref_line_count);
    unordered_set<string> ref_set;
    // ref_set.reserve(ref_line_count);
    
    if(reference_stream.is_open()){
        string line;
        while(getline(reference_stream, line)){
            if(ref_set.find(line) == ref_set.end()){
                ref_set.insert(line);
                istringstream iss(line);
                reference.emplace_back(istream_iterator<string>(iss), eos);
            }
        }
    }
    else{
        cerr << "data set not found for reference" << endl;
        throw exception();
    }
    // if(reference.empty()){
    //     cerr << "empty reference" << endl;
    //     throw exception();
    // }


    // ifstream hypothesis_stream_tmp(hypothesis_file);
    ifstream hypothesis_stream(hypothesis_file);
    // size_t hyp_line_count = count(istreambuf_iterator<char>(hypothesis_stream_tmp), eoc, '\n') + 1;
    
    vector<vector<string>> hypothesis;
    // hypothesis.reserve(hyp_line_count);
    vector<double> hyp_likelihood;
    // hyp_likelihood.reserve(hyp_line_count);
    
    if(hypothesis_stream.is_open()){
        string line;
        while(getline(hypothesis_stream, line)){
            size_t index;
            double likelihood = stod(line, &index);
            hyp_likelihood.push_back(likelihood);
            istringstream iss(line.substr(index + 2));
            hypothesis.emplace_back(istream_iterator<string>(iss), eos);
        }
    }
    else{
        cerr << "sample data file not found" << endl;
        throw exception();
    }
    // if(hypothesis.empty() || hyp_likelihood.empty()){
    //     cerr << "empty hypothesis or hypothesis likelihood" << endl;
    //     throw exception();
    // }

    // double timenow = (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000);

    double total_likelihood = 0.0;
    for_each(hyp_likelihood.begin(), hyp_likelihood.end(), [&total_likelihood](double tmp){total_likelihood += tmp;});


    double bleu = 0.0;
    long hyp_index = -1;
    for(const auto& hyp : hypothesis){
        ++hyp_index;

        // for brevity penalty, refering to brevity_penalty in https://www.nltk.org/_modules/nltk/translate/bleu_score.html
        size_t hyp_len = hyp.size();

        pair<unsigned, size_t> absDis_refSize_pair(abs((signed)reference[0].size() - (signed)hyp_len), reference[0].size());
        for(auto& it = ++reference.begin(); it != reference.end(); ++it){
            size_t ref_size = (*it).size();
            unsigned diff = abs((signed)ref_size - (signed)hyp_len);
            if(diff <= absDis_refSize_pair.first && ref_size < absDis_refSize_pair.second){
                absDis_refSize_pair.first = diff;
                absDis_refSize_pair.second = ref_size;
            }
        }

        double brevity_penalty;
        unsigned closed_ref_length = absDis_refSize_pair.second;
        if (hyp_len > closed_ref_length){
            brevity_penalty = 1.0;
        }
        else{
            brevity_penalty = exp(1 -  (double)closed_ref_length / hyp_len);
        }


        size_t len = hyp_len < 4 ? hyp_len : 4;
        
        double hyp_score = 0.0;
        //c++11: If there is no key equivalent to x in the map, inserts value_type(x, T()) into the map.

        for(size_t n_gram = 1; n_gram <= len; n_gram++){
            unordered_map<vector<string>, size_t, container_hash<vector<string> > > counts;
            for(size_t i = 0; i <= hyp_len - n_gram; i++){
                vector<string> word(n_gram);
                for(size_t j = 0; j < n_gram; j++){
                    word[j] = hyp[i + j];
                }
                counts[word]++;
            }
            unordered_map<vector<string> , size_t, container_hash<vector<string> > > clipped_counts;
            for(const auto& ref : reference){
                for(const auto& count : counts){
                    size_t current_clipped_count = 0;
                    auto it = ref.begin();
                    for(int i = 0; i <= (signed)ref.size() - (signed)n_gram; i++){ // must use signed type, ref.size() - n_gram can be -ve
                        if (equal(count.first.begin(), count.first.end(), it))
                            current_clipped_count++;
                        it++;
                    }
                    if(current_clipped_count > clipped_counts[count.first]){
                        clipped_counts[count.first] = current_clipped_count;
                    }
                }
            }

            size_t numerator = 0;
            size_t denominator = 0;
            
            for(const auto& count : counts){
                denominator += count.second;
                unsigned clipped_count = clipped_counts[count.first]; 
                if(clipped_count > count.second)
                    numerator += count.second;
                else
                    numerator += clipped_count;
            }

            if(denominator == 0){
                denominator = 1; // to avoid devision by 0 error
            }

            // if(numerator == 0){
            //     hyp_score = 0.0;
            //     break;
            // }
            // else{
                hyp_score += (1.0/len) * log(((double)numerator / denominator));
            // }
        }
        
        // if(hyp_score != 0.0){
            hyp_score = brevity_penalty * exp(hyp_score);
        // }
        bleu += hyp_likelihood[hyp_index] * hyp_score;
    }

    reference_stream.close();
    hypothesis_stream.close();

    bleu /= total_likelihood;

    return bleu;
}

int main(int argc, char* argv[]){

    if (argc == 3){
        double bleu = Bleu(argv[1], argv[2]);
        cout << "bleu: " << bleu << endl;
    }
    else{
        cerr << "wrong argument format" << endl;
        exit(1);
    }
}