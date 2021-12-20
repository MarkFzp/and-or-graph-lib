#include<iostream>
#include <string>
#include "Metrics.h"
using namespace std;
using namespace AOG_LIB;

//argb[1]: function name
//argv[2]: the learned grammar's path
//argv[3]: the groundtruth grammar's path / test set's path
int main(int argc, char* argv[])
{
    string option = argv[1];
    if(option == "perplexity" && argc == 5){
        double perplexity = Perplexity_cross_entropy(argv[2],argv[3], stod(argv[4]));
        cout << "perplexity: " << perplexity << endl;
    }
    else if(option == "precision" && argc == 5){
        double precision = Precision(argv[2],argv[3],stoi(argv[4]));
        cout << "precision: " << precision << endl;
    }
    else if(option == "recall" && argc == 4){
        double recall = Recall(argv[2],argv[3]);
        cout << "recall: " << recall << endl;
    }
    else if(option == "sample" && argc == 5){
        Sample_Viterbi(argv[2], stoul(argv[3]), (bool)stoi(argv[4]));
    }
    else if(option == "sample_total" && argc == 4){
        AOG_LIB_UTIL::Sample_Total(argv[2], stoul(argv[3]));
    }
    else if(option == "kl" && argc == 4){
        double kl = KL_Divergence_Ground_Truth_Tree(argv[2], argv[3]);
        cout << "kl divergence: " << kl << endl;
    }
    else if(option == "kl2" && argc == 4){
        double kl = KL_Divergence_Train_Set(argv[2], argv[3]);
        cout << "kl divergence: " << kl << endl;
    }
    else if(option == "js" && argc == 5){
        double js = JS_Divergence_By_Sampling(argv[2], argv[3], stoul(argv[4]));
        cout << "js divergence: " << js << endl;
    }
    else if(option == "bleu" && argc == 4){
        double bleu = Bleu(argv[2], argv[3]);
        cout << "bleu: " << bleu << endl;
    }
    else if(option == "parse_tree" && argc == 4)
    {
        GetBestParsingTree(argv[2],argv[3]);
    }
    else{
        cerr << "format: ./executable [perplexity/precision/recall] [file1(learned_tree file)] [file2(test set file/learned_tree of ground truth) / # of samples] [perplexity penalty / precision sample times]" << endl;
    }

}
