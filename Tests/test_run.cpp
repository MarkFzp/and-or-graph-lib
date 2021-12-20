//
// Created by Luyao Yuan on 18/4/8.
// Compile this file with options -lboost_system -lboost_filesystem

#include <vector>
#include <string>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include <chrono>
#include "../T-AOG/T_AOG.h"
#include "../Learner/Online_Learner.h"
// #include "../Utils/Metrics.h"
#define fs boost::filesystem

using namespace std;
using namespace AOG_LIB;

// extern std::vector<int> NUM_OF_BIGRAM;
// extern std::vector<double> NUM_OF_VARIATION;
// extern std :: vector<double> BEST_POSTERIOR;



bool copyFile(const char *SRC, const char *DEST)
{
	std::ifstream src(SRC, std::ios::binary);
	std::ofstream dest(DEST, std::ios::binary);
	dest << src.rdbuf();
	return src && dest;
}

int main(int argc, char *argv[])
{
	
	// check if directory exists
	if (argc > 2)
	{
		if (!boost::filesystem::exists(std::string(argv[1])))
		{
			std::cout << "Directory does not exist:" << argv[1];
			throw std::exception();
		}
		else
		{
			std::string config_path = std::string(argv[1]) + "/config.txt";
			if (!boost::filesystem::exists(config_path))
			{
				std::cout << "Config file need to be in the directory!";
				throw std::exception();
			}
			else
			{
				PATH = std::string(argv[1]);
			}
		}
	}
	else
	{
		std::cout << "Argument needed!" << std::endl;
		throw std::exception();
	}

	string FILEPATH = argv[2];
	// FILEPATH = "../../Grammar_Example/Output/"+ FILEPATH;
	// FILEPATH = "../../Grammar_Example/Train_test_split/"+FILEPATH;
	// FILEPATH = "../../Unique_Grammar_Example/Train_test_split/"+FILEPATH;

	cout << PATH + FILEPATH << endl;
	string filename(PATH + FILEPATH);
	vector<SequenceType<string>> dataset = FileParser<string>(filename);

	std::string str = PATH + "all_posteriors.txt";
	const char* file_name = str.c_str();
	int ret_code = std::remove(file_name);
	if (ret_code == 0) {
		std::cout << "File was successfully deleted\n";
	} else {
		std::cerr << "Error during the deletion: " << ret_code << '\n';
	}

	cout << "Read in dataset with " << dataset.size() << " data sequences\n";
	cout << "size of first data" << dataset[0].size() << "\n";
	shared_ptr<T_AOG<string>> graph_ptr(new T_AOG<string>());
	//shared_ptr<T_AOG<string>> graph_ptr = ReconstructTree(path + "learned_tree.txt");
	AOG_LIB::Online_Learner<string> learner = AOG_LIB::Online_Learner<string>(graph_ptr, PATH + "/config.txt");
	
	auto start = std::chrono::system_clock::now();

	DATA_SET_SIZE = dataset.size();

	for (int i = 0; i < DATA_SET_SIZE; ++i)
	{
		cout << "[Round] " << i << endl;
		learner.Learn(dataset[i]);
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> learning_time = end-start;
	std::cout << "The learning time: "<<learning_time.count()<<std::endl;

	if(OFFLINE){
		std::shared_ptr<T_AOG<std::string> > best_graph_ptr = learner.GetGraphPtr();
		best_graph_ptr->OutputLearnedTree(PATH);
	}

	// ofstream stats_file;
	// stats_file.open("../Statistics/stats.txt", ofstream::out | ofstream::trunc);
	// if(stats_file.is_open())
	// {
	// 	for(auto num_of_bigram : NUM_OF_BIGRAM)
	// 		stats_file << num_of_bigram<<" ";
	// 	stats_file << std::endl;
		
	// 	for(auto num_of_variation : NUM_OF_VARIATION)
	// 		stats_file << num_of_variation << " ";
	// 	stats_file << std::endl;

	// 	for(auto best_posterior : BEST_POSTERIOR)
	// 		stats_file << best_posterior << " ";
	// 	stats_file << std::endl;
		
	// 	for(auto kl_div : KL_DIVERGENCE)
	// 		stats_file << kl_div << " ";
	// 	stats_file << std::endl;
		
	// 	stats_file << TOTAL_VARIATIONS<<std::endl;
	// }

	// stats_file.close();
	
	return 0;
}
