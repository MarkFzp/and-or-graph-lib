//
// Created by Hao Wang on 18/4/16.
//

#include <vector>
#include <string>

#include "../T-AOG/T_AOG.h"
#include "../Learner/Online_Learner.h"

using namespace std;
using namespace AOG_LIB;

int main()
{
	vector<SequenceType<string> > dataset = FileParser<string>("../Grammar_Example/g1_data.txt");
	cout << "Read in dataset with " << dataset.size() << " data sequences\n";
	cout<<"size of first data"<<dataset[0].size() << "\n";
	shared_ptr<T_AOG<string> > graph_ptr(new T_AOG<string>());
	AOG_LIB::Online_Learner<string> learner = AOG_LIB::Online_Learner<string>(graph_ptr);

	for(auto d: dataset)
	{
		learner.Learn(d);
	}

	return 0;
}