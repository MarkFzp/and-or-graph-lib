#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <unordered_map>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        cerr << "Wrong usage, ./cpp learned_file groundtruth_file\n";
        exit(1);
    }
    string filename1 = string(argv[1]);
    string filename2 = string(argv[2]);
    ifstream f1(filename1, std::ifstream::in);
    ifstream f2(filename2, std::ifstream::in);

    unordered_map<string, int> f1_data;
    string s = "";
    int s1_size = 0, s2_size = 0;
    double kl = 0.0;
    while (getline(f1, s))
    {
        if (f1_data.find(s) == f1_data.end())
            f1_data.insert(make_pair(s, 0));
        f1_data[s]++;
        s1_size++;
        s = "";
    }
    unordered_map<string, int> f2_data;
    while (getline(f2, s))
    {
        if (f2_data.find(s) == f2_data.end())
            f2_data.insert(make_pair(s,0));
        f2_data[s]++;
        s2_size++;
        s = "";
    }    

    // cout << "-------------------------------\n";
    // cout << "Data in learned tree\n";
    for (auto it = f1_data.begin(); it != f1_data.end(); it++)
    {
        // cerr << it->first << "," << it->second;
        auto it_2 = f2_data.find(it->first);
        // cerr << "," << it->second - it_2->second << endl;
        if ( it_2 == f2_data.end())
        {
            cerr << "Error, Cannot find learned data in ground truth\n";
            cerr << "the data is: \n"<< it->first<<std::endl;

            exit(1);
        }
        // cerr << "Find ground truth, the data is: \n"<< it->first<<std::endl;
        
        double px = 1.0 * it->second/s1_size;
        double qx = 1.0 * it_2->second/s2_size;
        kl += (1.0 * it->second/s1_size) * log(px/qx);
    }

    // cout << "-------------------------------\n";
    // cout << "Data in ground truth\n";
    
    // for (auto it = f2_data.begin(); it != f2_data.end(); it++)
    //     cerr << it->first << "," << it->second << endl; 



    // cout << "-------------------------------\n";
    // cout << "Data in ground truth but not in learned tree\n";
    

    for (auto it = f2_data.begin(); it != f2_data.end(); it++)
    {
        auto it_1 = f1_data.find(it->first);
        if (it_1 == f1_data.end())
        {
            // cout << it->first << endl;
        }
    }
    cout << kl << endl;


    // ofstream output_file;
    // output_file.open("./KL_stats.txt",ofstream :: out | ofstream:: trunc);
    // if(output_file.is_open())
    // {
    //     output_file << kl << " ";
    // }
    f1.close();
    f2.close();
    // output_file.close();
}