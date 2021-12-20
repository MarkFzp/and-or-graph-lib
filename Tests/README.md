## Requirements
* C++ Boost Libraries

## To compile the test, run:
`g++ test_run.cpp -g -std=c++11 -lboost_system -lboost_filesystem -o test_run`

## Run the test:
To run the test, please make sure that:
1. You created a folder in Share_Results(Dropbox folder)
2. You have the config file inside that folder

The test takes in two arguments:
1. Absolute path to the folder that you created in Share_Results
2. file name of the text file you want to learn(these files are under Share_Results/Grammar_Example/Output)

### For example:
`./test_run ~/Dropbox/Share_Results/GRAMMAR5_ONLINE_TEST/ GRAMMAR_EXAMPLE5_OUTPUT.txt`


# To run online learning, modify START_EXPAND in the config.txt


