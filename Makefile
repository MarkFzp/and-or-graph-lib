run:
	g++ -O3 -std=c++11 Tests/test_run.cpp -o test_run -lboost_system -lboost_filesystem

eval:
	g++ -O3 -fuse-ld=gold -fsanitize=undefined -std=c++11 Evaluation/evaluation.cpp -o Evaluation/test_eval -lboost_system -lboost_filesystem

construct:
	g++ -O3 -std=c++11 Utils/construct_tree_from_benchmark.cpp -o Utils/test_construct -lboost_system -lboost_filesystem

clean:
	if [ -f test_run ] ; \
	then \
		rm test_run ; \
	fi;
