all: minimize map search classify
	
minimize:
	g++ minimize.cpp -std=c++11 -O3 -o minimize

map:
	g++ consult_map.cpp -std=c++11 -O3 -o consult_map

search:
	g++ consult_search.cpp -Wno-unused-result -std=c++11 -fopenmp -O3 -o consult_search

classify:
	g++ consult_classify.cpp -std=c++11 -fopenmp -O3 -o consult_classify

clean:
	rm -f ./minimize ./consult_map ./consult_search ./consult_classify
	@echo "Succesfully cleaned."
