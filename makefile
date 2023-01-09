all: minimize map search
	
minimize:
	g++ minimize.cpp -std=c++11 -o minimize

map:
	g++ consult_map.cpp -std=c++11 -O3 -o consult_map

search:
	g++ consult_search.cpp -std=c++11 -fopenmp -O3 -o consult_search

clean:
	rm -f ./minimize ./consult_map ./consult_search
	@echo "Succesfully cleaned."
