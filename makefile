all: minimization map search
	
minimization:
	g++ minimization.cpp -std=c++11 -o minimization

map:
	g++ main_map.cpp -std=c++11 -O3 -o main_map


search:
	g++ main_search.cpp -std=c++11 -fopenmp -O3 -o main_search

clean:
	rm -f ./minimization ./main_map ./main_search
	@echo "Succesfully cleaned."
