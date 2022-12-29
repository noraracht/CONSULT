all: minimize map search
	
minimize:
	g++ minimize.cpp -std=c++11 -o minimize

map:
	g++ map.cpp -std=c++11 -O3 -o map

search:
	g++ search.cpp -std=c++11 -fopenmp -O3 -o search

clean:
	rm -f ./minimize ./map ./search
	@echo "Succesfully cleaned."
