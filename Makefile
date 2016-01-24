INCLUDE_PATHS = -I lib/ -I src/

build: 
	g++ src/*.cpp -O3 -o fluid_solver $(INCLUDE_PATHS) -L /usr/lib64/atlas -lcblas -fopenmp
