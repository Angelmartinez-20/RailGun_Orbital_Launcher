CXX = g++
CXXFLAGS = -Wall -std=c++11 -fopenmp

all: railGunLauncher_omp

railGunLauncher_omp: railGunLauncher_omp.o
	$(CXX) $(CXXFLAGS) -o $@ $^

railGunLauncher_omp.o: railGunLauncher_omp.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f railGunLauncher_omp railGunLauncher_omp.o
