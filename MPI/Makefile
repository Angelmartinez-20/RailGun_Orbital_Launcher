CXX = mpic++
CXXFLAGS = -Wall -std=c++11

all: railGunLauncher_mpi

railGunLauncher_mpi: railGunLauncher_mpi.o
	$(CXX) $(CXXFLAGS) -o $@ $^

railGunLauncher_mpi.o: railGunLauncher_mpi.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f railGunLauncher_mpi railGunLauncher_mpi.o
