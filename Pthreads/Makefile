CXX = g++
CXXFLAGS = -Wall -std=c++11 -lpthread

all: railGunLauncher_pthread

railGunLauncher_pthread: railGunLauncher_pthread.o
	$(CXX) $(CXXFLAGS) -o $@ $^

railGunLauncher_pthread.o: railGunLauncher_pthread.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f railGunLauncher_pthread railGunLauncher_pthread.o
