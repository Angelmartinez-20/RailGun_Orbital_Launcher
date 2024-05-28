CXX = g++
CXXFLAGS = -Wall -std=c++11

all: railGunLauncher

railGunLauncher: railGunLauncher.o
	$(CXX) $(CXXFLAGS) -o $@ $^

railGunLauncher.o: railGunLauncher.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f railGunLauncher railGunLauncher.o
