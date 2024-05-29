CXX = nvcc
CXXFLAGS = -std=c++11 -Xcompiler -Wall

all: railGunLauncher_cuda

railGunLauncher_cuda: railGunLauncher_cuda.o
	$(CXX) $(CXXFLAGS) -o $@ $^

railGunLauncher_cuda.o: railGunLauncher_cuda.cu
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f railGunLauncher_cuda railGunLauncher_cuda.o
