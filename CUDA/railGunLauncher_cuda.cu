#include <iostream>
#include <cuda.h>
#include <time.h>
#include <cmath>
using namespace std;

#define G 9.8

struct launcher {
    int length;
    int angle;
    double stepSize;
    int iterations;
} Launcher;

__constant__ launcher d_launcher;
__device__ double GetVelocity(double a, double b);
__device__ double IntegrateAccel(double time);
__device__ double SIN(double x);

__global__ void computeVelocity(double* velocities) {
    int i = threadIdx.x + blockIdx.x * blockDim.x + 1;
    if (i <= d_launcher.length) {
        double local_a = (double)(i - 1);
        double local_b = (double)(i);
        velocities[i-1] = GetVelocity(local_a, local_b);
    }
}

int main(int argc, char* argv[]) {
    /* ensures correct input */
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " <Ramp Length> <Angle> <Step Size> <Sin() Iterations>" << endl;
        exit(-1);
    }

    struct timespec start, end;
    double finalVel = 0.0, altitude, time_taken;
    Launcher.length = atoi(argv[1]);
    Launcher.angle = atoi(argv[2]);
    Launcher.stepSize = atof(argv[3]); 
    Launcher.iterations = atoi(argv[4]); // Corrected to atoi for integer conversion

    /* Copy Launcher struct from host to device */
    cudaMemcpyToSymbol(d_launcher, &Launcher, sizeof(launcher));

    /* gets device properties */
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    dim3 threadsPerBlock(prop.maxThreadsPerBlock);
    dim3 numBlocks((Launcher.length + threadsPerBlock.x - 1) / threadsPerBlock.x);
    cout << "Number of Blocks: " << numBlocks.x << endl;
    cout << "Threads per Block: " << threadsPerBlock.x << endl;  
      
    clock_gettime(CLOCK_MONOTONIC, &start); 

    /* allocate mem and run on GPU */
    double* d_velocities;
    cudaMalloc((void**)&d_velocities, Launcher.length * sizeof(double));
    computeVelocity<<<numBlocks, threadsPerBlock>>>(d_velocities);

    /* allocate mem and compute on CPU */
    double* h_velocities = (double*)malloc(Launcher.length * sizeof(double));
    cudaMemcpy(h_velocities, d_velocities, Launcher.length * sizeof(double), cudaMemcpyDeviceToHost);
    for (int i = 0; i < Launcher.length; ++i) {
        finalVel += h_velocities[i];
    }

    cudaFree(d_velocities);
    free(h_velocities);

    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    altitude = Launcher.length*sin(Launcher.angle*(M_PI/180));

    printf("Final Velocity: %15.14fm/s, Altitude: %fm, Angled at: %d degrees\n", 
            finalVel, altitude, Launcher.angle);
    printf("Execution time: %f seconds\n", time_taken);

    return 0;
}


/* Implements --> [V0 + sigma A(t)dt] using Runge-Kutta 4 method */
__device__ double GetVelocity(double a, double b) {
    double k1, k2, k3, k4;
    double h = d_launcher.stepSize;
    double y = 0.0;
    int numSteps = (int)((b - a) / h);

    for (int i = 0; i < numSteps; i++) {
        double x0 = a + i * h;
        k1 = h * IntegrateAccel(x0);
        k2 = h * IntegrateAccel(x0 + 0.5 * h);
        k3 = h * IntegrateAccel(x0 + 0.5 * h);
        k4 = h * IntegrateAccel(x0 + h);
        y += (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    }

    return y;
}

__device__ double IntegrateAccel(double time) {
    double ascale = 10.236589076381454;                     // increase amplitude of sine curve
    double tscale = d_launcher.length/(M_PI/2);             // gets first 1/4 sine curve
    double flatAccel = SIN(time/tscale)*ascale;             // acceleration when flat ground
    double rampAccel = flatAccel-G*SIN(d_launcher.angle);   // acceleration on ramp (against earth gravity)
    return rampAccel;
}

__device__ double SIN(double x) {
    double term = x, sum = x;
    /* finds sin() with taylor series */
    for (int i = 2; i < d_launcher.iterations*2; i += 2) {
        term = -term * (x*x) / ((double)(i + 1) * i);
        sum += term;
    }
    return sum;
}