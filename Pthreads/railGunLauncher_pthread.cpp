#include <iostream>
#include <time.h>
#include <cmath>
#include <pthread.h>
using namespace std;

#define G 9.8
#define NUM_THREADS 4

/* data only read by functions */
struct launcher {
    int length;
    int angle;
    double stepSize;
    int iterations;
} Launcher;

/* struct for each thread */
struct ThreadData {
    int left_idx;
    int right_idx;
    double result;
};

void* threadCompute(void* arg);
double GetVelocity(double a, double b);
double IntegrateAccel(double time);
double SIN(double x);

int main(int argc, char* argv[]) {
    if (argc != 6) {
        cerr << "Usage: " << argv[0] << " <Ramp Length> <Angle> <Step Size> <Sin() Iterations> <num threads>" << endl;
        exit(-1);
    }

    struct timespec start, end;
    double finalVel = 0.0, altitude, time_taken;
    Launcher.length = atof(argv[1]);
    Launcher.angle = atoi(argv[2]);
    Launcher.stepSize = atof(argv[3]);
    Launcher.iterations = atof(argv[4]);
    int thread_count = atoi(argv[5]);

    pthread_t threads[thread_count];
    ThreadData threadData[thread_count];
    int local_n = (Launcher.length / thread_count);
    int remain_n = (Launcher.length % thread_count);

    clock_gettime(CLOCK_MONOTONIC, &start);

    for (int i = 0; i < thread_count; i++) {
        threadData[i].left_idx = i * local_n + 1;
        threadData[i].right_idx = threadData[i].left_idx + local_n - 1;
        if (i == thread_count - 1) threadData[i].right_idx += remain_n;
        pthread_create(&threads[i], NULL, threadCompute, (void*)&threadData[i]);
    }

    for (int i = 0; i < thread_count; i++) {
        pthread_join(threads[i], NULL);
        finalVel += threadData[i].result;
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    altitude = Launcher.length * SIN(Launcher.angle * (M_PI / 180));

    printf("Final Velocity: %15.14fm/s, Altitude: %fm, Angled at: %d degrees\n",
        finalVel, altitude, Launcher.angle);
    printf("Execution time: %f seconds\n", time_taken);

    return 0;
}


void* threadCompute(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    double localVel = 0.0;
    for (int i = data->left_idx; i <= data->right_idx; i++) {
        double local_a = (double)(i - 1);
        double local_b = (double)(i);
        localVel += GetVelocity(local_a, local_b);
    }
    data->result = localVel;
    pthread_exit(NULL);
}

double GetVelocity(double a, double b) {
    double k1, k2, k3, k4;
    double h = Launcher.stepSize;
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

double IntegrateAccel(double time) {
    static double ascale = 10.236589076381454;
    static double tscale = Launcher.length / (M_PI / 2);
    double flatAccel = SIN(time / tscale) * ascale;
    double rampAccel = flatAccel - G * SIN(Launcher.angle);
    return rampAccel;
}

double SIN(double x) {
    double term = x, sum = x;
    for (int i = 2; i < Launcher.iterations * 2; i = (i + 2)) {
        term = -term * (x * x) / ((double)(i + 1) * i);
        sum += term;
    }
    return sum;
}
