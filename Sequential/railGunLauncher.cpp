#include <iostream>
#include <time.h>
#include <cmath>
using namespace std; 

#define G 9.8

struct launcher {
    double length;
    int angle;
    double stepSize;
    int iterations;
} Launcher;

double GetVelocity();
double IntegrateAccel(double time);
double SIN(double x);

int main(int argc, char* argv[]) {     
    /* ensures correct input */
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " <Ramp Length> <Angle> <Step Size> <Sin() Iterations>" << endl;
        exit(-1);
    }

    struct timespec start, end;
    double finalVel = 0.0, altitude, time_taken;
    Launcher.length = atof(argv[1]);
    Launcher.angle = atoi(argv[2]);
    Launcher.stepSize = atof(argv[3]); 
    Launcher.iterations = atof(argv[4]);

    clock_gettime(CLOCK_MONOTONIC, &start); 
    finalVel = GetVelocity();
    clock_gettime(CLOCK_MONOTONIC, &end);
 
    time_taken = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    altitude = Launcher.length*SIN(Launcher.angle*(M_PI/180));

    printf("Final Velocity: %15.14fm/s, Altitude: %fm, Angled at: %d degrees\n", 
            finalVel, altitude, Launcher.angle);
    printf("Execution time: %f seconds\n", time_taken);
    
    return 0;
} 

/* Implements --> [V0 + sigma A(t)dt] using Runge-Kutta method */
double GetVelocity(){
    double k1, k2, k3, k4;
    double y = 0;
    double x = Launcher.length;
    double h = Launcher.stepSize;
    int numSteps = (int)(x / h);

    for (int i = 0; i < numSteps; i++){
        double x0 = i * h;
        k1 = h * IntegrateAccel(x0);
        k2 = h * IntegrateAccel(x0 + 0.5 * h);
        k3 = h * IntegrateAccel(x0 + 0.5 * h);
        k4 = h * IntegrateAccel(x0 + h);

        y += (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    return y;
}

double IntegrateAccel(double time) {
    static double ascale = 10.236589076381454;              // increase ampltiude of sine curve
    static double tscale = Launcher.length/(M_PI/2);        // gets first 1/4 sine curve
    double flatAccel = SIN(time/tscale)*ascale;             // acceleration when flat ground
    double rampAccel = flatAccel-G*SIN(Launcher.angle);     // acceleration on ramp (against earth gravity)
    return rampAccel;
}

double SIN(double x) {
    double term = x, sum = x;
    /* finds sin() with taylor series */
    for(int i = 2; i < Launcher.iterations*2; i = (i + 2)) {
        term = -term * (x*x)/ ((double)(i+1)*i);
        sum += term;
    }
    return sum;
}
