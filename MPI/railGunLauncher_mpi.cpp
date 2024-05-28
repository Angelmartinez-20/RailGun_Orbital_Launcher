#include <iostream>
#include <mpi.h>
#include <cmath>
using namespace std; 

#define G 9.8

struct launcher {
    int length;
    int angle;
    double stepSize;
    int iterations;
} Launcher;

double GetVelocity(double a, double b);
double IntegrateAccel(double time);
double SIN(double x);

int main(int argc, char* argv[]) {    
    /* ensures correct input */
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << "<Ramp Length> <Angle> <Step Size> <Sin() Iterations>" << endl;
        exit(-1);
    }

    Launcher.length = atof(argv[1]);
    Launcher.angle = atoi(argv[2]);
    Launcher.stepSize = atof(argv[3]); 
    Launcher.iterations = atof(argv[4]);
    int my_rank, comm_sz, local_n, remain_n;  
    double localVel = 0.0, finalVel = 0.0;
    double start_time, end_time; 
    double left_idx, right_idx;

    MPI_Init(NULL, NULL);                    // systems starts MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // ges my process rank
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); // finds how #processes being used

    /* local integration bounds */
    local_n = (Launcher.length / comm_sz);
    remain_n = (Launcher.length % comm_sz);
    left_idx = my_rank * local_n + 1;
    right_idx = left_idx + local_n - 1;
    if (my_rank == comm_sz - 1) right_idx += remain_n;  

    start_time = MPI_Wtime();
    /* each treads finds its local velocity */
    for (int i = left_idx; i <= right_idx; i++) {
        double local_a = (double)(i - 1);
        double local_b = (double)(i);
        localVel += GetVelocity(local_a, local_b);   
    }
    MPI_Reduce(&localVel, &finalVel, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    /* Print the result */
    if (my_rank == 0) {
        double altitude = Launcher.length*SIN(Launcher.angle*(M_PI/180));
        printf("Final Velocity: %15.14fm/s, Altitude: %fm, Angled at: %d degrees\n", 
                finalVel, altitude, Launcher.angle);
        printf("Execution time: %f seconds\n", end_time - start_time);
    }

    MPI_Finalize();

    return 0;
} 

/* Implements --> [V0 + sigma A(t)dt] using Runge-Kutta method */
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
