# RailGun Orbital Launcher Simulator

## Introduction

For this project, I am simulating a Railgun Orbital Launcher. This is a method to launch items 
into space without using rocket launcher but instead a rail track system using electromagnetic 
forces to launch an item from a incline rail track. The important values are the final velocity 
of the object when it has left the ramp as that velocity determines if the object can break 
through the earth's atmosphere. Additionally, altitude is calculated so that the path of the 
object can be determined later. The acceleration function that this simulator uses is the 1st 
quarter of a sine wave which has been modified to adjust its length and amplitude depending on 
input. Also the gravitational force being applied to the object as it goes uphill is calculated 
considered aswell. 

Numerical methods used include the 4th order Runge-Kutta numerical integration to find the 
antiderivative of the acceleration function to get velocity, and Taylor series to calculate the 
sine function value which is part of the acceleration function. The different parallel methods 
utilized are Pthreads, OpenMP, MPI, and CUDA.

## Instructions

To use the simulator, each directory contains its own `Makefile` which is used to compile onto 
your system. The arguments needed for the `Sequential` and `CUDA` version are
```bash
./railGunLauncher <Ramp Length> <Angle> <Step Size> <Sin() Iterations>
```
Then for the `Pthreads` or `OpenMP` versions
```bash
./railGunLauncher <Ramp Length> <Angle> <Step Size> <Sin() Iterations> <Num Threads>
```
Lastly to run `MPI` on only your host machine, use
```bash
mpiexec -n <Num Threads> ./railGunLauncher <Ramp Length> <Angle> <Step Size> <Sin() Iterations>
```
`<Ramp Length>` is the length of the rail system you are simulating 

`<Angle>` is the incline angle in degrees of your ramp

`<Step Size>` is the size you want to break down each meter for the numeric integration

`<Sin() Iterations>` is the number of iterations to run the taylor series SIN() function

`<Num Threads>` is the number of threads you want to split the workload into

## Analysis

An analysis was conducted on the different parallel methods using 2, 4, & 8 threads on a 
machine with 12 cores. The **Average Time** of each script was calculated by running it 5 times 
and getting the average. The **Accuracy** was determined by comparing the result of the script 
to the antiderivative of the acceleration function using WolframAlpha. **SU** is the speed up of 
the parallel script calculated by the equation $`SU = T_{sequential} / T_{parallel}`$ . Lastly **P** is 
the percentage of the script that is run in parallel which is calculated using the equation 
$`P = \frac{S(1 - \frac{1}{SU})}{S - 1}`$ . 

<div align="center">
  
| Script     | Threads | Average Time | Accuracy      | SU       | P        |
|------------|---------|--------------|---------------|----------|----------|
| Sequential | 1       | 28.470565    | -0.0355033746 | 1        | 0        |
| Pthreads   | 2       | 15.525461    | -0.0355033741 | 1.833798 | 0.909368 |
| OpenMP     | 2       | 16.420991    | -0.0355033741 | 1.733791 | 0.846458 |
| MPI        | 2       | 15.356219    | -0.0355033741 | 1.854009 | 0.921256 |
| Pthreads   | 4       | 9.783226     | -0.0355033740 | 2.910141 | 0.875165 |
| OpenMP     | 4       | 10.609861    | -0.0355033740 | 2.683406 | 0.836452 |
| MPI        | 4       | 10.92802     | -0.0355033740 | 2.605281 | 0.821552 |
| Pthreads   | 8       | 6.262513     | -0.0355033740 | 4.546189 | 0.891469 |
| OpenMP     | 8       | 5.970608     | -0.0355033740 | 4.768454 | 0.903187 |
| MPI        | 8       | 6.796102     | -0.0355033740 | 4.189249 | 0.87005  |

</div>

**Note:** Because CUDA isn't meant to run on the same number of threads as CPU cores,
I did not do the same analysis. I did however ran a thread for each meter inwhich the
SpeedUp ended being 101 times faster and the Accuracy ended up being of by -0.0355033741m











