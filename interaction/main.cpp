#include <iostream>
#include <random>
#include <cmath>
#include "system.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include "conjugategradient.h"
#include <chrono>
#include <string>

using namespace std;

/* Notes:
 *
 */


int main(int argc, char* argv[]){

//    if (argc < 3){

//        cout << "Put dimensions and particles in cmd line" << endl;
//        exit(EXIT_FAILURE);
//    }
//    int numberOfParticles   = atoi(argv[2]);
//    int numberOfDimensions  = atoi(argv[1]);




    int numberOfParticles   = 10;
    int numberOfDimensions  = 3;
    double alpha            = 0.5;      // Variational parameter.

    double timeStep         = 0.01;        // Importance sampling time step

    int numberOfSteps       = (int) 1e+5;
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 2.82843;     // Oscillator frequency z-direction


    double a_ho             = 1.5e-4;
    double beta             = 2.82843;            // beta 2.82843
    double interactionSize  = 0.00043;            // trap size
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = 0.7;          // Amount of the total steps used for equilibration.



    double bucketSize = 0.01;
    int bins = ceil(4 / bucketSize);

    string filename         = "alpha_" + to_string(alpha) + "_dim_" + to_string(numberOfDimensions) + "_particles_" + to_string(numberOfParticles)  + ".dat";

    // Set filename to "0" to stop from writing to file

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));

    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, interactionSize, timeStep, bins, bucketSize));

    system->openDataFile                (filename);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);




    //system->setConjugateGradient(new conjugateGradient(system, alphaZero, CJsteps));

    auto start = std::chrono::system_clock::now();
    system->runMetropolisSteps          (numberOfSteps);
    system->printOut();

//    system->printOneBodyDensity();

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << " Computation time = " << diff.count() << " s\n" << endl; //display run time

    return 0;
}