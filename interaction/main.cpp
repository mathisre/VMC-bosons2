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
//#include "conjugategradient.h"
#include <chrono>
#include <string>

using namespace std;



int main(int argc, char* argv[]){

//    if (argc < 3){

//        cout << "Put dimensions and particles in cmd line" << endl;
//        exit(EXIT_FAILURE);
//    }
//    int numberOfParticles   = atoi(argv[2]);
//    int numberOfDimensions  = atoi(argv[1]);

    int numberOfParticles   = 10;
    int numberOfDimensions  = 3;
    double alpha            = 0.50018;      // Variational parameter.
    double beta             = 2.82843;            // beta 2.82843

    int numberOfSteps       = (int) 1e+6;
    double timeStep         = 0.01;        // Importance sampling time step
    double interactionSize  = 0;            // 0.00043
    double stepLength       = 1.0;          // Metropolis step length.
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z = beta;
    double equilibration    = 0.7;          // Amount of the total steps used for equilibration.

    // Parameters for onebody density histogram
    double bucketSize = 0.01;
    int bins = ceil(4 / bucketSize);

//    string filename = "0";
    string filename         = "alpha_" + to_string(alpha) + "_dim_" + to_string(numberOfDimensions) + "_particles_" + to_string(numberOfParticles)  + "_optimal.dat";
    // Set filename to "0" to stop from writing to file

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));
    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, interactionSize, timeStep, bins, bucketSize));
    system->openDataFile                (filename);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);



    // Gradient descent method to find energy minimum

//    int maxIterations = 200;
//    double initialAlpha = 0.49;
//    string minFilename = "find_minimum_" + to_string(initialAlpha) +"_N_" + to_string(numberOfParticles) + "_non_inter.dat";
//    alpha = system->gradientDescent(initialAlpha, minFilename, maxIterations);
//    cout << "Optimal alpha found by steepest descent: " << alpha << endl;
//    vector<double> parameters(3);
//    parameters[0] = alpha;
//    parameters[1] = alpha;
//    parameters[2] = alpha*beta;
//    system->getWaveFunction()->setParameters(parameters);


    auto start = std::chrono::system_clock::now();
    system->runMetropolisSteps          (numberOfSteps);
    system->printOut();


    string densityFilename = "density_alpha_" + to_string(alpha) + "_beta_" + to_string(beta) + "_inter.dat";
//    system->printOneBodyDensity(densityFilename);

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << " Computation time = " << diff.count() << " s\n" << endl; //display run time
    return 0;
}
