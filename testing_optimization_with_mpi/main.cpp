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
#include "mpi.h"

using namespace std;

/* Notes:
 *
 */


int main(int argc, char* argv[]){
    int numprocs, myrank, numberOfParticles, numberOfDimensions;
    double alpha;

    if (argc < 2) {
        cout << "-------------------------------------------------------" << endl
             << "Type dimensions, particles and alpha." << endl
             << "-------------------------------------------------------" << endl;
        exit(EXIT_FAILURE);
    }
    if ( argc >= 2) {
        numberOfParticles   = atoi(argv[2]);
        numberOfDimensions  = atoi(argv[1]);
        alpha = atof(argv[3]);
    }


    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);




    double timeStep         = 0.01;        // Importance sampling time step

    int numberOfSteps       = (int) 1e+6;
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 2.82843;     // Oscillator frequency z-direction


    double a_ho             = 1.5e-4;
//    double alpha            = 0.5;//0.775;    //*a_ho*a_ho);          // Variational parameter.
    double beta             = 2.82843;            // beta
    double interactionSize   =0.0043;            // trap size
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = 0.7;          // Amount of the total steps used for equilibration.

    string filename         = "alpha_" + to_string(alpha) + "_dim_" + to_string(numberOfDimensions) + "_particles_" + to_string(numberOfParticles) + "_" + to_string(myrank) + ".dat";
    // Set filename to "0" to stop from writing to file

    System* system = new System();
    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));

    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, interactionSize, timeStep));

    system->openDataFile                (filename);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);




    //system->setConjugateGradient(new conjugateGradient(system, alphaZero, CJsteps));

    auto start = std::chrono::system_clock::now();
    system->runMetropolisSteps          (numberOfSteps);
    system->printOut(myrank, numprocs);


    auto end = std::chrono::system_clock::now();

    MPI_Finalize();

    if (myrank == 0) {
        std::chrono::duration<double> diff = end-start;
        std::cout << "Time " << diff.count() << " s\n"; //display run time
    }

    return 0;
}
