#include <iostream>
#include <random>
#include <cmath>
#include "system.h"
#include "mpi/mpi.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "WaveFunctions/simplegaussian.h"
#include "Hamiltonians/hamiltonian.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "InitialStates/initialstate.h"
#include "InitialStates/randomuniform.h"
#include "Math/random.h"
#include "conjugategradient.h"
#include <time.h>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
    int numprocs, myrank, numberOfParticles, numberOfDimensions;

    if (argc < 2) {
        cout << "-------------------------------------------------------" << endl
             << "Write number of dimensions (<= 3) and particles." << endl
             << "For example; mpirun -n 4 ./prog.x 2 5" << endl
             << "-------------------------------------------------------" << endl;
        exit(EXIT_FAILURE);
    }
    if ( argc >= 2) {
        numberOfParticles   = atoi(argv[2]);
        numberOfDimensions  = atoi(argv[1]);
    }


    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    int numberOfSteps       = (int) 1e5 / (int) numprocs;
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 1.0;          // Oscillator frequency z-direction
    double timeStep         = 0.001;        // Importance sampling time step

    double a_ho             = 1-2e-4;
    double alpha            = 1.0/(2.0);    //*a_ho*a_ho);          // Variational parameter.
    double beta             = 1;            // beta
    double trapSize         = 0;            // trap size
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = 0.4;          // Amount of the total steps used for equilibration.

    string filename         = "0";          // Set equal to "0" if you don't want any data

    // Optimalization of alpha using steepest descent method
    int CJsteps       = (int) 1e4;    // Number of steps MC steps
    double alphaZero        = 0.25;         // Initial guess

    System* system = new System();

    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));

    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, trapSize, timeStep));
    system->openDataFile                (filename);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    double startTime = MPI_Wtime();



    //system->setConjugateGradient(new conjugateGradient(system, alphaZero, CJsteps));


    //cout << "Rank = " << myrank << endl;

    system->runMetropolisSteps(numberOfSteps);

    system->printOut(myrank, numprocs);

    double endTime = MPI_Wtime();

    MPI_Finalize();

    if (myrank == 0) cout << " Computation time = " << endTime - startTime << " s" << endl;
    return 0;
}
