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
#include <time.h>
#include <string>

using namespace std;

/* Notes:
 *
 */


int main() {
    int numberOfDimensions  = 1;
    int numberOfParticles   = 10;
    int numberOfSteps       = (int) 1e5;
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 1.0;          // Oscillator frequency z-direction
    double timeStep         = 0.001;        // Importance sampling time step

    double a_ho             = 1-2e-4;
    double alpha            = 1.0/(2.0);    //*a_ho*a_ho);          // Variational parameter.
    double beta             = 1;            // beta
    double trapSize         = 0;            // trap size
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = 0.5;          // Amount of the total steps used for equilibration.

    string filename         = "E_0part.txt";          // Set equal to "0" if you don't want any data

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

    clock_t startTime, endTime;



    //system->setConjugateGradient(new conjugateGradient(system, alphaZero, CJsteps));

    startTime = clock();
    system->runMetropolisSteps          (numberOfSteps);
    system->printOut();

    endTime = clock();

    cout << " Computation time = " << difftime(endTime, startTime) / CLOCKS_PER_SEC << " s" << endl;
    return 0;
}
