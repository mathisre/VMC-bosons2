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

using namespace std;

/* Notes:
 * Program is not using all the steps? When using 20 MC steps the program
 * says it is using 10^1.3 for some reason.
 * For 10 particles and 20 steps, the numerical is ~15 times faster than the analytical method.
 * Analytical derivative gives 4.7 energy and numerical gives 3.2
 * Analytical 75000 ms and numerical 4860 ms
 */


int main() {
    int numberOfDimensions  = 1;
    int numberOfParticles   = 3;
    int numberOfSteps       = (int) 1e6;
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 1.0;          // Oscillator frequency z-direction

    double a_ho             = 1-2e-4;
    double alpha            = 1.0/(2.0); //*a_ho*a_ho);          // Variational parameter.
    double beta             = 1;          // beta
    double trapSize         = 0;            // trap size
    double stepLength       = 0.1;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used
    // for equilibration.

    System* system = new System();

    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));
    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, trapSize));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (numberOfSteps);
    return 0;
}
