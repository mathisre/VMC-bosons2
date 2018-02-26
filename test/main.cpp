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

<<<<<<< HEAD
int main() {
    int numberOfDimensions  = 1;
    int numberOfParticles   = 40;
    int numberOfSteps       = (int) 1e5;
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 1.0;          // Oscillator frequency z-direction
    double a_ho             = 1-2e-4;
    double alpha            = 1.0/(2.0); //*a_ho*a_ho);          // Variational parameter.
    double beta             = 1;          // beta
    double trapSize         = 0;            // trap size
    double stepLength       = 0.1;          // Metropolis step length.
    double equilibration    = 0.1;          // Amount of the total steps used
    // for equilibration.
=======
/* Notes:
 * For 100 particles the analytical derivative gives 52.5 energy and the numerical one gives 7ish.
 * For this, the numerical one is ~10^6 times faster than the analytical ????
 *
 * Program only uses 11 MC steps for analytical but
 */


int main() {
    int numberOfDimensions  = 1;
    int numberOfParticles   = 100;
    int numberOfSteps       = (int) 20;
    double omega            = 1.0;             // Oscillator frequency.<<<<<<< Updated upstream
    double omega_z          = 1.0;             // Oscillator frequency z-direction
    double a_ho             = 1-2e-4;
    double alpha            = 1/(2*a_ho*a_ho); // Variational parameter.
    double beta             = 1;               // beta=======
    double trapSize         = 0;               // trap size
    double stepLength       = 0.1;             // Metropolis step length.
    double equilibration    = 0.1;             // Amount of the total steps used for equilibration.
>>>>>>> 0f9d2e4f2ae1ceead194dde21c5b4e2377e5cdc8

    System* system = new System();

    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));
    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));
    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, trapSize));
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);
    system->runMetropolisSteps          (numberOfSteps);
    return 0;
}
