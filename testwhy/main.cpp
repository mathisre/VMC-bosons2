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
    int numberOfParticles   = 1;
    int numberOfDimensions  = 1;
    int numberOfSteps       = (int) 1e5;
    double omega            = 1.0;          // Oscillator frequency.
    double omega_z          = 1.0;          // Oscillator frequency z-direction
    double timeStep         = 0.001;        // Importance sampling time step

    double a_ho             = 1-2e-4;
    double alpha            = 1.0/(2.0);    //*a_ho*a_ho);          // Variational parameter.
    double beta             = 1;            // beta
    double trapSize         = 0;            // trap size
    double stepLength       = 1.0;          // Metropolis step length.
    double equilibration    = 0.7;          // Amount of the total steps used for equilibration.

    string filename         = "0";          // Set equal to "0" if you don't want any data

    // Optimalization of alpha using steepest descent method
    int CJsteps       = (int) 1e5;    // Number of steps MC steps
    double alphaZero        = 0.5;         // Initial guess

    System* system = new System();

    system->setHamiltonian              (new HarmonicOscillator(system, omega, omega_z));

    system->setWaveFunction             (new SimpleGaussian(system, alpha, beta));

    system->setInitialState             (new RandomUniform(system, numberOfDimensions, numberOfParticles, trapSize, timeStep));
    system->openDataFile                (filename);
    system->setEquilibrationFraction    (equilibration);
    system->setStepLength               (stepLength);

    clock_t startTime, endTime;

    startTime = clock();
    system->runMetropolisSteps          (numberOfSteps);
    system->printOut();

    endTime = clock();

   // cout << " Computation time = " << difftime(endTime, startTime) / CLOCKS_PER_SEC << " s" << endl;
    return 0;
}

/*
 *     double energyDerivative, meanEnergy;
    double energyDerivative_2;
    double iterations = 0;
    double tolerance =1e-10;
    double max_iterations = 1e3;
    double matrix=1;
   // double tol = (double) 1e-10;

    system->setWaveFunction(new SimpleGaussian(system, alphaZero, beta));
    system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles, trapSize, timeStep));
    system->runMetropolisSteps(CJsteps);
    energyDerivative_2 =system->findEnergyDerivative(CJsteps);

    double P_0=-energyDerivative_2;
   // cout<<P_0<<endl;
    double P=0;
    double psycho=0;
    double s=P_0*1;
  alphaZero+=s;


    double y=0;
energyDerivative=1;
    while (fabs(energyDerivative) > tolerance  && iterations<max_iterations)
    {

        system->setWaveFunction(new SimpleGaussian(system, alphaZero, beta));
        system->setInitialState(new RandomUniform(system, numberOfDimensions, numberOfParticles, trapSize, timeStep));
        system->runMetropolisSteps(CJsteps);
        energyDerivative=system->findEnergyDerivative(CJsteps);

        y=energyDerivative-energyDerivative_2;
        matrix=(y)/s;
        P=- energyDerivative/matrix;
        s=P/(y*s);
      // cout<<s<<endl;
        alphaZero+=s;
        energyDerivative_2=energyDerivative;
        iterations++;
        cout<<alphaZero<<endl;
    }
    //cout<<iterations<<endl;
    */