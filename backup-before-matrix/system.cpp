#include "system.h"
#include <cassert>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "particle.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/hamiltonian.h"
#include "InitialStates/initialstate.h"
#include "Math/random.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include "conjugategradient.h"

using namespace std;

bool System::metropolisStep() {
    /* Perform the actual Metropolis step: Choose a particle at random and
     * change it's position by a random amount, and check if the step is
     * accepted by the Metropolis test (compare the wave function evaluated
     * at this new position with the one at the old position).
     */
    int randparticle=Random::nextInt(m_numberOfParticles);

    vector <double> r_old=m_particles.at(randparticle).getPosition();
    double psi_old=m_waveFunction->evaluate(m_particles);
    vector <double> r_new(m_numberOfDimensions);
    vector <double> QuantumForce = m_waveFunction->QuantumForce(m_particles);

    for(int d = 0 ; d < m_numberOfDimensions;d++){
        //r_new[d] = r_old[d] + m_stepLength*(Random::nextDouble()-0.5);
        //r_new[d] = r_old[d] +  0.5 * QuantumForce[d]*m_timeStep +  m_sqrtTimeStep*(Random::nextDouble()-0.5);

        r_new[d] = r_old[d] + QuantumForce[d]*m_timeStep +  m_sqrtTimeStep*(Random::nextDouble()-0.5);
    }
    m_particles.at(randparticle).setPosition(r_new);
    double psi_new=m_waveFunction->evaluate(m_particles);

    double random_number=Random::nextDouble();

    if (random_number <= psi_new * psi_new / (psi_old * psi_old)){
        return true;
    }
    else m_particles.at(randparticle).setPosition(r_old); return false;
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    getSampler()->setStepNumber(0);
    getSampler()->setAcceptedNumber(0);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);


    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();

        /* Here you should sample the energy (and maybe other things using
         * the m_sampler instance of the Sampler class. Make sure, though,
         * to only begin sampling after you have let the system equilibrate
         * for a while. You may handle this using the fraction of steps which
         * are equilibration steps; m_equilibrationFraction.
         */
        //if(m_sampler->getStepNumber()/m_sampler->getNumberOfMetropolisSteps() > 1.0 - m_equilibrationFraction){
            m_sampler->sample(acceptedStep);
            //if (i % 1000 == 0)
                m_sampler->writeToFile();


            //cout << i+1 << endl;

    }

}

void System::runConjugateGradient(){
    m_conjugateGradient->conjugateGradientSteps();
}

void System::printOut()
{
    m_sampler->computeAverages();
    m_sampler->printOutputToTerminal();
}



double System::computedistance(int i){
    double temp=0;
    for(int j=0;j<m_numberOfDimensions;j++){
        temp+=m_particles.at(i).getPosition()[j]*m_particles.at(i).getPosition()[j];
    }
    return sqrt(temp);
}


std::vector<vector<double>> System::computematrixdistance(std::vector<class Particle> &particles){
    vector<vector<double>> distancematrix(m_numberOfParticles, vector<double>(m_numberOfParticles));

    double temp=0;
    int j=0;
    while(j<m_numberOfParticles){
        for(int i=0; i<j; i++){
            for(int k=0;k<m_numberOfDimensions;k++){
                temp+=(m_particles.at(i).getPosition()[k] - m_particles.at(j).getPosition()[k]) *
                      (m_particles.at(i).getPosition()[k] - m_particles.at(j).getPosition()[k]);
            }
            distancematrix[i][j]=sqrt(temp);
            j++;
        }
    }
    while(j<m_numberOfParticles){
        for(int i=0; i<j; i++){
            distancematrix[j][i]=distancematrix[i][j];
            j++;
        }
    }
    return distancematrix;
}



double System::computedistanceABS(int i, int j){
    double temp=0;
    for(int k=0;k<m_numberOfDimensions;k++){
        temp+=(m_particles.at(i).getPosition()[k] - m_particles.at(j).getPosition()[k]) *
                (m_particles.at(i).getPosition()[k] - m_particles.at(j).getPosition()[k]);
    }
    return sqrt(temp);{
    }
}

void System::openDataFile(string filename){
    m_sampler->openDataFile(filename);
}

double System::getTrapSize() const
{
    return m_trapSize;
}

void System::setTrapSize(double trapSize)
{
    m_trapSize = trapSize;
}

double System::getTimeStep() const
{
    return m_timeStep;
}

void System::setTimeStep(double timeStep)
{
    m_timeStep = timeStep;
}

double System::getSqrtTimeStep() const
{
    return m_sqrtTimeStep;
}

void System::setSqrtTimeStep(double sqrtTimeStep)
{
    m_sqrtTimeStep = sqrtTimeStep;
}



void System::setNumberOfParticles(int numberOfParticles) {
    m_numberOfParticles = numberOfParticles;
}

void System::setNumberOfDimensions(int numberOfDimensions) {
    m_numberOfDimensions = numberOfDimensions;}

void System::setStepLength(double stepLength) {
    assert(stepLength >= 0);
    m_stepLength = stepLength;
}

void System::setEquilibrationFraction(double equilibrationFraction) {
    assert(equilibrationFraction >= 0);
    m_equilibrationFraction = equilibrationFraction;
}

void System::setHamiltonian(Hamiltonian* hamiltonian) {
    m_hamiltonian = hamiltonian;
}

void System::setWaveFunction(WaveFunction* waveFunction) {
    m_waveFunction = waveFunction;
}

void System::setInitialState(InitialState* initialState) {
    m_initialState = initialState;
}


void System::setConjugateGradient(conjugateGradient* conjugateGradient)
{
    m_conjugateGradient = conjugateGradient;
}



