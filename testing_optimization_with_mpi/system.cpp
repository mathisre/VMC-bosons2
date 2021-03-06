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
#//include "conjugategradient.h"

using namespace std;

bool System::metropolisStep() {
    int randparticle=Random::nextInt(m_numberOfParticles);

    vector <double> r_old=m_particles.at(randparticle).getPosition();
    vector <double> r_new(m_numberOfDimensions);
    vector <vector<double>> QF_old(m_numberOfDimensions,vector<double>(m_numberOfParticles));
vector <vector<double>> QF_new(m_numberOfDimensions,vector<double>(m_numberOfParticles));
QF_old=m_QuantumForce;

    for(int d = 0; d < m_numberOfDimensions; d++){
//        r_new[d] = r_old[d] + m_stepLength*(Random::nextDouble()-0.5);
        r_new[d] = r_old[d] + m_QuantumForce[d][randparticle]*m_timeStep +  m_sqrtTimeStep*(Random::nextGaussian(0,1));
   //cout<<r_new[d]<<endl;
    }


    // Remove energy and quantum force from that single particle and add back when either accept or decline move
    m_particles.at(randparticle).setPosition(r_new);
    getSampler()->updateEnergy(-getHamiltonian()->LocalEnergySingleParticle(m_particles, randparticle));
    setQuantumForce(m_waveFunction->QuantumForce(m_particles));
QF_new=m_QuantumForce;
    bool bosonsTooClose = updateDistanceMatrix(m_particles, randparticle);
    if ( bosonsTooClose == true){ // Then don't accept
        m_particles.at(randparticle).setPosition(r_old);
        updateDistanceMatrix(m_particles, randparticle);
        cout << "hello bosons" << endl;
        getSampler()->updateEnergy(getHamiltonian()->LocalEnergySingleParticle(m_particles, randparticle));
        setQuantumForce(QF_old);
        return false;
    }


    double GreensFunction =0.0;

    for(int d = 0; d < m_numberOfDimensions; d++){
        for(int j=0; j<m_numberOfParticles; j++){
            GreensFunction += 0.5*(QF_old[d][j]+QF_new[d][j])*(0.5*m_timeStep*0.5*(QF_old[d][j]+QF_new[d][j]-r_new[d]+r_old[d]));
        }
    }
    GreensFunction = exp(GreensFunction);


    double psi_new = m_waveFunction->evaluate(m_particles);
    if (Random::nextDouble() <= GreensFunction*psi_new * psi_new / (m_psiOld * m_psiOld)){ // Accept
        m_psiOld = psi_new;
   //    getSampler()->updateEnergy(getHamiltonian()->LocalEnergySingleParticle(m_particles, randparticle));
        //updateQuantumForce(m_waveFunction->QuantumForceSingleParticle(m_particles, randparticle), false);
        return true;
    }
    else{ // Don't accept
        m_particles.at(randparticle).setPosition(r_old);
       // getSampler()->updateEnergy(getHamiltonian()->LocalEnergySingleParticle(m_particles, randparticle));
    setQuantumForce(QF_old);
    //    updateQuantumForce(m_waveFunction->QuantumForceSingleParticle(m_particles, randparticle), false);
        updateDistanceMatrix(m_particles, randparticle);
        return false;
    }
}

void System::runMetropolisSteps(int numberOfMetropolisSteps) {
    m_particles                 = m_initialState->getParticles();
    m_sampler                   = new Sampler(this);
    m_numberOfMetropolisSteps   = numberOfMetropolisSteps;
    //cout<<getSampler()->getEnergy()<<endl;
    getSampler()->setStepNumber(0);
    getSampler()->setAcceptedNumber(0);
    m_sampler->setNumberOfMetropolisSteps(numberOfMetropolisSteps);

    // Initial values
    setDistanceMatrix(computematrixdistance(m_particles));
    m_psiOld = m_waveFunction->evaluate(m_particles);
    getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(m_particles));
    cout<<getSampler()->getEnergy()<<endl;
    setQuantumForce(m_waveFunction->QuantumForce(m_particles));
    for (int i=0; i < numberOfMetropolisSteps; i++) {
        bool acceptedStep = metropolisStep();
        //cout<<m_sampler->getEnergy()<<endl;
            m_sampler->sample(acceptedStep);
            m_sampler->writeToFile();
    }

}

//void System::runConjugateGradient(){
//    m_conjugateGradient->conjugateGradientSteps();
//}

void System::printOut(int myrank, int numprocs)
{
    m_sampler->computeAverages(myrank, numprocs);
    if (myrank == 0) m_sampler->printOutputToTerminal(numprocs);
}



double System::computedistance(int i){
    double temp=0;
    for(int j=0;j<m_numberOfDimensions;j++){
        temp+=m_particles.at(i).getPosition()[j]*m_particles.at(i).getPosition()[j];
    }
    return sqrt(temp);
}

bool System::updateDistanceMatrix( std::vector<class Particle> &particles, int randparticle){
    double temp = 0;
    for (int j = 0; j < randparticle; j++){
        temp = 0;
        for (int d = 0; d < m_numberOfDimensions; d++){
            temp += (particles.at(randparticle).getPosition()[d] - particles.at(j).getPosition()[d]) *
                    (particles.at(randparticle).getPosition()[d] - particles.at(j).getPosition()[d]);
        }
        m_distanceMatrix[randparticle][j] = sqrt(temp);
        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
            return true;
        }
        m_distanceMatrix[j][randparticle] = m_distanceMatrix[randparticle][j];
    }
    for (int j = randparticle+1; j < m_numberOfParticles; j++){
        temp = 0;
        for (int d = 0; d < m_numberOfDimensions; d++){
            temp += (particles.at(randparticle).getPosition()[d] - particles.at(j).getPosition()[d]) *
                    (particles.at(randparticle).getPosition()[d] - particles.at(j).getPosition()[d]);
        }
        m_distanceMatrix[randparticle][j] = sqrt(temp);
        if (m_distanceMatrix[randparticle][j] < getinteractionSize()){
            return true;
        }
        m_distanceMatrix[j][randparticle] = m_distanceMatrix[randparticle][j];

    }
    return false;
}

std::vector<vector<double>> System::computematrixdistance(std::vector<class Particle> &particles){

    vector<vector<double>> distancematrix(m_numberOfParticles, vector<double>(m_numberOfParticles));
    double temp=0;
    int j=0;
    while(j < m_numberOfParticles){
        temp = 0;
        for(int i = 0; i < j; i++){

            for(int k=0;k<m_numberOfDimensions;k++){
                temp+=(particles.at(i).getPosition()[k] - particles.at(j).getPosition()[k]) *
                      (particles.at(i).getPosition()[k] - particles.at(j).getPosition()[k]);
            }
            distancematrix[i][j]=sqrt(temp);
            distancematrix[j][i]=distancematrix[i][j];
        }

        j++;
    }

    return distancematrix;
}

void System::updateQuantumForce(std::vector<std::vector<double> > deltaQuantumForce, bool subtract){
    if (subtract == false){
        for (int d = 0; d < m_numberOfDimensions; d++){
        for(int i=0; i<=m_numberOfParticles; i++){
            m_QuantumForce[d][i] += deltaQuantumForce[d][i];
        }
        }
    }
    else {
        for (int d = 0; d < m_numberOfDimensions; d++){
        for(int i=0; i<=m_numberOfParticles; i++){
            m_QuantumForce[d][i] += deltaQuantumForce[d][i];
        }
        }
    }
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

double System::getinteractionSize() const
{
    return m_interactionSize;
}

void System::setinteractionSize(double interactionSize)
{
    m_interactionSize = interactionSize;
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

std::vector<vector<double> > System::getDistanceMatrix() const
{
    return m_distanceMatrix;
}

double System::getDistanceMatrixij(int i, int j) const
{
    return m_distanceMatrix[i][j];
}

double System::getPsiOld() const
{
    return m_psiOld;
}

void System::setPsiOld(double psiOld)
{
    m_psiOld = psiOld;
}

std::vector<vector<double>> System::getQuantumForce() const
{
    return m_QuantumForce;
}

void System::setQuantumForce(const std::vector<vector<double>> &QuantumForce)
{
    m_QuantumForce = QuantumForce;
}

void System::setDistanceMatrix(const std::vector<vector<double> > &distanceMatrix)
{
    m_distanceMatrix = distanceMatrix;
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


//void System::setConjugateGradient(conjugateGradient* conjugateGradient)
//{
//    m_conjugateGradient = conjugateGradient;
//}

double System::findEnergyDerivative(double CJsteps)
{

    double meanEnergy      = getSampler()->getEnergy() / CJsteps;
    double meanWFderiv     = getSampler()->getWFderiv() / CJsteps;
    double meanWFderivEloc = getSampler()->getWFderivMultELoc() / CJsteps;




    // Make the sampler sample the wavefunc deriv and things like that
    // Then just find the mean and we are good
//cout<<meanWFderivEloc<<"meanwfderiveloc"<<endl;
//cout<<-meanEnergy*meanWFderiv<<"other"<<endl;
//if(2 * (meanWFderivEloc - meanEnergy*meanWFderiv)<0) return -2 * (meanWFderivEloc - meanEnergy*meanWFderiv);
 return 2 * (meanWFderivEloc - meanEnergy*meanWFderiv);
}

