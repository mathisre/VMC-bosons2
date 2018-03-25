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

bool System::metropolisStepBruteForce() {
    int randparticle=Random::nextInt(m_numberOfParticles);

    vector <double> r_old=m_particles.at(randparticle).getPosition();
    vector <double> r_new(m_numberOfDimensions);

   for(int d = 0; d < m_numberOfDimensions; d++){
        r_new[d] = r_old[d] + m_stepLength*(Random::nextDouble()-0.5);
//        cout << "R_old = " << r_old[d] << endl;
//        cout << "R_new = " << r_new[d] << endl;
    }

    m_particles.at(randparticle).setPosition(r_new);
    updateDistanceMatrix(m_particles, randparticle);

    double psi_new = m_waveFunction->evaluate(m_particles);
    if (Random::nextDouble() <= psi_new * psi_new / (m_psiOld * m_psiOld)){ // Accept
        m_psiOld = psi_new;

        getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(getParticles()));
        return true;
    }
    else{ // Don't accept
        m_particles.at(randparticle).setPosition(r_old);
        updateDistanceMatrix(m_particles, randparticle);
        return false;
    }
}



bool System::metropolisStepImportance() {
    int randparticle=Random::nextInt(m_numberOfParticles);
    vector <double> r_old=m_particles.at(randparticle).getPosition();
    vector <double> r_new(m_numberOfDimensions);

    vector <vector<double>> QF_old(m_numberOfDimensions,vector<double>(m_numberOfParticles));
    vector <vector<double>> QF_new(m_numberOfDimensions,vector<double>(m_numberOfParticles));
    QF_old=m_QuantumForce;

   for(int d = 0; d < m_numberOfDimensions; d++){
        r_new[d] = r_old[d] + 0.5*m_QuantumForce[d][randparticle]*m_timeStep + m_sqrtTimeStep*(Random::nextGaussian(0, 1));
//        cout << r_new[d] - r_old[d] << endl;
//        cout << "R_old = " << r_old[d] << endl;
//        cout << "R_new = " << r_new[d] << endl;
    }

    m_particles.at(randparticle).setPosition(r_new);
    setQuantumForce(m_waveFunction->QuantumForce(m_particles));
    QF_new = m_QuantumForce;
    updateDistanceMatrix(m_particles, randparticle);

    double GreensFunction =0.0;
    for(int d = 0; d < m_numberOfDimensions; d++){
        int j = randparticle;
        GreensFunction += 0.5*(QF_old[d][j ] + QF_new[d][j])
                *(0.5 * 0.5 * m_timeStep * (QF_old[d][j] - QF_new[d][j]) - r_new[d] + r_old[d]);
        }

    GreensFunction = exp(GreensFunction);
    double psi_new = m_waveFunction->evaluate(m_particles);

    // Accept
    if (Random::nextDouble() <= GreensFunction*psi_new * psi_new / (m_psiOld * m_psiOld)){
        m_psiOld = psi_new;
        getSampler()->setEnergy(getHamiltonian()->computeLocalEnergy(getParticles()));
        return true;
    }

    // Don't accept
    else{
        m_particles.at(randparticle).setPosition(r_old);
        setQuantumForce(QF_old);
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
    setQuantumForce(m_waveFunction->QuantumForce(m_particles));

    setHistogram();

    for (int i=0; i < numberOfMetropolisSteps; i++) {
//        bool acceptedStep = metropolisStepBruteForce();
        bool acceptedStep = metropolisStepImportance();
        //cout<<m_sampler->getEnergy()<<endl;
            m_sampler->sample(acceptedStep);
            m_sampler->writeToFile();

    }

}




void System::oneBodyDensity(){
    vector<int> histogram(m_bins);
    double r2 = 0;
    for (int j = 0; j < getNumberOfParticles(); j++){
        r2 = 0;
        for (int d = 0; d < getNumberOfDimensions(); d++){
            r2 += m_particles.at(j).getPosition()[d]*m_particles.at(j).getPosition()[d];
        }
        r2 = sqrt(r2);
        int bucket = (int)floor(r2 / m_bucketSize);
        histogram[bucket] += 1;
    }

    // Update histogram
    for (int k = 0; k < m_bins; k++){
       m_histogram[k] += histogram[k];
    }
}

void System::printOneBodyDensity(string filename){
    ofstream myFile;
    myFile.open(filename);
    for (int j = 0; j < m_bins; j++){
        myFile <<(double) m_histogram[j] /(getNumberOfParticles()*getNumberOfMetropolisSteps()*getEquilibrationFraction()*m_bucketSize) << "    " << j*m_bucketSize<< endl;
    }
    cout << "Printed ob density! " << endl;
}
void System::setHistogram()
{
    vector<int> histogram(getBins());
    m_histogram = histogram;
}

double System::gradientDescent(double initialAlpha, string filename, int maxIterations){

    int steepestDescentSteps = (int) 1e+5;
    double alpha = initialAlpha;
    double beta = getWaveFunction()->getParameters()[2] / getWaveFunction()->getParameters()[0];
    double lambda = -0.001;
    int iterations = 0;
    double energyDerivative = 100;
    double cumulativeAlpha = 0;
    double tol = 1e-10;
    double percentAlphasToSave = 0.3;
    ofstream myFile;
    myFile.open(filename);

    while (iterations < maxIterations && fabs(energyDerivative) > tol){
        vector<double> parameters(3);
        parameters[0] = alpha;
        parameters[1] = alpha;
        parameters[2] = alpha*beta;
        getWaveFunction()->setParameters(parameters);
        runMetropolisSteps(steepestDescentSteps);
        printOut();
        energyDerivative = findEnergyDerivative();

        // Make sure we accept enough moves (with interaction can get stuck)
        if ((double)m_sampler->getAcceptedNumber() / steepestDescentSteps > 0.90){
            alpha += lambda*energyDerivative;
            iterations++;
        }

        cout << " New alpha = "  << alpha <<  endl;
        cout << " Energy derivative = " << energyDerivative << endl;
        cout << " Iterations = " << iterations << endl;

        // Write alpha, mean local energy and st dev to file
        myFile << alpha << "   "  << getSampler()->getEnergy() << "  " <<
                  sqrt(getSampler()->getCumulativeEnergySquared() - getSampler()->getEnergy()*getSampler()->getEnergy()) << endl;

        if ((double) iterations / maxIterations > 1-percentAlphasToSave){
            cumulativeAlpha += alpha;
        }
    }
    myFile.close();

    alpha = cumulativeAlpha / (maxIterations*percentAlphasToSave);
    return alpha;
}


int System::getBins() const
{
    return m_bins;
}

void System::setBins(int bins)
{
    m_bins = bins;
}

double System::getBucketSize() const
{
    return m_bucketSize;
}

void System::setBucketSize(double bucketSize)
{
    m_bucketSize = bucketSize;
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


double System::findEnergyDerivative()
{

    double meanEnergy      = getSampler()->getCumulativeEnergy() / (m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderiv     = getSampler()->getCumulativeWFderiv() / (m_numberOfMetropolisSteps*getEquilibrationFraction());
    double meanWFderivEloc =  getSampler()->getCumulativeWFderivMultEloc() / (m_numberOfMetropolisSteps*getEquilibrationFraction());




    return 2 * (meanWFderivEloc - meanEnergy*meanWFderiv);
}

