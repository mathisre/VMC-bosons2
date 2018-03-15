#include <iostream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include <string>
#include "mpi.h"
#include <fstream>
using std::cout;
using std::endl;
std::ofstream ofile;


Sampler::Sampler(System* system) {
    m_system = system;
    m_stepNumber = 0;
}

void Sampler::setNumberOfMetropolisSteps(int steps) {
    m_numberOfMetropolisSteps = steps;
}

void Sampler::sample(bool acceptedStep) {
    // Make sure the sampling variable(s) are initialized at the first step.
    if (m_stepNumber == 0) {
        m_cumulativeEnergy = 0;
    }

    /* Here you should sample all the interesting things you want to measure.
     * Note that there are (way) more than the single one here currently.
     */

        if (acceptedStep==true){
            m_acceptedNumber++;
            m_energy = m_system->getHamiltonian()->
                       computeLocalEnergy(m_system->getParticles());
            for (int i = 0; i < m_system->getNumberOfDimensions(); i++){
                for (int d = 0; d < m_system->getNumberOfParticles(); d++){
                    m_WFderiv -= m_system->getParticles().at(i).getPosition()[d] * m_system->getParticles().at(i).getPosition()[d];
                    //Remember to include (1,1,beta) vector
                }
            }
            m_WFderivMultELoc = m_WFderiv * m_energy;
        }
        if (m_energy == 0) cout << m_energy << endl;

    if ((double)getStepNumber()/getNumberOfMetropolisSteps() >= 1.0 - m_system->getEquilibrationFraction()){
        m_cumulativeEnergy          += m_energy;
        m_cumulativeEnergySquared   += m_energy*m_energy;
        m_cumulativeWFderiv         += m_WFderiv;
        m_cumulativeWFderivMultEloc += m_WFderivMultELoc;
    }
    m_stepNumber++;
}

void Sampler::printOutputToTerminal(int numprocs) {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps()*numprocs;
    int     p  = m_system->getWaveFunction()->getNumberOfParameters();
    double  ef = m_system->getEquilibrationFraction();
    std::vector<double> pa = m_system->getWaveFunction()->getParameters();
    ofile.close();

    cout << endl;
    cout << "  -- System info -- " << endl;
    cout << " Number of particles  : " << np << endl;
    cout << " Number of dimensions : " << nd << endl;
    cout << " Number of Metropolis steps run : 10^" << std::log10(ms) << endl;
    cout << " Number of equilibration steps  : 10^" << std::log10(std::round(ms*ef)) << endl;
    cout << endl;
    cout << "  -- Wave function parameters -- " << endl;
    cout << " Number of parameters : " << p << endl;
    for (int i=0; i < p; i++) {
        cout << " Parameter " << i+1 << " : " << pa.at(i) << endl;
    }
    cout << endl;
    cout << "  -- Results -- " << endl;
    cout << " Energy : " << m_energy << endl;
    cout << " St. dev: " << sqrt(m_cumulativeEnergySquared - m_energy*m_energy) << endl;
    cout << " Number of accepted steps: " << m_acceptedNumber << endl;
    cout << endl;
}



void Sampler::computeAverages(int myrank, int numprocs)
{
    /* Compute the averages of the sampled quantities. You need to think
     * thoroughly through what is written here currently; is this correct?
     *//*
//    m_energy = m_cumulativeEnergy /(m_system->getNumberOfMetropolisSteps()*(1.0-m_system->getEquilibrationFraction()));
//    m_cumulativeEnergySquared /=(m_system->getNumberOfMetropolisSteps()*(1.0-m_system->getEquilibrationFraction()));*/

    MPI_Reduce(&m_cumulativeEnergy, &m_totalCumulativeEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_cumulativeEnergySquared, &m_totalCumulativeEnergySqr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_acceptedNumber, &m_totalAcceptedNumber, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (myrank == 0){
        m_energy                  = m_totalCumulativeEnergy    / (m_system->getNumberOfMetropolisSteps()*m_system->getEquilibrationFraction()*numprocs);
        m_cumulativeEnergySquared = m_totalCumulativeEnergySqr / (m_system->getNumberOfMetropolisSteps()*m_system->getEquilibrationFraction()*numprocs);
        m_acceptedNumber = m_totalAcceptedNumber;
    }
}


void Sampler::openDataFile(std::string filename){
    if (filename != "0") ofile.open(filename);

}


void Sampler::writeToFile(){
    if (ofile.is_open()) ofile << m_energy << endl;
}

int Sampler::getStepNumber() const
{
    return m_stepNumber;
}

int Sampler::getNumberOfMetropolisSteps() const
{
    return m_numberOfMetropolisSteps;
}

double Sampler::getWFderivMultELoc() const
{
    return m_WFderivMultELoc;
}

void Sampler::setWFderivMultELoc(double WFderivMultELoc)
{
    m_WFderivMultELoc = WFderivMultELoc;
}

double Sampler::getCumulativeWF() const
{
    return m_cumulativeWF;
}

void Sampler::setCumulativeWF(double cumulativeWF)
{
    m_cumulativeWF = cumulativeWF;
}

double Sampler::getWFderiv() const
{
    return m_WFderiv;
}

void Sampler::setWFderiv(double WFderiv)
{
    m_WFderiv = WFderiv;
}

double Sampler::getCumulativeWFderiv() const
{
    return m_cumulativeWFderiv;
}

void Sampler::setCumulativeWFderiv(double cumulativeWFderiv)
{
    m_cumulativeWFderiv = cumulativeWFderiv;
}

double Sampler::getCumulativeWFderivMultEloc() const
{
    return m_cumulativeWFderivMultEloc;
}

void Sampler::setCumulativeWFderivMultEloc(double cumulativeWFderivMultEloc)
{
    m_cumulativeWFderivMultEloc = cumulativeWFderivMultEloc;
}

int Sampler::getAcceptedNumber() const
{
    return m_acceptedNumber;
}

void Sampler::setAcceptedNumber(int acceptedNumber)
{
    m_acceptedNumber = acceptedNumber;
}

void Sampler::setStepNumber(int stepNumber)
{
    m_stepNumber = stepNumber;
}

double Sampler::getCumulativeEnergy() const
{
    return m_cumulativeEnergy;
}

void Sampler::setCumulativeEnergy(double cumulativeEnergy)
{
    m_cumulativeEnergy = cumulativeEnergy;
}
