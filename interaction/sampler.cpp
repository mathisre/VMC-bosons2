#include <iostream>
#include <cmath>
#include <vector>
#include "sampler.h"
#include "system.h"
#include "particle.h"
#include "Hamiltonians/hamiltonian.h"
#include "WaveFunctions/wavefunction.h"
#include <string>
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

    if (m_stepNumber == 0) {
        m_cumulativeEnergy          = 0;
        m_cumulativeEnergySquared   = 0;
        m_cumulativeWFderiv         = 0;
        m_cumulativeWFderivMultEloc = 0;
    }
    if (acceptedStep==true){
        m_energy=m_system->getHamiltonian()->computeLocalEnergy(m_system->getParticles());
        // Sampling of energy moved to metropolisstep
        m_acceptedNumber++;
        m_WFderiv = 0;
        for (int i = 0; i < m_system->getNumberOfParticles(); i++){
            for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
                m_WFderiv -= m_system->getParticles().at(i).getPosition()[d] * m_system->getParticles().at(i).getPosition()[d];
                //Remember to include (1,1,beta) vector
            }
        }
        m_WFderivMultELoc = m_WFderiv * m_energy;
    }

    //cout<<m_energy<<endl;
     if (((double)getStepNumber()/getNumberOfMetropolisSteps() > 1.0 - m_system->getEquilibrationFraction())||fabs((double)getStepNumber()/getNumberOfMetropolisSteps() -( 1.0 - m_system->getEquilibrationFraction()))<1e-10){

        m_cumulativeEnergy          += m_energy;
        m_cumulativeEnergySquared   += m_energy*m_energy;
        m_cumulativeWFderiv         += m_WFderiv;
        m_cumulativeWFderivMultEloc += m_WFderivMultELoc;
}
    //cout<<m_cumulativeEnergy<<endl;
    m_stepNumber++;
}

void Sampler::printOutputToTerminal() {
    int     np = m_system->getNumberOfParticles();
    int     nd = m_system->getNumberOfDimensions();
    int     ms = m_system->getNumberOfMetropolisSteps();
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
    cout << " Number of accepted steps: " << m_acceptedNumber/m_stepNumber << endl;
    cout << endl;
}



void Sampler::computeAverages() {


    m_energy = m_cumulativeEnergy / (m_system->getNumberOfMetropolisSteps()*m_system->getEquilibrationFraction());
    m_cumulativeEnergySquared /= m_system->getNumberOfMetropolisSteps()*m_system->getEquilibrationFraction();
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

void Sampler::setEnergy(double energy)
{
    m_energy = energy;
}

void Sampler::updateEnergy(double dE){
    m_energy += dE;
}

double Sampler::getCumulativeEnergy() const
{
    return m_cumulativeEnergy;
}
