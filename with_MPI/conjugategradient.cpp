#include "conjugategradient.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "system.h"
#include "particle.h"
#include "sampler.h"

conjugateGradient::conjugateGradient(System *system, double alphaZero, double beta, int CJsteps)
{
    setAlpha(alphaZero);
    m_CJparameters.reserve(3);
    m_CJparameters.push_back(alphaZero);
    m_CJparameters.push_back(alphaZero);
    m_CJparameters.push_back(alphaZero*beta);
    setCJsteps(CJsteps);
}


void conjugateGradient::conjugateGradientSteps()
{
    double energyDerivative, meanEnergy;
    m_iterations = 0;
    double tol = (double) 1e-10;
    m_system->getWaveFunction()->setParameters(m_CJparameters);
    m_system->runMetropolisSteps(m_CJsteps);
    energyDerivative =findEnergyDerivative();
    while (energyDerivative > tol){
        m_system->runMetropolisSteps(m_CJsteps);

        energyDerivative = findEnergyDerivative();
        // Find new alpha

        m_system->getWaveFunction()->setParameters(m_CJparameters);
        m_iterations++;
    }

}

double conjugateGradient::findEnergyDerivative()
{
    double meanEnergy      = m_system->getSampler()->getCumulativeEnergy() / getCJsteps();
    double meanWFderiv     = m_system->getSampler()->getCumulativeWFderiv() / getCJsteps();
    double meanWFderivEloc = m_system->getSampler()->getCumulativeWFderivMultEloc() / getCJsteps();

    return 2 * (meanWFderivEloc - meanEnergy*meanWFderiv);
}


double conjugateGradient::alpha() const
{
    return m_alpha;
}

void conjugateGradient::setAlpha(double alpha)
{
    m_alpha = alpha;
}

int conjugateGradient::getCJsteps() const
{
    return m_CJsteps;
}

void conjugateGradient::setCJsteps(int value)
{
    m_CJsteps = value;
}

std::vector<double> conjugateGradient::getCJparameters() const
{
    return m_CJparameters;
}

void conjugateGradient::setCJparameters(const std::vector<double> &CJparameters)
{
    m_CJparameters = CJparameters;
}
