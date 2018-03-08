#include "conjugategradient.h"
#include "WaveFunctions/simplegaussian.h"
#include "WaveFunctions/wavefunction.h"
#include "Hamiltonians/harmonicoscillator.h"
#include "system.h"
#include "particle.h"
#include "sampler.h"

conjugateGradient::conjugateGradient(System *system, double alphaZero, int CJsteps)
{
    setAlpha(alphaZero);
    setCJsteps(CJsteps);
}


void conjugateGradient::conjugateGradientSteps()
{
    double energyDerivative, meanEnergy;
    double iterations = 0;
    double tol = (double) 1e-10;
    energyDerivative =findEnergyDerivative();
    while (energyDerivative > tol)
    {
        m_system->runMetropolisSteps(m_CJsteps);
        meanEnergy = m_system->m_sampler->m_cumulativeEnergy / getCJsteps();
        energyDerivative = findEnergyDerivative();
        // Find new alpha
        iterations++;
    }

}

double conjugateGradient::findEnergyDerivative()
{
    // Are we supposed to do analytical or numerical here?
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
