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

        energyDerivative = findEnergyDerivative();
        // Find new alpha
        iterations++;
    }

}

double conjugateGradient::findEnergyDerivative()
{

    double meanEnergy      = m_system->getSampler()->getEnergy() / getCJsteps();
    double meanWFderiv     = m_system->getSampler()->getWFderiv() / getCJsteps();
    double meanWFderivEloc = m_system->getSampler()->getWFderivMultELoc() / getCJsteps();




    // Make the sampler sample the wavefunc deriv and things like that
    // Then just find the mean and we are good

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
