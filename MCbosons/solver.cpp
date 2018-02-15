#include "solver.h"

solver::solver()
{


}

void solver::Metropolis(double rng, double WRatio, )
{
    if (rng <= Wratio){
        for (int k = 0; k<m_Dimensions; k++){
            r_Old(i,k) = r_New(i,k);
            Psi_Old = Psi_New
        }
    }
}

double solver::E_local()
{}

int solver::getNumerOfAtoms() const
{
    return m_NumerOfAtoms;
}

void solver::setNumerOfAtoms(int value)
{
    m_NumerOfAtoms = value;
}

int solver::getDimensions() const
{
    return m_Dimensions;
}

void solver::setDimensions(int value)
{
    m_Dimensions = value;
}

int solver::getNumberOfSteps() const
{
    return m_NumberOfSteps;
}

void solver::setNumberOfSteps(int value)
{
    m_NumberOfSteps = value;
}

double solver::getStepSize() const
{
    return m_StepSize;
}

void solver::setStepSize(double value)
{
    m_StepSize = value;
}
