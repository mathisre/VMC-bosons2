#include "solver.h"
#include <vector>
#include <math.h>
using namespace std;
using std::vector;
solver::solver()
{


}

void solver::Metropolis(double rng, double &Psi_Old, double r_New, vector<double> &r_Old)
{
    double Psi_New = WaveFunction(R_New);

    if (rng <= Psi_New*Psi_New/(Psi_Old*Psi_Old)){
        for (int k = 0; k < m_Dimensions; k++){
            r_Old(i,k) = r_New(i,k);
            Psi_Old = Psi_New;
        }
    }
}

double solver::WaveFunction(vector<double> R, double alpha)
{
    double Psi = 1;
    double g_sum = 0;
    double g = 0;
    double f = 1;
    double r_kj;
    for (int k = 0; k < m_NumerOfParticles; k++){
        g = 0;
        g_sum = 0;
        for (int d = 0; d<m_Dimensions; d++){
            g_sum += R(k,d)*R(k,d); //Add some beta
        }
        g = exp(-alpha * g_sum);
        for (int j = 0; j < k; j++){
            r_kj = sqrt((r(k,0)-r(j,0))*(r(k,0)-r(j,0)) + (r(k,1)-r(j,1))*(r(k,1)-r(j,1)) + (r(k,2)-r(j,2))*(r(k,2)-r(j,2)));
            if (r_kj < m_a){
                f *= 0;
                j = i;
            }
            else {
                f *= 1 - m_a / r_kj;
            }
        }
    }
        Psi *=
}

double solver::E_local()
{}

int solver::getNumerOfParticles() const
{
    return m_NumerOfParticles;
}

void solver::setNumerOfParticles(int value)
{
    m_NumerOfParticles = value;
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

int solver::getNumerOfParticles() const
{
    return m_NumerOfParticles;
}

void solver::setNumerOfParticles(int NumerOfParticles)
{
    m_NumerOfParticles = NumerOfParticles;
}
