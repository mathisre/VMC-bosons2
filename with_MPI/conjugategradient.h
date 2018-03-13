#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H
#include "WaveFunctions/wavefunction.h"


class conjugateGradient{

public:
    conjugateGradient(class System* system, double alphaZero, int CJsteps);
    void conjugateGradientSteps();
    double findEnergyDerivative();

    double alpha() const;
    void setAlpha(double alpha);

    int getCJsteps() const;
    void setCJsteps(int value);

private:
    double m_alpha = 0;
    int m_CJsteps = 0;

    class System* m_system = nullptr;
};

#endif // CONJUGATEGRADIENT_H
