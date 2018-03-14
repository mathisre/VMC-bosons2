#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H
#include "WaveFunctions/wavefunction.h"


class conjugateGradient{

public:
    conjugateGradient(class System* system, double alphaZero, double beta, int CJsteps);
    void conjugateGradientSteps();
    double findEnergyDerivative();

    double alpha() const;
    void setAlpha(double alpha);

    int getCJsteps() const;
    void setCJsteps(int value);

    std::vector<double> getCJparameters() const;
    void setCJparameters(const std::vector<double> &CJparameters);

private:
    double m_alpha = 0;
    int m_CJsteps = 0;
    int m_iterations = 0;

    class System* m_system = nullptr;

    std::vector<double> m_CJparameters = std::vector<double>();
};

#endif // CONJUGATEGRADIENT_H
