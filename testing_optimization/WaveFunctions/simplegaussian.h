#pragma once
#include "wavefunction.h"
#include "vector"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha, double beta);
    double evaluate(std::vector<Particle>& particles);

    double computeDoubleDerivative(std::vector<class Particle>& particles);
    double computeDoubleDerivativeSingleParticle  (std::vector<class Particle>& particles, int singParticle);
    std::vector<double> QuantumForce              (std::vector<class Particle>& particles);
    std::vector<double> QuantumForceSingleParticle(std::vector<class Particle>& particles, int singParticle);

protected:
    class WaveFunction* m_wavefunction = nullptr;
};
