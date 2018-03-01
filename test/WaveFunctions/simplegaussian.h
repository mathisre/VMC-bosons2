#pragma once
#include "wavefunction.h"
#include "vector"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha, double beta);
    double evaluate(std::vector<Particle>& particles);
<<<<<<< HEAD
    double computeDoubleDerivative(std::vector<class Particle>& particles);
    double QuantumForce(std::vector<class Particle>& particles);
=======
    double computeDoubleDerivative(std::vector<class Particle>& particles);    
    std::vector<double> QuantumForce(std::vector<class Particle>& particles);
>>>>>>> 8110118d3b39098718ea19524c2e3fd5ea5612d5

protected:
    class WaveFunction* m_wavefunction = nullptr;
};
