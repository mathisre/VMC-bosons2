#pragma once
#include "wavefunction.h"
#include "vector"

class SimpleGaussian : public WaveFunction {
public:
    SimpleGaussian(class System* system, double alpha, double beta);
    double evaluate(std::vector<Particle>& particles);
<<<<<<< HEAD
    double computeDoubleDerivative(std::vector<class Particle>& particles);    
    std::vector<double> QuantumForce(std::vector<class Particle>& particles);

=======

//    double QuantumForce(std::vector<class Particle>& particles);


    double computeDoubleDerivative(std::vector<class Particle>& particles);
    std::vector<double> QuantumForce(std::vector<class Particle>& particles);
>>>>>>> e3866af46ea3d05141cab421ed87bbd892d8cac1

protected:
    class WaveFunction* m_wavefunction = nullptr;
};
