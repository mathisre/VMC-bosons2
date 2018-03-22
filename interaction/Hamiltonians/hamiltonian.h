#pragma once
#include <vector>

class Hamiltonian {
public:
    Hamiltonian(class System* system);
    virtual double computeLocalEnergy(std::vector<class Particle>& particles) = 0;
    virtual double computeNumericalDoubleDerivative (std::vector<class Particle>& particles);
    virtual double LocalEnergySingleParticle(std::vector<class Particle>& particles, int singParticle) = 0;
    double computeNumericalDoubleDerivativeSingleParticle(std::vector<class Particle>& particles, int singParticle);

protected:
    class System* m_system = nullptr;
};

