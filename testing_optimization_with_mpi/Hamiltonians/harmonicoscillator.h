#pragma once
#include "hamiltonian.h"
#include <vector>

class HarmonicOscillator : public Hamiltonian {
public:
    HarmonicOscillator(System* system, double omega, double omega_z);
    double computeLocalEnergy(std::vector<Particle>& particles);
    double LocalEnergySingleParticle(std::vector<Particle>& particles, int singParticle);

private:
    std::vector<double> m_omega = std::vector<double>();
};

