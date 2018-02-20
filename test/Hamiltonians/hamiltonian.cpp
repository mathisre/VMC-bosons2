#include "hamiltonian.h"
#include "WaveFunctions/simplegaussian.h"
#include "../system.h"
#include "particle.h"
#include "vector"

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeNumericalDoubleDerivative (std::vector<Particle*> particles){

    double h=0.01;
    double h_squared=h*h;
    double wf=0;

    std::vector<Particle*> particles_backward;
    std::vector<Particle*> particles_forward;
    std::vector <double> r_backward(m_system->getNumberOfDimensions());
    std::vector <double> r_forward(m_system->getNumberOfDimensions());

    for(int j=0; j < m_system->getNumberOfParticles(); j++){

        for(int d=0; d < m_system->getNumberOfDimensions(); d++){

            r_backward[d] = particles.at(j)->getPosition()[d]-h;

            r_forward[d] = particles.at(j)->getPosition()[d]+h;

        }

        particles_backward.at(j)->setPosition(r_backward);

        particles_forward.at(j)->setPosition(r_forward);

        wf -= ( m_system->getWaveFunction()->evaluate(particles_backward)+m_system->getWaveFunction()->evaluate(particles_forward)
             -2*m_system->getWaveFunction()->evaluate(particles) )/h_squared;

    }

    return wf;
}
