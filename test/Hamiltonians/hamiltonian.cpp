#include "hamiltonian.h"
#include "WaveFunctions/simplegaussian.h"
#include "../system.h"
#include "particle.h"
#include "vector"
#include <iostream>

using namespace std;

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeNumericalDoubleDerivative (std::vector<class Particle*> particles){

    double h=0.001;
    double h_squared=h*h;
    double wf=0;
    double dimension=m_system->getNumberOfDimensions()*m_system->getNumberOfParticles();

    double backward=0;
    double forward =0;
    double present=0;

    std::vector<Particle*> particles_backward(dimension); //write in a more fancy way??
    std::vector<Particle*> particles_forward(dimension);

    for(int j=0; j < m_system->getNumberOfParticles(); j++){
        std::vector <double> r_backward(m_system->getNumberOfDimensions());
        std::vector <double> r_forward(m_system->getNumberOfDimensions());
        std::vector <double> r(m_system->getNumberOfDimensions());

        for(int d=0; d < m_system->getNumberOfDimensions(); d++){
            r[d]=particles.at(j)->getPosition()[d];
            r_backward[d] = r[d]-h;
            r_forward[d] = r[d]+h;

        }

      //  particles_backward.at(j)->setPosition(r_backward);
        //  particles_forward.at(j)->setPosition(r_forward);

        present= m_system->getWaveFunction()->evaluate(particles);

        particles.at(j)->setPosition(r_backward);

        backward = m_system->getWaveFunction()->evaluate(particles);

        particles.at(j)->setPosition(r_forward);

        forward=m_system->getWaveFunction()->evaluate(particles);

        wf += ( backward+forward-2*present )/h_squared;

        particles.at(j)->setPosition(r);

    }
    return wf;
}
