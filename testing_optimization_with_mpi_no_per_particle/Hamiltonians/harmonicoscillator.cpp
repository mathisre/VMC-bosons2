#include "harmonicoscillator.h"
#include <cassert>
#include <iostream>
#include "../system.h"
#include "../particle.h"
#include "../WaveFunctions/wavefunction.h"

using std::cout;
using std::endl;

HarmonicOscillator::HarmonicOscillator(System* system, double omega, double omega_z) :
        Hamiltonian(system) {
    assert(omega > 0);
    assert(omega_z > 0);
    m_omega.reserve(3);
    m_omega.push_back(omega);
    m_omega.push_back(omega);
    m_omega.push_back(omega_z);

}

double HarmonicOscillator::computeLocalEnergy(std::vector<Particle> &particles) {

    double potentialEnergy = 0;
    double kineticEnergy   = 0;




    kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);
    //kineticEnergy = -0.5*m_system->getHamiltonian()->computeNumericalDoubleDerivative(particles);



    for (int k = 0; k < m_system->getNumberOfParticles(); k++ ){
        for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
            potentialEnergy += m_omega[d]*m_omega[d]*particles.at(k).getPosition()[d]*particles.at(k).getPosition()[d];

        }


    }
    potentialEnergy *= 0.5;

    //cout<<"Kinetic energy = "<<kineticEnergy<<endl;
    return kineticEnergy + potentialEnergy;
}

double HarmonicOscillator::LocalEnergySingleParticle(std::vector<Particle> &particles, int singParticle) {


    double potentialEnergy = 0;
    double kineticEnergy;

   //kineticEnergy = m_system->getWaveFunction()->computeDoubleDerivativeSingleParticle(particles, singParticle); //-0.5*m_system->getWaveFunction()->computeDoubleDerivativeSingleParticle(particles, singParticle);
   // kineticEnergy = -0.5*m_system->getHamiltonian()->computeNumericalDoubleDerivativeSingleParticle(particles, singParticle);



    for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
        potentialEnergy += m_omega[d]*m_omega[d]*particles.at(singParticle).getPosition()[d]*particles.at(singParticle).getPosition()[d];

    }


    potentialEnergy *= 0.5;

    //cout<<potentialEnergy<<"|||||||||"<<kineticEnergy<<endl;
    return kineticEnergy + potentialEnergy;
}




