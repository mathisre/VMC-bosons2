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
    /* Here, you need to compute the kinetic and potential energies. Note that
     * when using numerical differentiation, the computation of the kinetic
     * energy becomes the same for all Hamiltonians, and thus the code for
     * doing this should be moved up to the super-class, Hamiltonian.
     *
     * You may access the wave function currently used through the
     * getWaveFunction method in the m_system object in the super-class, i.e.
     * m_system->getWaveFunction()...
     */

    double potentialEnergy = 0;
    double kineticEnergy   = 0;




    kineticEnergy = -0.5*m_system->getWaveFunction()->computeDoubleDerivative(particles);
    //kineticEnergy=-0.5*m_system->getHamiltonian()->computeNumericalDoubleDerivative(particles);


    //cout<<"kinetic"<<kineticEnergy<<endl;

    for (int k = 0; k < m_system->getNumberOfParticles(); k++ ){
        for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
            potentialEnergy += m_omega[d]*m_omega[d]*particles.at(k).getPosition()[d]*particles.at(k).getPosition()[d];
            //cout<<"---- "<<m_omega[d]<<endl;
        }

//        for (int j = 0; j < m_system->getNumberOfParticles();j++){ // Some better way?
//            potentialEnergy += (int)(1e10) * (m_system->computedistanceABS(k,j) > m_system->getTrapSize());
//        }

    }
    potentialEnergy *= 0.5;
    //2potentialEnergy /= m_system->getWaveFunction()->evaluate(particles);

   // cout<<"pot"<<potentialEnergy;
    //cout << "\t kin" << kineticEnergy<<" \n";
    return kineticEnergy + potentialEnergy;
}




