#include "hamiltonian.h"
#include "../WaveFunctions/simplegaussian.h"
#include "../system.h"
#include "../particle.h"
#include "vector"
#include <iostream>

using namespace std;

Hamiltonian::Hamiltonian(System* system) {
    m_system = system;
}

double Hamiltonian::computeNumericalDoubleDerivative (std::vector<Particle> &particles){

    double h=0.001;
    double h_squared=h*h;
    double wf=0;

    double backward=0;
    double forward =0;
    double present=0;


    present=m_system->getPsiOld();

 std::vector <double> r(m_system->getNumberOfDimensions());
 for(int j=0; j < m_system->getNumberOfParticles(); j++){

 for(int d=0; d<m_system->getNumberOfDimensions(); d++){
     r[d]=particles.at(j).getPosition()[d];
 }
     for(int d=0; d < m_system->getNumberOfDimensions(); d++){

          r[d]-=h;
          particles.at(j).setPosition(r);
          backward += m_system->getWaveFunction()->evaluate(particles);

          r[d]+=2*h;
          particles.at(j).setPosition(r);
          forward += m_system->getWaveFunction()->evaluate(particles);

          r[d]-=h;
          particles.at(j).setPosition(r);
      }
  }

 wf=(backward+forward-2*present*m_system->getNumberOfDimensions()*m_system->getNumberOfParticles())/(h_squared);
 return wf/(present);
}

double Hamiltonian::computeNumericalDoubleDerivativeSingleParticle(std::vector<class Particle>& particles, int singParticle){
    double h=0.001;
    double h_squared=h*h;
    double wf=0;

    double backward=0;
    double forward =0;
    double present=0;


    std::vector <double> r(m_system->getNumberOfDimensions());


    present =m_system->getPsiOld();

    for(int d=0; d<m_system->getNumberOfDimensions(); d++){
        r[d]=particles.at(singParticle).getPosition()[d];
    }

        for(int d=0; d < m_system->getNumberOfDimensions(); d++){


            // r[d]=particles.at(singParticle).getPosition()[d];



             r[d]-=h;
             particles.at(singParticle).setPosition(r);

             backward += m_system->getWaveFunction()->evaluate(particles);

             r[d]+=2*h;
             particles.at(singParticle).setPosition(r);
             forward += m_system->getWaveFunction()->evaluate(particles);

             r[d]-=h;

             particles.at(singParticle).setPosition(r);


         }


    wf=(backward+forward-2*present*m_system->getNumberOfDimensions())/(h_squared);
   return wf/present;
}


