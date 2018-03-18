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

    double h=0.0001;
    double h_squared=h*h;
    double wf=0;
    double dimension=m_system->getNumberOfDimensions()*m_system->getNumberOfParticles();

    double backward=0;
    double forward =0;
    double present=0;

//    std::vector<Particle> particles_backward = particles;
//    particles_backward.reserve(m_system->getNumberOfDimensions());//write in a more fancy way??
//    std::vector<Particle> particles_forward;
//    particles_forward.reserve(m_system->getNumberOfDimensions());

 std::vector <double> r(m_system->getNumberOfDimensions());

    for(int j=0; j < m_system->getNumberOfParticles(); j++){
       // std::vector <double> r_backward(m_system->getNumberOfDimensions());
        //std::vector <double> r_forward(m_system->getNumberOfDimensions());


        for(int d=0; d < m_system->getNumberOfDimensions(); d++){
            r[d]=particles.at(j).getPosition()[d];
            r[d]-=h;
            //r_backward[d] = r[d]-h;
            //r_forward[d] = r[d]+h;

        }
        //particles_backward.at(j).setPosition(r_backward);
        //particles_forward.at(j).setPosition(r_forward);
        particles.at(j).setPosition(r);


        backward = m_system->getWaveFunction()->evaluate(particles);
        for(int d=0; d < m_system->getNumberOfDimensions(); d++){
            //r[d]=particles.at(j).getPosition()[d];
            r[d]+=2*h;
            //r_backward[d] = r[d]-h;
            //r_forward[d] = r[d]+h;

        }

        particles.at(j).setPosition(r);
        forward = m_system->getWaveFunction()->evaluate(particles);

        for(int d=0; d < m_system->getNumberOfDimensions(); d++){
            //r[d]=particles.at(j).getPosition()[d];
            r[d]-=h;
            //r_backward[d] = r[d]-h;
            //r_forward[d] = r[d]+h;

        }

        particles.at(j).setPosition(r);
        present = m_system->getWaveFunction()->evaluate(particles);
    wf += ( backward+forward-2*present )/h_squared;
}

return wf/present;
}


/*

 backward = m_system->getWaveFunction()->evaluate(particles);

    for(int j=0; j < m_system->getNumberOfParticles(); j++){
        for(int d=0; d<m_system->getNumberOfDimensions(); d++){
            r[d]+=2*h;
            //cout<<"2"<<r[d]<<endl;
        }
        particles.at(j).setPosition(r);
    }

forward=m_system->getWaveFunction()->evaluate(particles);

    for(int j=0; j < m_system->getNumberOfParticles(); j++){
        for(int d=0; d<m_system->getNumberOfDimensions(); d++){
            r[d]-=h;
           // cout<<"3"<<r[d]<<endl;
        }
        particles.at(j).setPosition(r);
    }

present= m_system->getWaveFunction()->evaluate(particles);

        //present= m_system->getWaveFunction()->evaluate(particles);

        //particles.at(j)->setPosition(r_backward);

       // backward = m_system->getWaveFunction()->evaluate(particles_backward);

      //  particles.at(j)->setPosition(r_forward);

        //forward=m_system->getWaveFunction()->evaluate(particles_forward);

        wf = ( backward+forward-2*present ) /h_squared;

        cout<<(backward+forward-2*present) /h_squared<<endl;
        //particles.at(j).setPosition(r);

    return wf/present;
}

*/
