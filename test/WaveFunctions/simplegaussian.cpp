#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include <vector>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
using std::vector;

SimpleGaussian::SimpleGaussian(System* system, double alpha) :
        WaveFunction(system) {
    assert(alpha >= 0);
    m_numberOfParameters = 1;
    m_parameters.reserve(1);
    m_parameters.push_back(alpha);
}

double SimpleGaussian::evaluate(std::vector<class Particle*> particles) {
    /* You need to implement a Gaussian wave function here. The positions of
     * the particles are accessible through the particle[i].getPosition()
     * function.
     *
     * For the actual expression, use exp(-alpha * r^2), with alpha being the
     * (only) variational parameter.
     */
    double temp=0;
    double temp2=0;
    double u=0;
    double f=0;
    double a=0;
    int k=0;
    for(int i=0;i<m_system->getNumberOfParticles();i++){
        for(int j=0;j<m_system->getNumberOfDimensions();j++){
            temp+=particles.at(i)->getPosition()[j]*particles.at(i)->getPosition()[j];
            while (i<k&&k<m_system->getNumberOfParticles()){
                temp2+=(particles.at(i)->getPosition()[j]-particles.at(k)->getPosition()[j])*(particles.at(i)->getPosition()[j]-particles.at(k)->getPosition()[j]);
                k++;
            }
        }
        double r_abs=sqrt(temp);
        if(r_abs<=a) f=0;
        else f=1-a/(r_abs);
        u+=log(f);
    }

    double r_squared=temp;

//    for(int j=0; j<m_system->getNumberOfParticles();j++){
//        for(int i=0; i<j; i++){
//            if(m_system->computedistanceABS(i,j)) f=0;
//            else f=1-a/(m_system->computedistanceABS(i,j));
//            u+=log(f);
//        }
//    }
    return exp(-m_parameters[1]*r_squared+u);
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle*> particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    double first,second,third, fourth=0;
    double a=0;
    double temp=0;
    double temp2=0;
    double temp3=0;
    for(int k=0; k<m_system->getNumberOfParticles();k++){
        first=m_parameters[1]*m_parameters[1]*(4*particles.at(k)->getPosition()[1]*particles.at(k)->getPosition()[1]+4*particles.at(k)->getPosition()[2]*particles.at(k)->getPosition()[2]+
                4*m_parameters[2]*particles.at(k)->getPosition()[3]*particles.at(k)->getPosition()[3])-m_parameters[1]*(4+2*m_parameters[2]);
        for(int j=0; j<m_system->getNumberOfParticles(); j++){
            if(j!=k) temp+=a/(m_system->computedistanceABS(k,j)*(m_system->computedistanceABS(k,j)-a)), temp2+=-temp*temp;
            for(int i=0; i<m_system->getNumberOfParticles(); i++){
            if(j!=k&&i!=k) temp3+=a*a/(m_system->computedistanceABS(k,j)*m_system->computedistanceABS(k,i)*(m_system->computedistanceABS(k,j)-a)*(m_system->computedistanceABS(k,i)-a));
            }
        }
        second=temp;
        fourth=temp2;
        third=temp3;
    }
    return first+second+third+fourth;
}
