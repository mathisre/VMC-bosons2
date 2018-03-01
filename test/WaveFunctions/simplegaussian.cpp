#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include <vector>
#include "wavefunction.h"
#include "../system.h"
#include "../particle.h"
#include <iostream>

using namespace std;
using std::vector;

SimpleGaussian::SimpleGaussian(System* system, double alpha, double beta) :
        WaveFunction(system) {
    assert(alpha >= 0);
    assert(beta >= 0);
    m_numberOfParameters = 3;
    m_parameters.reserve(3);
    m_parameters.push_back(alpha);
    m_parameters.push_back(alpha);
    m_parameters.push_back(alpha*beta);
}

double SimpleGaussian::evaluate(std::vector<class Particle> &particles) {
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


    for(int i=0;i<m_system->getNumberOfParticles();i++){

        for(int j=0;j<m_system->getNumberOfDimensions();j++){
            temp+=particles.at(i).getPosition()[j]*particles.at(i).getPosition()[j];
               int k=0;
            while (i < k && k < m_system->getNumberOfParticles()){
                temp2 += (particles.at(i).getPosition()[j] - particles.at(k).getPosition()[j])
                       * (particles.at(i).getPosition()[j] - particles.at(k).getPosition()[j]);
                k++;
            }
        }
        double r_abs=sqrt(temp);
        if(r_abs <= m_system->getTrapSize()) f=0;
        else f=1-m_system->getTrapSize()/(r_abs);
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
    return exp(-m_parameters[0]*r_squared+u);
}

double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle>& particles) {
    /* All wave functions need to implement this function, so you need to
     * find the double derivative analytically. Note that by double derivative,
     * we actually mean the sum of the Laplacians with respect to the
     * coordinates of each particle.
     *
     * This quantity is needed to compute the (local) energy (consider the
     * Schrödinger equation to see how the two are related).
     */

    /* Thoughts:
     * computeddistanceABS is for all dimensions yet we use the distance between
     * each particle in each dimension.
     *
     * Precalculate matrix R(k,i) that has all distances between particles?
     * Some way to prevent calculating both R_ki and R_ik
     *
     * If r_ij <= a then u' and u'' need to be set to zero.
     * Only the phi is dependent on the number of dimensions. The u's are just dependent on the distance between
     *
     */

    double first,second,third, fourth=0;
    double a=0;
    double temp=0;
    double temp2=0;
    double temp3=0;
    double temp4= 0;
    double temp5 = 0;
    double temp6 = 0;
    double R_kj = 0;
    double R_ki = 0;

    for(int k=0; k<m_system->getNumberOfParticles();k++){
        for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
            first += (4*particles.at(k).getPosition()[d]*particles.at(k).getPosition()[d]*m_parameters[d]*m_parameters[d])
                    -m_parameters[d]*2;
                    //-m_parameters[0]*(2*beta_vector[d]);
            temp5 -= m_parameters[d]*particles.at(k).getPosition()[d];
        }
        //first -= 2*(m_system->getNumberOfDimensions())*m_parameters[0];// 2*m_parameters[2];

        for(int j=0; j<m_system->getNumberOfParticles(); j++){
            if (j != k){
                R_kj = m_system->computedistanceABS(k,j);
                temp4 = a / (R_kj * (R_kj - a));
                temp2 -= temp4*temp4;

                for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
                    temp6 -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);
                    //temp -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);

                }

                temp += 4 * temp6 * temp4;
                //temp *= 4*temp4;
                for(int i=0; i<m_system->getNumberOfParticles(); i++){
                    if(i != k){
                        R_ki = m_system->computedistanceABS(k,i);
                        temp3 += temp4 * a / (R_ki * (R_ki  - a));
                        //temp3 += a*a / (R_kj * R_ki * (R_kj  - a) * (R_ki -a));
                    }
                }
            }





            /*
             *
             *
             *
            for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
                if (j != k){
                    R_kj = m_system->computedistanceABS(k,j);
                    temp4 = a / (R_kj * (R_kj - a));
                    temp -= 4 * m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]) * temp4;
                    temp2 -= temp4*temp4;
                }
                for(int i=0; i<m_system->getNumberOfParticles(); i++){
                    if(j!=k && i!=k){
                        R_ki = m_system->computedistanceABS(k,i);
                        temp3 += a*a / (R_kj * R_ki * (R_kj  - a) * (R_ki -a));
                    }
                }
            } */
        }
    }

    second = temp;
    third  = temp3;
    fourth = temp2;

    return first+second+third+fourth;
}


std::vector<double> SimpleGaussian::QuantumForce(std::vector<class Particle>& particles) {
    double u_p = 0;
    double first = 0;
    double a = 0;
    double R_kj;

    std::vector<double> QuantumForce = std::vector<double>();


    for (int k = 0; k < m_system->getNumberOfParticles(); k++){
        for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
            first -= m_parameters[d]*particles.at(k).getPosition()[d];
        }
        QuantumForce.at(k) += first;
        for (int j = 0; j < m_system->getNumberOfParticles(); j++){
            if (j != k){
                R_kj = m_system->computedistanceABS(k,j);
                QuantumForce += (particles.at(k).getPosition() - particles.at(j).getPosition()) *
                        a / (R_kj*R_kj*(R_kj-a));


            }

        }
    }
    first *= 2;
    return QuantumForce;
}
