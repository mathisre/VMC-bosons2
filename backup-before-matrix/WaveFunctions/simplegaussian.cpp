#include "simplegaussian.h"
#include <cmath>
#include <cassert>
#include <vector>
#include <algorithm>
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
    double r_squared = 0;
    double temp2=0;
    double u=0;
    double f=0;
    double trap_size_sqr = m_system->getTrapSize()*m_system->getTrapSize();


    for(int i=0;i < m_system->getNumberOfParticles();i++){

        for(int d=0; d < m_system->getNumberOfDimensions();d++){
            r_squared += particles.at(i).getPosition()[d]*particles.at(i).getPosition()[d];
            int k=0;
            while (i < k && k < m_system->getNumberOfParticles()){
                temp2 +=  (particles.at(i).getPosition()[d] - particles.at(k).getPosition()[d])
                        * (particles.at(i).getPosition()[d] - particles.at(k).getPosition()[d]);
                k++;
            }
        }/*
        if(r_squared <= trap_size_sqr) f = 0;
        else f = 1 - m_system->getTrapSize() / (sqrt(r_squared));
        u += log(f);
*/
        /*
        double r_abs=sqrt(temp);
        if(r_abs <= m_system->getTrapSize()) f = 0;
        else f = 1 - m_system->getTrapSize() / (r_abs);
        u+=log(f);
        */
    }


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
    //    double temp=0;
    //    double temp2=0;
    //    double temp3=0;
    double temp4= 0;
    //double temp5 = 0;
    double temp6 = 0;
    double R_kj = 0;
    double R_ki = 0;
    vector <double> m_parameters_squared(3);

    transform(m_parameters.begin(),m_parameters.end(),m_parameters.begin(),  std::bind1st (std::multiplies <double> () , 2.0));
    transform(m_parameters.begin(), m_parameters.end(), m_parameters.begin(),m_parameters_squared.begin(), multiplies<double>());
    transform(m_parameters.begin(),m_parameters.end(),m_parameters.begin(),  std::bind1st (std::multiplies <double> () , 1/2.0));

    for(int k=0; k<m_system->getNumberOfParticles();k++){
        for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
            first += (m_parameters_squared[d]*particles.at(k).getPosition()[d]*particles.at(k).getPosition()[d])
                    -m_parameters_squared[d]; //*2
            //-m_parameters[0]*(2*beta_vector[d]);
            //temp5 -= m_parameters[d]*particles.at(k).getPosition()[d];
        }
    }
    if(a==0) return first;
    int i=0;
    int j=0;
    while(i<m_system->getNumberOfParticles()){
        while(j<m_system->getNumberOfParticles()){
            for(int k=0; k<i;k++){
                //first -= 2*(m_system->getNumberOfDimensions())*m_parameters[0];// 2*m_parameters[2];
                //if(a==0) second =0, fourth = 0, third = 0;
                //else {

                R_kj = m_system->computedistanceABS(k,j);
                temp4 = a / (R_kj * (R_kj - a));
                fourth -= temp4*temp4;

                for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
                    temp6 -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);
                    //temp -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);

                }

                second += temp6 * temp4;
                //temp *= 4*temp4;

                R_ki = m_system->computedistanceABS(k,i);
                third += temp4 * a / (R_ki * (R_ki  - a));
                //temp3 += a*a / (R_kj * R_ki * (R_kj  - a) * (R_ki -a));


            }
            j++;
        }
        i++;
    }



    /*
    vector<vector<double>>distancematrix=m_system->computematrixdistance(particles);
    for(int k=0; k<m_system->getNumberOfParticles();k++){
        for(int j=0; j<m_system->getNumberOfParticles(); j++){
            if (j != k){
                //R_kj = m_system->computedistanceABS(k,j);
                temp4 = a / (distancematrix[k][j] * (distancematrix[k][j] - a));
                fourth -= temp4*temp4;

                for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
                    temp6 -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);
                    //temp -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);

                }

                second += temp6 * temp4;
                //temp *= 4*temp4;
                for(int i=0; i<m_system->getNumberOfParticles(); i++){
                    if(i != k){
                        //R_ki = m_system->computedistanceABS(k,i);
                        third += temp4 * a / (distancematrix[k][i] * (distancematrix[k][i]  - a));
                        //temp3 += a*a / (R_kj * R_ki * (R_kj  - a) * (R_ki -a));
                    }
                }
            }
        }
    }
    */



   return first+4*second*2+third*4+fourth*2;
}

    /*
     for(int k=0; k<j;k++){
        //first -= 2*(m_system->getNumberOfDimensions())*m_parameters[0];// 2*m_parameters[2];
<<<<<<< HEAD
        //if(a==0) second =0, fourth = 0, third = 0;
        //else {
=======


        if(a==0) second =0, fourth = 0, third = 0;
        else {
>>>>>>> e3866af46ea3d05141cab421ed87bbd892d8cac1
            for(int j=0; j<m_system->getNumberOfParticles(); j++){
                if (j != k){
                    R_kj = m_system->computedistanceABS(k,j);
                    temp4 = a / (R_kj * (R_kj - a));
                    fourth -= temp4*temp4;

                    for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
                        temp6 -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);
                        //temp -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);

                    }

                    second += temp6 * temp4;
                    //temp *= 4*temp4;
                    for(int i=0; i<m_system->getNumberOfParticles(); i++){
                        if(i != k){
                            R_ki = m_system->computedistanceABS(k,i);
                            third += temp4 * a / (R_ki * (R_ki  - a));
                            //temp3 += a*a / (R_kj * R_ki * (R_kj  - a) * (R_ki -a));
                        }
                    }
                }



            }
    }
    return first+4*second+third+fourth;
}
*/
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

       // } //else loop

    //    second = temp;
    //    third  = temp3;
    //    fourth = temp2;




std::vector<double> SimpleGaussian::QuantumForce(std::vector<class Particle>& particles) {

    double a = 0;
    double R_kj;
    double constant;

    std::vector<double> QuantumForce(3);
    std::vector<double> u_deriv(3);
    for (int k = 0; k < m_system->getNumberOfParticles(); k++){

        for (int j = 0; j < k; j++){
                R_kj = m_system->computedistanceABS(k,j);
                constant = a / (R_kj*R_kj*(R_kj-a));
                for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
                    u_deriv[d] += (particles.at(k).getPosition()[d] - particles.at(j).getPosition()[d]) * constant;
                }
            }

        for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
            QuantumForce[d] += -2 * (m_parameters[d]*particles.at(k).getPosition()[d])
                               + u_deriv[d];
        }
    }



    return QuantumForce;

}
