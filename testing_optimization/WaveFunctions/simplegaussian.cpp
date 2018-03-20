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
    m_parameters.push_back(alpha); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    m_parameters.push_back(alpha*beta);
    //m_parameters.push_back(beta);
}

double SimpleGaussian::evaluate(std::vector<class Particle> &particles) {

    // WF is just the product of the individual wavefunctions
    // Divide by the wf of the old particle and multiply in the wf of the new one

    double r_squared = 0;
    double f=1;

    for(int j=0; j<m_system->getNumberOfParticles();j++){

        if(m_system->getNumberOfParticles()==1) break;
        for(int i=0; i<j; i++){
            f = 1-m_system->getTrapSize()/(m_system->getDistanceMatrixij(i,j));
    }
}

for(int i=0;i < m_system->getNumberOfParticles();i++){

        for(int d=0; d < m_system->getNumberOfDimensions();d++){
            r_squared += particles.at(i).getPosition()[d]*particles.at(i).getPosition()[d]*m_parameters[d];
                  }
}

return exp(-r_squared)*f;

//if(m_system->getNumberOfDimensions()==1){
//    for(int i=0;i < m_system->getNumberOfParticles();i++){
//            r_squared += particles.at(i).getPosition()[0]*particles.at(i).getPosition()[0];
//    }

//    return exp(-m_parameters[0]*r_squared+u);

//}

//if(m_system->getNumberOfDimensions()==2){
//    for(int i=0;i < m_system->getNumberOfParticles();i++){
//            r_squared += particles.at(i).getPosition()[0]*particles.at(i).getPosition()[0]+particles.at(i).getPosition()[1]*particles.at(i).getPosition()[1];
//    }
//    return exp(-m_parameters[0]*r_squared+u);

//}

//if(m_system->getNumberOfDimensions()==3){
//    for(int i=0;i < m_system->getNumberOfParticles();i++){
//            r_squared += particles.at(i).getPosition()[0]*particles.at(i).getPosition()[0]+particles.at(i).getPosition()[1]*particles.at(i).getPosition()[1]+
//                    m_parameters[2]*particles.at(i).getPosition()[2]*particles.at(i).getPosition()[2]/m_parameters[0];
//    }
//    return exp(-m_parameters[0]*r_squared)*f;

//}

}



//        for(int j=0; j<m_system->getNumberOfParticles();j++){
//            for(int i=0; i<j; i++){
//                if(m_system->computedistanceABS(i,j)<= m_system->getTrapSize()) f=0;
//                else f=1-m_system->getTrapSize()/(m_system->computedistanceABS(i,j));
//                u+=log(f);
//            }
//        }






//for(int i=0;i < m_system->getNumberOfParticles();i++){

//        for(int d=0; d < m_system->getNumberOfDimensions();d++){
//            r_squared += particles.at(i).getPosition()[d]*particles.at(i).getPosition()[d];
//                  }
//}


//            while (i < k && k < m_system->getNumberOfParticles()){
//                temp2 +=  (particles.at(i).getPosition()[d] - particles.at(k).getPosition()[d])
//                        * (particles.at(i).getPosition()[d] - particles.at(k).getPosition()[d]);
//                k++;
//            }
        /*
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


//        for(int j=0; j<m_system->getNumberOfParticles();j++){
//            for(int i=0; i<j; i++){
//                if(m_system->computedistanceABS(i,j)) f=0;
//                else f=1-a/(m_system->computedistanceABS(i,j));
//                u+=log(f);
//            }
//        }
//    return exp(-m_parameters[0]*r_squared+u);

//}


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

    double first = 0;
    double second = 0, third = 0, fourth=0;
    double a=0;
    //    double temp=0;
    //    double temp2=0;
    //    double temp3=0;
    double temp4= 0;
    //double temp5 = 0;
    double temp6 = 0;
    vector <double> m_parameters_squared(3);

    transform(m_parameters.begin(),m_parameters.end(),m_parameters.begin(),  std::bind1st (std::multiplies <double> () , 2.0));
    transform(m_parameters.begin(), m_parameters.end(), m_parameters.begin(),m_parameters_squared.begin(), multiplies<double>());
    transform(m_parameters.begin(),m_parameters.end(),m_parameters.begin(),  std::bind1st (std::multiplies <double> () , 1/2.0));

    for(int k=0; k<m_system->getNumberOfParticles();k++){
        for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
            first += (m_parameters_squared[d]*particles.at(k).getPosition()[d]*particles.at(k).getPosition()[d])
                    -m_parameters_squared[d];
            //*2
            //-m_parameters[0]*(2*beta_vector[d]);
            //temp5 -= m_parameters[d]*particles.at(k).getPosition()[d];
        }
    }
    if(a == 0) return first;
    int i=0;
    int j=0;
    while(i<m_system->getNumberOfParticles()){
        while(j<m_system->getNumberOfParticles()){
            for(int k=0; k<i;k++){
                //first -= 2*(m_system->getNumberOfDimensions())*m_parameters[0];// 2*m_parameters[2];
                //if(a==0) second =0, fourth = 0, third = 0;
                //else {

                temp4 = a / (m_system->getDistanceMatrix()[k][j] * (m_system->getDistanceMatrix()[k][j] - a));
                fourth -= temp4*temp4;

                for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
                    temp6 -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);
                    //temp -= m_parameters[d] * (particles.at(j).getPosition()[d] - particles.at(k).getPosition()[d]);

                }

                second += temp6 * temp4;
                //temp *= 4*temp4;

                third += temp4 * a / (m_system->getDistanceMatrix()[k][i] * (m_system->getDistanceMatrix()[k][i]  - a));
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
    double constant;
    double R_kj;

    std::vector<double> QuantumForce(3);
    std::vector<double> u_deriv(3);
    for (int k = 0; k < m_system->getNumberOfParticles(); k++){

        for (int j = 0; j < k; j++){
                R_kj = m_system->getDistanceMatrixij(k,j);
                constant = 2*a / (R_kj*R_kj*(R_kj-a));
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


std::vector<double> SimpleGaussian::QuantumForceSingleParticle(std::vector<class Particle>& particles, int singParticle) {

    double a = m_system->getTrapSize();
    double R_kj;
    double constant;

    std::vector<double> QuantumForce(3);
    std::vector<double> u_deriv(3);
    for (int j = 0; j < singParticle; j++){
        R_kj = m_system->getDistanceMatrixij(singParticle,j);
        constant = a / (R_kj*R_kj*(R_kj-a));
        for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
            u_deriv[d] += (particles.at(singParticle).getPosition()[d] - particles.at(j).getPosition()[d]) * constant;
        }
    }

    for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
        QuantumForce[d] += -2 * (m_parameters[d]*particles.at(singParticle).getPosition()[d])
                               + u_deriv[d];
    }




    return QuantumForce;

}

double SimpleGaussian::computeDoubleDerivativeSingleParticle(std::vector<class Particle>& particles, int singParticle) {


    double first = 0;
    double a=0;
    //    double temp=0;
    //    double temp2=0;
    //    double temp3=0;
    double temp4= 0;
    //double temp5 = 0;
    double temp6 = 0;
    vector <double> m_parameters_squared(3);


    transform(m_parameters.begin(),m_parameters.end(),m_parameters.begin(),  std::bind1st (std::multiplies <double> () , 2.0));
    transform(m_parameters.begin(), m_parameters.end(), m_parameters.begin(),m_parameters_squared.begin(), multiplies<double>());
    transform(m_parameters.begin(),m_parameters.end(),m_parameters.begin(),  std::bind1st (std::multiplies <double> () , 1/2.0));

    for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
        first += (m_parameters_squared[d]*particles.at(singParticle).getPosition()[d]*particles.at(singParticle).getPosition()[d])
                -m_parameters_squared[d];
//        cout << "Particle " << singParticle << ": " << (m_parameters_squared[d]*particles.at(singParticle).getPosition()[d]*particles.at(singParticle).getPosition()[d])-m_parameters_squared[d] << endl;
        //*2
            //-m_parameters[0]*(2*beta_vector[d]);
            //temp5 -= m_parameters[d]*particles.at(k).getPosition()[d];
    }
//    cout << "First = " << first << endl;

    return first;
}
