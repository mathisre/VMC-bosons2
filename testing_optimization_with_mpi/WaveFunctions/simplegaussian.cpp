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
    //m_parameters.push_back(beta);
}

double SimpleGaussian::evaluate(std::vector<class Particle> &particles) {

    double r_squared = 0;
    double f=1;

    for(int j=0; j<m_system->getNumberOfParticles();j++){

        if(m_system->getNumberOfParticles()==1) break;
        for(int i=0; i<j; i++){
            f *= 1-m_system->getinteractionSize()/(m_system->getDistanceMatrixij(i,j));
    }
}

for(int i=0;i < m_system->getNumberOfParticles();i++){

        for(int d=0; d < m_system->getNumberOfDimensions();d++){
            r_squared += particles.at(i).getPosition()[d]*particles.at(i).getPosition()[d]*m_parameters[d];
                  }
}

return exp(-r_squared)*f;

}


double SimpleGaussian::computeDoubleDerivative(std::vector<class Particle>& particles) {

    double one=0;
    double interaction=0;

    double a=m_system->getinteractionSize();

    for(int i=0; i<m_system->getNumberOfParticles(); i++){
        double r_i_square=0;
        for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
         r_i_square += particles.at(i).getPosition()[d]*
                       particles.at(i).getPosition()[d];
         }
         int d = m_system->getNumberOfDimensions()-1;
         r_i_square += particles.at(i).getPosition()[d]*
                       particles.at(i).getPosition()[d]*m_parameters[2]/(m_parameters[0]);

       double second=0;
       double third=0;
       double fourth=0;
       double fifth=0;
       double temp;

       for(int j=0; j < i; j++) {

           double r_ij = m_system->getDistanceMatrixij(i,j);

           temp= a / ( (r_ij-a) * r_ij );

           second += temp;

           double r_ir_j = 0;
           for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
               r_ir_j += particles.at(i).getPosition()[d]*
                         particles.at(j).getPosition()[d];
           }
           int d = m_system->getNumberOfDimensions() - 1;
               r_ir_j += particles.at(i).getPosition()[d]*
                         particles.at(j).getPosition()[d]*
                         m_parameters[2]/(m_parameters[0]);

           fourth-= temp * temp;

           fifth -= 4 * m_parameters[0]  * (r_i_square - r_ir_j) * temp/
                   ( r_ij );

       }
       for(int j = j+1; j < m_system->getNumberOfParticles(); j++){

           double r_ij = m_system->getDistanceMatrixij(i,j);

           temp = a / ( (r_ij-a) * r_ij );

           second += temp;

           double r_ir_j = 0;
           for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
               r_ir_j += particles.at(i).getPosition()[d]*
                         particles.at(j).getPosition()[d];
           }
           int d = m_system->getNumberOfDimensions() - 1;
               r_ir_j += particles.at(i).getPosition()[d]*
                         particles.at(j).getPosition()[d]*
                         m_parameters[2]/(m_parameters[0]);

           fourth-= temp * temp;

           fifth -= 4 * m_parameters[0]  * (r_i_square - r_ir_j) * temp/
                   ( r_ij );
       }

       third=second*second;

       interaction+=second+third+fourth+fifth;

    }



    for(int i = 0; i < m_system->getNumberOfParticles(); i++){
        for(int d = 0; d < m_system->getNumberOfDimensions(); d++){
            one += m_parameters[d]*m_parameters[d]*
                   particles.at(i).getPosition()[d]*
                   particles.at(i).getPosition()[d];
        }
    }

    one*=4.0;
    one-= 2 * ( (m_system->getNumberOfDimensions() - 1) * m_parameters[0] + m_parameters[2])
            * m_system->getNumberOfParticles(); //constant term

    return one+interaction;
}

double SimpleGaussian::computeDoubleDerivativeSingleParticle(std::vector<class Particle>& particles, int singParticle) {

    double one=0;
    double interaction=0;
    int i = singParticle;

    double a=m_system->getinteractionSize();
    double r_i_square=0;
    for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
     r_i_square += particles.at(i).getPosition()[d]*
                   particles.at(i).getPosition()[d];
     }
     int d = m_system->getNumberOfDimensions()-1;
     r_i_square += particles.at(i).getPosition()[d]*
                   particles.at(i).getPosition()[d]*m_parameters[2]/(m_parameters[0]);

   double second=0;
   double third=0;
   double fourth=0;
   double fifth=0;
   double temp;

   for(int j=0; j < i; j++) {

       double r_ij = m_system->getDistanceMatrixij(i,j);

       temp= a / ( (r_ij-a) * r_ij );

       second += temp;

       double r_ir_j = 0;
       for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
           r_ir_j += particles.at(i).getPosition()[d]*
                     particles.at(j).getPosition()[d];
       }
       int d = m_system->getNumberOfDimensions() - 1;
           r_ir_j += particles.at(i).getPosition()[d]*
                     particles.at(j).getPosition()[d]*
                     m_parameters[2]/(m_parameters[0]);

       fourth-= temp * temp;

       fifth -= 4 * m_parameters[0]  * (r_i_square - r_ir_j) * temp/
               ( r_ij );

   }
   for(int j = j+1; j < m_system->getNumberOfParticles(); j++){

       double r_ij = m_system->getDistanceMatrixij(i,j);

       temp = a / ( (r_ij-a) * r_ij );

       second += temp;

       double r_ir_j = 0;
       for(int d = 0; d < m_system->getNumberOfDimensions() - 1; d++){
           r_ir_j += particles.at(i).getPosition()[d]*
                     particles.at(j).getPosition()[d];
       }
       int d = m_system->getNumberOfDimensions() - 1;
           r_ir_j += particles.at(i).getPosition()[d]*
                     particles.at(j).getPosition()[d]*
                     m_parameters[2]/(m_parameters[0]);

       fourth-= temp * temp;

       fifth -= 4 * m_parameters[0]  * (r_i_square - r_ir_j) * temp/
               ( r_ij );
   }

   third=second*second;

   interaction+=second+third+fourth+fifth;


   for(int i = 0; i < m_system->getNumberOfParticles(); i++){
       for(int d = 0; d < m_system->getNumberOfDimensions(); d++){
           one += m_parameters[d]*m_parameters[d]*
                  particles.at(i).getPosition()[d]*
                  particles.at(i).getPosition()[d];
       }
   }

   one*=4.0;
   one-= 2 * ( (m_system->getNumberOfDimensions() - 1) * m_parameters[0] + m_parameters[2])
           * m_system->getNumberOfParticles(); //constant term

   return one+interaction;
}



std::vector<vector<double>> SimpleGaussian::QuantumForce(std::vector<class Particle>& particles) {
//    // CREATE THE MATRIX
//    m_quantumForce.resize(m_system->getNumberOfDimensions());
//    for(int i = 0; i < m_system->getNumberOfDimensions(); i++)
//        m_quantumForce[i].resize(2*m_system->getNumberOfParticles());

//    // ASSIGN THE VALUES
//    double a = m_parameters[2];
//    for(int i = 0; i < m_system->getNumberOfDimensions(); i++)
//        for(int j = 0; j < m_system->getNumberOfParticles(); j++){
//            m_quantumForce[i][j] = -2*particles[j]->getPosition()[i]*m_parameters[0];
//            for(int k = 0; k < j; k++){
//                double rjk= m_system->getDistance()[j][k];
//                m_quantumForce[i][j]+= a* (particles[j]->getPosition()[i]-particles[k]->getPosition()[i])/((rjk-a)*rjk*rjk);
//            }
//            for(int k = j+1; k < m_system->getNumberOfParticles(); k++){
//                double rjk= m_system->getDistance()[k][j];
//                m_quantumForce[i][j]+= a*(particles[j]->getPosition()[i]-particles[k]->getPosition()[i])/((rjk-a)*rjk*rjk);
//            }
//            m_quantumForce[i][j]*=2;
//        }



    double a = m_system->getinteractionSize() ;
    double constant;
    double R_kj;
    double dimension=m_system->getNumberOfDimensions();
    double number =m_system->getNumberOfParticles();
    std::vector<std::vector<double>> QuantumForce(dimension,vector<double>(number));
    std::vector<double> u_deriv(3);
    for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
        for (int k = 0; k < m_system->getNumberOfParticles(); k++){
    QuantumForce[d][k] = -2 * (m_parameters[d]*particles.at(k).getPosition()[d]);
        for (int j = 0; j < k; j++){
                R_kj = m_system->getDistanceMatrixij(k,j);
                constant = 2*a / (R_kj*R_kj*(R_kj-a));

                    QuantumForce[d][k] += (particles.at(k).getPosition()[d] - particles.at(j).getPosition()[d]) * constant;

            }
        }
    }
    return QuantumForce;
}


std::vector<vector<double>> SimpleGaussian::QuantumForceSingleParticle(std::vector<class Particle>& particles, int singParticle) {

    double a = m_system->getinteractionSize();
    double R_kj;
    double constant;

    double dimension=m_system->getNumberOfDimensions();
    double number =m_system->getNumberOfParticles();
    std::vector<std::vector<double>> QuantumForce(dimension,vector<double>(number));
    //std::vector<double> u_deriv(3);
    int k=singParticle;
    for (int d = 0; d < m_system->getNumberOfDimensions(); d++){

    QuantumForce[d][k] = -2 * (m_parameters[d]*particles.at(k).getPosition()[d]);
        for (int j = 0; j < k; j++){
                R_kj = m_system->getDistanceMatrixij(k,j);
                constant = 2*a / (R_kj*R_kj*(R_kj-a));

                    QuantumForce[d][k] += (particles.at(k).getPosition()[d] - particles.at(j).getPosition()[d]) * constant;

            }

    }
    return QuantumForce;

}


