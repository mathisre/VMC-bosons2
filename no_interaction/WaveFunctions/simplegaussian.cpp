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

    double r_squared = 0;
    double f=1;

    for(int j=0; j<m_system->getNumberOfParticles();j++){
        if(m_system->getNumberOfParticles()==1) break;
        for(int i=0; i<j; i++){
            f = 1-m_system->getInteractionSize()/(m_system->getDistanceMatrixij(i,j));
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


    double first,second,third, fourth=0;
    first = 0;
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
        }
    }
    double r_squared = 0;
    for (int i=0; i < m_system->getNumberOfParticles(); i++) {
        for (int d=0; d < m_system->getNumberOfDimensions() ; d++) {
          r_squared += m_system->getParticles().at(i).getPosition()[d]*m_system->getParticles().at(i).getPosition()[d];
        }
    }
    first = 4*r_squared*m_parameters[0]*m_parameters[0] - 2*m_system->getNumberOfDimensions()*m_parameters[0]*m_system->getNumberOfParticles();

    return first;
}



std::vector<vector<double>> SimpleGaussian::QuantumForce(std::vector<class Particle>& particles) {

    double a = m_system->getInteractionSize() ;
    double constant;
    double R_kj;
    double dimension=m_system->getNumberOfDimensions();
    double number =m_system->getNumberOfParticles();
    std::vector<std::vector<double>> QuantumForce(dimension,vector<double>(number));

    for (int d = 0; d < m_system->getNumberOfDimensions(); d++){
        for (int k = 0; k < m_system->getNumberOfParticles(); k++){

            QuantumForce[d][k] = -4 * (m_parameters[d]*particles.at(k).getPosition()[d]);

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

    double a = m_system->getInteractionSize();
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
    }
    int i = singParticle;

    double r_squared = 0;
        for (int d=0; d < m_system->getNumberOfDimensions() ; d++) {
          r_squared += m_system->getParticles()[i].getPosition()[d]*m_system->getParticles()[i].getPosition()[d];
        }

    first = 4*r_squared*m_parameters[0]*m_parameters[0] -
            2*m_system->getNumberOfDimensions()*m_parameters[0]*m_system->getNumberOfParticles();

    return first;
}
