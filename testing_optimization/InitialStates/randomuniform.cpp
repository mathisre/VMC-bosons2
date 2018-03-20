#include "randomuniform.h"
#include <iostream>
#include <cassert>
#include "../Math/random.h"
#include "../particle.h"
#include "../system.h"
#include <cmath>

using std::cout;
using std::endl;

using namespace std;

RandomUniform::RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double interactionSize, double timeStep)  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    m_system->setinteractionSize(interactionSize);
    m_system->setTimeStep(timeStep);
    m_system->setSqrtTimeStep(sqrt(timeStep));

    setupInitialState();
}

void RandomUniform::setupInitialState() {
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (int j=0; j < m_numberOfDimensions; j++) {

            position.push_back(Random::nextDouble()-0.5);
        }
        Particle p;
        m_particles.push_back(p);
        m_particles.at(i).setNumberOfDimensions(m_numberOfDimensions);
        m_particles.at(i).setPosition(position);
    }
}

void RandomUniform::setupInitialStateWithInteraction() {
    int placedParticles = 0;
    double R_ki = 0;
    vector<vector<double>> distancematrix(m_numberOfParticles, vector<double>(m_numberOfParticles));

    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();
        for (int d=0; d < m_numberOfDimensions; d++) {
            position.push_back(Random::nextDouble()-0.5);
        }
        for (int k = 0; k < placedParticles + 1; k++){
            for (int d=0; d < m_numberOfDimensions; d++) {
                R_ki += (m_particles.at(k).getPosition()[d] - m_particles.at(i).getPosition()[d]) *
                        (m_particles.at(k).getPosition()[d] - m_particles.at(i).getPosition()[d]);
            }
            R_ki = sqrt(R_ki);
            if (R_ki < m_system->getinteractionSize()){
                Particle p;
                m_particles.push_back(p);
                m_particles.at(i).setNumberOfDimensions(m_numberOfDimensions);
                m_particles.at(i).setPosition(position);

                distancematrix[k][i] = R_ki;
                distancematrix[i][k] = R_ki;
                placedParticles++;
            }
        }
    }
    m_system->setDistanceMatrix(distancematrix);
}
