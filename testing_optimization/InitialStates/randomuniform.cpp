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

RandomUniform::RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double trapSize, double timeStep)  :
        InitialState(system) {
    assert(numberOfDimensions > 0 && numberOfParticles > 0);
    m_numberOfDimensions = numberOfDimensions;
    m_numberOfParticles  = numberOfParticles;

    /* The Initial State class is in charge of everything to do with the
     * initialization of the system; this includes determining the number of
     * particles and the number of dimensions used. To make sure everything
     * works as intended, this information is passed to the system here.
     */
    m_system->setNumberOfDimensions(numberOfDimensions);
    m_system->setNumberOfParticles(numberOfParticles);
    m_system->setTrapSize(trapSize);
    m_system->setTimeStep(timeStep);
    m_system->setSqrtTimeStep(sqrt(timeStep));

    setupInitialState();
}

void RandomUniform::setupInitialState() {
    for (int i=0; i < m_numberOfParticles; i++) {
        std::vector<double> position = std::vector<double>();

        for (int j=0; j < m_numberOfDimensions; j++) {
            /* This is where you should actually place the particles in
             * some positions, according to some rule. Since this class is
             * called random uniform, they should be placed randomly according
             * to a uniform distribution here. However, later you will write
             * more sub-classes of the InitialState class in which the
             * particles are placed in other configurations.
             *
             * Note: For now, the particles are simply placed in positions
             * according to their index in the particles list (this is
             * obviously NOT a good idea).
             */

            position.push_back(Random::nextDouble()-0.5);
            //position.push_back(i);
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
            if (R_ki < m_system->getTrapSize()){
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
