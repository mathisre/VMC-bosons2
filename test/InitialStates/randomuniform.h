#pragma once
#include "initialstate.h"

class RandomUniform : public InitialState {
public:
    RandomUniform(System* system, int numberOfDimensions, int numberOfParticles, double trapSize, double timeStep);
    void setupInitialState();
    void setupInitialStateWithInteraction();
};

