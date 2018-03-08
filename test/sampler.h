#pragma once
#include <string>
class Sampler {
public:
    Sampler(class System* system);
    void setNumberOfMetropolisSteps(int steps);
    void sample(bool acceptedStep);
    void printOutputToTerminal();

    void openDataFile(std::string filename);
    void writeToFile();

    void computeAverages();
    double getEnergy()          { return m_energy; }
    int getStepNumber() const;

    int getNumberOfMetropolisSteps() const;

private:
    int     m_numberOfMetropolisSteps  = 0;
    int     m_stepNumber               = 0;
    double  m_energy                   = 0;
    double  m_cumulativeEnergy         = 0;
    double  m_cumulativeWF             = 0;
    std::string  m_filename;
    class System* m_system = nullptr;
};
