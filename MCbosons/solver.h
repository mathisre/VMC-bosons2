#ifndef SOLVER_H
#define SOLVER_H


class solver
{
public:
    solver();



    void Metropolis(double rng, double w, double &Psi_Old);
    double Psi_local();
    double E_local();
    double WaveFunction(vector<double> R, double alpha);




    int getDimensions() const;
    void setDimensions(int value);

    int getNumberOfSteps() const;
    void setNumberOfSteps(int value);

    double getStepSize() const;
    void setStepSize(double value);

    int getNumerOfParticles() const;
    void setNumerOfParticles(int NumerOfParticles);

private:
    int m_NumerOfParticles;
    int m_Dimensions;
    int m_NumberOfSteps;
    double m_StepSize;
    double m_a;
};


#endif // SOLVER_H
