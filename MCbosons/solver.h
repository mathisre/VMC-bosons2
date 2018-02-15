#ifndef SOLVER_H
#define SOLVER_H


class solver
{
public:
    solver();



    void Metropolis(double rng, double w);
    double Psi_local();
    double E_local();



    int getNumerOfAtoms() const;
    void setNumerOfAtoms(int value);

    int getDimensions() const;
    void setDimensions(int value);

    int getNumberOfSteps() const;
    void setNumberOfSteps(int value);

    double getStepSize() const;
    void setStepSize(double value);

private:
    int m_NumerOfAtoms;
    int m_Dimensions;
    int m_NumberOfSteps;
    double m_StepSize;
};


#endif // SOLVER_H
