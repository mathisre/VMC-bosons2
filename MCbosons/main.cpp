#include <iostream>
#include <armadillo>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <time.h>
#include <random>
#include <omp.h>
#include <vector>
#include "solver.h"


using namespace arma;
using namespace std;
using std::vector;
ofstream ofile;


double E_local(vec R){

}


double psi_local(vec R){

}

void update_avg(){

}

void metropolis(int N_accepted){
    N_accepted += 1;

}

int main(int argc, char* argv[])
{
    string filename;
    filename = "ising_result.txt";
    ofile.open(filename);


    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> RNG(0.0,1.0);

    int nDimensions = 3;

    int nParticles = 10;
    int MC_cycles = 100;
    vector<double> r_Old (nParticles, nDimensions);
    vec r_Old = randu<vec>(N,nDimension);
    vec r_New = zeros<vec>(N,nDimension);

    double rng;
    int N_accepted = 0;
    vec alpha = linspace<vec>(0.5,1.5,100);

    double local_E;
    double local_E_new;
    double step;

    double WaveFunctionOld;
    double WaveFunctionNew;

    //maybe an alpha loop?

    for (int i = 0; i<MC_cycles; i++){ //MC loop

        for (int k = 0; k<nParticles; k++){ // particle loop(??)

            rng = RNG(gen);
            r_New = r_Old + rng*step;
            WaveFunctionNew = solver.WaveFunction(r_New, alpha);

            solver.Metropolis(rng);
            }
        }



    }
    return 0;
}


