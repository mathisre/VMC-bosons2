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
using namespace arma;
using namespace std;
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

    int N = 10;
    int MC_cycles = 100;
    vec r_Old = randu<vec>(N,nDimension);
    vec r_New = zeros<vec>(N,nDimension);

    double rng;
    double step;
    int N_accepted = 0;
    vec alpha = linspace<vec>(0.5,1.5,100);

    double local_E;
    double local_E_new;
    double W;

    //maybe an alpha loop?

    for (int i = 0; i<MC_cycles; i++){ //MC loop

        for (int k = 0; k<N; k++){ // particle loop(??)

            rng = RNG(gen);

            if (RNG(gen) <= W){
                metropolis(N_accepted);
            }
        }



    }
    return 0;
}


