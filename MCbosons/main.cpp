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

    // Initialize uniform RNG un [0,1]
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Set up the uniform distribution for x \in [[0, 1]
    std::uniform_real_distribution<double> RNG(0.0,1.0);

    //maybe write a program that can take dimensions as argument, instead of writing three programs for 1d, 2d and 3d

    int N = 10;
    int MC_cycles = 100;
    vec R = randu<vec>(N);

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
            R(k) += rng*step;
            //Energy as func of R
            local_E_new = E_local(R);
            W = local_E_new / local_E;
            if (RNG(gen) <= W){
                metropolis(N_accepted);
            }
        }



    }
    return 0;
}


