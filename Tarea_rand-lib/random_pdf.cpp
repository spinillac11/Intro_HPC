#include <random>
#include <iostream>
#include <cstdlib>
#include <vector>

void compute_pdf(int seed, int nsamples, double mu, double sigma, double xmin, double xmax, int nbins);

int main(int argc, char **argv)
{
    const int SEED = std::atoi(argv[1]);
    const int NSAMPLES = std::atoi(argv[2]);
    const double MU = std::atof(argv[3]);
    const double SIGMA = std::atof(argv[4]);
    const double XMIN = std::atof(argv[5]);
    const double XMAX = std::atof(argv[6]);
    const int NBINS = std::atoi(argv[7]);

    compute_pdf(SEED, NSAMPLES, MU, SIGMA, XMIN, XMAX, NBINS);
}

void compute_pdf(int seed, int nsamples, double mu, double sigma, double xmin, double xmax, int nbins)
{
    // random stuff
    std::mt19937 gen(seed);
    std::normal_distribution<double> dis{mu, sigma};
    // TODO: histogram stuff

    //Longitud de los intevalos
    double l = (xmax-xmin)/nbins;

    //numero de bin
    int nbin = 0;

    // Vector de conteo
    std::vector<double> counts(nbins, 0);

    for(int n = 0; n < nsamples; ++n) {
        double r = dis(gen);
        // TODO: fill here the counting histogram stuff

        // calculo del bin
        nbin = int((r-xmin)/l);
        // comprobar que esta en el intervalo
        if(nbin >= 0 && nbin < nbins){
            counts[nbin] += 1;
        }
    }
    // TODO: compute and print the pdf
    
    for(int ii = 0; ii < nbins; ii++){
        // valor del bin ii
        double xval = xmin + ii*l;
        //valor de la funcion densidad de probabilidad
        double pdf = counts[ii]/(l*nsamples);
        // imprimir datos
        std::cout << xval << '\t' << pdf << '\n';
    }
}
    

