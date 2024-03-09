#include <random>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>

typedef std::vector<double> vec; // Nombre para vector de doubles

using fptr=double(double, vec, int);



void phi_i(vec & phi, const double dx, double x, const int nbins);
void Dphi_i(vec & dphi, const double dx, double x, const int nbins);
void global_m(vec & eps, vec & phi, const int nbins, const double dx);
double M(double x, vec phi, int n);
void simpson(fptr func, vec & phi, vec & I_phi, const int nbins, const double dx, const double a, const double b, const int npoint);

double F = 1;
double W = 1;
double L = 1;

int main(int argc, char **argv){
    //# puntos para elementos finitos
    int N = 4;
    // Tamano de los intervalos
    double DX = (L/2)*(1.0/(N-1));
    //# puntos para la integral
    int n_points = 100;
    //Modulo de Young y momento de inercia a la flexion
    double E = 1;
    double I = 1;
    
    //Condiciones de frontera    
    double Etha_o = 1;
    double Etha_l = 0;

    

    // Funciones phi
    vec PHI(N, 0.0);
    // Matriz global
    vec EPS(N*N, 0.0);
    // Vector de integrales
    vec I_PHI(N, 0.0);
    // Vector de condiciones de frontera
    vec RHS(N, 0.0);
    RHS[0] = Etha_o;
    RHS[N-1] = Etha_l;

    //Calcular matriz global
    global_m(EPS, PHI, N, DX);

    //Imprimir matriz
    std::ofstream f1out("matrix.txt");
    for (int ii = 0; ii < N; ++ii) {
        for (int jj = 0; jj < N; ++jj) {
            f1out << EPS[ii*N + jj] <<  "\t";
        }
        f1out << "\n";
    }
    f1out.close();

    //Calcular integrales de M*phi_i
    simpson(M, PHI, I_PHI, N, DX, 0, L/2, n_points);

    //Imprimir vector
    std::ofstream f2out("vec.txt");
    for (int ii = 0; ii < N; ++ii) {
        f2out << -I_PHI[ii] + RHS[ii] << "\n";
    }
    f2out.close();




    return 0;
} 

// Funciones phi_i por intervalos
void phi_i(vec & phi, const double dx, double x, const int nbins){
    for(int ii = 0; ii < nbins; ii++){

        //Recta creciente
        if(dx*(ii-1) <= x && x <= dx*ii){
            phi[ii] = (1.0/dx)*(x-dx*(ii-1));
        }
        //Recta decreciente 
        else if(dx*ii < x && x <= dx*(ii+1)){
            phi[ii] = -(1.0/dx)*(x-dx*ii)+1.0;
        }
        else{
            phi[ii] = 0.0;
        }
    }
}

// Derivadas de dphi_i por intervalos 
void Dphi_i(vec & dphi, const double dx, double x, const int nbins){
    for(int ii = 0; ii < nbins; ii++){

        //Recta creciente
        if(dx*(ii-1) <= x && x <= dx*ii){
            dphi[ii] = 1.0/dx;
        }
        //Recta decreciente 
        else if(dx*ii < x && x <= dx*(ii+1)){
            dphi[ii] = -1.0/dx;
        }
        else{
            dphi[ii] = 0.0;
        }
    }
}
//Calcular matriz global a partir de las derivadas de phi_i
void global_m(vec & eps, vec & dphi, const int nbins, const double dx){

    for(int n = 0; n < nbins-1; n++){
        //calcular x en la mitad del intervalo
        double xi = (2*n+1)*dx/2;
        Dphi_i(dphi, dx, xi, nbins);
        
        for(int i = 0; i < nbins; i++){
            for(int j = 0; j < nbins; j++){
                eps[i*nbins+j] += dx*(dphi[i])*(dphi[j]);
            }
        }
    }
}

//Momento flector
double M(double x, vec phi, int n){
    return (-W/(2*L)*x*x + (F+W)*x/2)*phi[n];
}
//Integracion por metodo simpson
void simpson(fptr func, vec & phi, vec & I_phi, const int nbins, const double dx, const double a, const double b, const int npoint){
    
    for(int jj = 0; jj < nbins; jj++){
        // check for even number of points
        int nlocal = npoint;
        if (nlocal%2 != 0) {
            nlocal += 1;
        }

        phi_i(phi, dx, a, nbins);
        double resulta = func(a, phi, jj);
        
        phi_i(phi, dx, b, nbins);
        double resultb = func(b, phi, jj);

        double result = resulta + resultb;

        double sum = 0;
        double x;
        const double h = (b-a)/nlocal;

        // first sum
        sum = 0;
        for(int ii = 1; ii <= nlocal/2 - 1; ++ii ) {
            x = a + 2*ii*h;
            phi_i(phi, dx, x, nbins);
            sum += func(x, phi, jj);
        }
        result += 2*sum;

        // second sum
        sum = 0;
        for(int ii = 1; ii <= nlocal/2; ++ii ) {
            x = a + (2*ii-1)*h;
            phi_i(phi, dx, x, nbins);
            sum += func(x, phi, jj);
        }
        result += 4*sum;

        I_phi[jj] = result*h/3;
    }
}