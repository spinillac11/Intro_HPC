#include<iostream>
#include<cstdlib>
#include<cmath>

float sum1(int N); //Declaracion suma 1
float sum2(int N); //Declaracion suma 2

int main(int argc, char **argv){
    float S1 = 0;
    float S2 = 0;
    float dif_rel = 0;
    float dif_por = 0;

    for(int ii = 1; ii <= 1000; ii++){
        int n = 1000*ii; //Pasos de a 1000
        S1 = sum1(n);
        S2 = sum2(n);
        dif_rel = std::fabs(1.0-(S1/S2));
        dif_por = dif_rel*100; //Diferencia porcentual

        //Imprimir N, diferencia relativa y porcentual
        std::cout << n << "\t" << dif_rel << "\t" << dif_por << "\n";
    }
    return 0;
}

float sum1(int N){
    float sum = 0;
    for(int k = 1; k <= N; k++){
        sum += 1.0/k;
    }
    return sum;
}

float sum2(int N){
    float sum = 0;
    for(int k = N; k >= 1; k--){
        sum += 1.0/k;
    }
    return sum;
}