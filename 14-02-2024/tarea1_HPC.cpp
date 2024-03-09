// Tarea 1  Santiago Pinilla Correa

#include<iostream>
#include<cmath>

int main(void){
    int n_min = 10;
    int n_max = 30;

    for(int ii = n_min; ii <= n_max; ii++){
        int count = 0; 
        for(int jj = 2; jj <= sqrt(ii); jj++){ 
            if(ii % jj == 0){
                count++;
                break; 
            }
        }
        if(count == 0){ 
            std::cout<< ii << "\n";
        }
    }

    return 0;
}