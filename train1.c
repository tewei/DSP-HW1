#include <stdio.h>
#include <math.h>
#include "hmm.h"

#ifndef MAX_TIME
#   define MAX_TIME 200
#endif

int load_model(){
    //load model
    return 0;
}

//Implement Baum–Welch algorithm

double 

int train_model(HMM *hmm, int itr){
    //Initialize
    double epsilon[MAX_TIME][MAX_STATE][MAX_STATE];
    double gamma[MAX_TIME][MAX_STATE];
    double alpha[MAX_TIME][MAX_STATE];
    double beta[MAX_TIME][MAX_STATE];
    double *pi = hmm->initial;
    double *a = hmm->transition;
    double *b = hmm->observation;
    

    //Init
    int T = ;//time
    int N = ;//state

    //Calculate alpha
    //alpha[0][i]
    int i = 0;
    for(; i < N; ++i){
        alpha[0][i] = pi[i]*b[ ??? ][i];
    }
    //alpha[t][j]
    int t = 0, j = 0, i = 0;
    for(; t < T-1; ++t){
        int ot = ?; //what?
        for(j = 0; j < N){
            double sigma_alpha = 0;

            for(i = 0; i < N; ++i){
                sigma_alpha += alpha[t][i]*a[i][j];
            }
            alpha[t+1][j] = sigma_alpha*b[ot][j];

        }
    }
    //Calculate beta
    //beta[T][i]
    int i = 0;
    for(; i < N; ++i){
        beta[T-1][i] = 1;
    }
    //beta[t][j]
    int t = T-2, i = 0; j = 0;
    for(; t >= 0, --t){
        int ot = ?; //what?
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < N; ++j){
                beta[t][i] += a[i][j]*b[ot][j]*beta[t+1][j];
            }
        }
    }
    //Calculate epsilon
    int t = 0;

    for(; t < T-1; ++t){
        int ot = ?; //what?
        for(int i = 0; i < N; ++i){
            for(int j = 0; j < N; ++j){
                double p = 0,  q = 0;
                p = alpha[t][i]*a[i][j]*b[ot][j]*beta[t+1][j];
                q = ???; //P(O|Lambda)
                epsilon[t][i][j] = p/q;
            }
        }

    }
    //Calculate gamma


    //Reestimate pi
    //Reestimate a
    //Reestimate b



}


int main(){
    //Initialize
    HMM this_model;
    //Load
    char *file_name;
    loadHMM(, );
    //Train

    //Results
    char *model_file = 
    dumpHMM(, );
    
    return 0;
}


