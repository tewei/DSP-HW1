#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "hmm.h"

#ifndef MAX_SAMPLES
#   define MAX_SAMPLES 10000
#endif

//Global Variables
double alpha[MAX_SAMPLES][MAX_SEQ][MAX_STATE];
double beta[MAX_SAMPLES][MAX_SEQ][MAX_STATE];
double epsilon[MAX_SAMPLES][MAX_SEQ][MAX_STATE][MAX_STATE];
double ggmma[MAX_SAMPLES][MAX_SEQ][MAX_STATE];

int sample_seq[MAX_SAMPLES][MAX_SEQ];
int sample_seq_len[MAX_SAMPLES];

int iter = 0; // Number of Training Iterations
int N; //Number of States

void printHMM(HMM *hmm){
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            printf("a%d%d %.4f ", i, j, hmm->a[i][j]);
        }
        printf("\n");
    }
}
void print_sample(int s){
    for(int i = 0; i < sample_seq_len[s]; ++i){
        printf("%d ", sample_seq[s][i]);
    }
    printf("\n");
}

void trainHMM(HMM *hmm, int iter, FILE *file){

    int num_samples = 0;
    int num_symbols = hmm->observ_num; // Number of Observed Symbols
    N = hmm->state_num; // Number of States


    //Process Raw Data
    int s = 0;
    char raw_sample[MAX_SEQ];
    while(fgets(raw_sample, MAX_SEQ, file) != NULL){
        int i = 0;
        sample_seq_len[s] = strlen(raw_sample);
        for(; i < sample_seq_len[s]; ++i){
            if(raw_sample[i] == '\n' || raw_sample[i] == '\0') continue;
            sample_seq[s][i] = raw_sample[i] - 'A';
        }
        s++;
    }
    num_samples = s;

    //Start Training
    int epoch = 0;
    for(epoch = 0; epoch < iter; ++epoch){

        int s = 0, i = 0, j = 0, k = 0;
        
        for(s = 0; s < num_samples; ++s){
            //print_sample(s);
            //Init
            int T = sample_seq_len[s], t = 0;

            //Calculate alpha
            //alpha[s][0][i]
            for(i = 0; i < N; ++i){
                int ot = sample_seq[s][0];
                alpha[s][0][i] = hmm->pi[i] * hmm->b[ot][i];
            }
            //alpha[s][t][j]

            for(t = 0; t < T-1; ++t){
                int ot = sample_seq[s][t+1];
                for(j = 0; j < N; ++j){
                    double sigma_alpha = 0;

                    for(i = 0; i < N; ++i){
                        sigma_alpha += alpha[s][t][i] * hmm->a[i][j];
                        //printf("%.10f \n", sigma_alpha);
                    }
                    alpha[s][t+1][j] = sigma_alpha * hmm->b[ot][j];
                }
            }

            //Calculate beta
            //beta[s][T][i]
            for(i = 0; i < N; ++i){
                beta[s][T-1][i] = 1;
            }
            //beta[s][t][j]
            for(t = T-2; t >= 0; --t){
                int ot = sample_seq[s][t];
                for(i = 0; i < N; ++i){
                    double sum_aij_bj_ot = 0;
                    for(j = 0; j < N; ++j){
                        sum_aij_bj_ot += hmm->a[i][j] * hmm->b[ot][j] * beta[s][t+1][j];
                    }
                    beta[s][t][i] = sum_aij_bj_ot;
                }
            }

            //Calculate epsilon

            for(t = 0; t < T; ++t){
                int ot = sample_seq[s][t];
                double prob_observ_seq = 0;
                for(i = 0; i < N; ++i){
                    for(j = 0; j < N; ++j){
                        double prob_ij_t = alpha[s][t][i] * hmm->a[i][j] * hmm->b[ot][j] * beta[s][t+1][j];
                        
                        prob_observ_seq += prob_ij_t;
                        epsilon[s][t][i][j] = prob_ij_t;
                    }
                }

                for(i = 0; i < N; ++i){
                    for(j = 0; j < N; ++j){
                        epsilon[s][t][i][j] /= prob_observ_seq;
                    }
                }
            }
            //Calculate gamma
            for(t = 0; t < T; ++t){
                double prob_observ_seq = 0;
                for(j = 0; j < N; ++j){
                    prob_observ_seq += alpha[s][t][j] * beta[s][t][j];
                }
                for(i = 0; i < N; ++i){
                    ggmma[s][t][i] = alpha[s][t][i] * beta[s][t][i] / prob_observ_seq;
                }
            }

        }
        
        //Update Pi
        
        for(i = 0; i < N; ++i){
            double avg_pi = 0;
            for(s = 0; s < num_samples; ++s){
                avg_pi += ggmma[s][0][i];
            }
            avg_pi /= num_samples;
            hmm->pi[i] = avg_pi;
        }

        //Update Aij
        for(i = 0; i < N; ++i){
            for(j = 0; j < N; ++j){

                double avg_aij = 0, sum_P = 0, sum_Q = 0;
                for(s = 0; s < num_samples; ++s){
                    int T = sample_seq_len[s], t = 0;
                    for(t = 0; t < T-1; ++t){
                        sum_P += epsilon[s][t][i][j];
                        sum_Q += ggmma[s][t][i];
                    }
                }
                avg_aij = sum_P/sum_Q;
                //printf("Debug %f %f\n", sum_P, sum_Q);
                hmm->a[i][j] = avg_aij;

            }
        }

        //Update Bj(k)
        for(j = 0; j < N; ++j){
            for(k = 0; k < num_symbols; ++k){

                double avg_bjk = 0, sum_P = 0, sum_Q = 0;
                int s = 0;
                for(s = 0; s < num_samples; ++s){
                    int T = sample_seq_len[s], t = 0;
                    for(t = 0; t < T; ++t){
                        int ot = sample_seq[s][t];
                        if(ot == k){
                            sum_P += ggmma[s][t][j];
                        }
                        sum_Q += ggmma[s][t][j];
                    }
                }

                avg_bjk = sum_P/sum_Q;
                
                hmm->b[k][j] = avg_bjk;
            }
        }
        printHMM(hmm);
        printf("Training iteration %d completed!\n", epoch+1);

    }
}


int main(int argc, char *argv[]){
    
    //Handle Command Line Arguments
    if(argc < 5){
        perror("Too few arguments, expected 4");
    }else if(argc > 5){
        perror("Too many arguments, expected 4");
    }
    
    iter = atoi(argv[1]);
    
    //Open Files
    FILE *training_data_file = open_or_die(argv[3], "r");
    FILE *trained_model_file = open_or_die(argv[4], "w");
    //Load Initial Model
    HMM hmm;
    loadHMM(&hmm, argv[2]);

    //Training
    trainHMM(&hmm, iter, training_data_file);

    //Results
    dumpHMM(trained_model_file, &hmm);

    return 0;
}