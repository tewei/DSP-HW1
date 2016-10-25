#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "hmm.h"

#ifndef MAX_SAMPLES
#   define MAX_SAMPLES 10000
#endif

#ifndef MAX_MODELS
#   define MAX_MODELS 10
#endif
//Global variables

HMM hmms[MAX_MODELS];

int sample_seq[MAX_SAMPLES][MAX_SEQ];
int sample_seq_len[MAX_SAMPLES];
double delta[MAX_SAMPLES][MAX_SEQ][MAX_STATE];
int psi[MAX_SAMPLES][MAX_SEQ][MAX_STATE];

void printHMM(HMM *hmm){
    int N = 6;
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            printf("a%d%d %.4f ", i, j, hmm->a[i][j]);
        }
        printf("\n");
    }
}

double accuracy(FILE *testing_result_file, FILE *ground_truth_file){
    char prediction[MAX_LINE], truth[MAX_LINE];
    int num_samples = 0, correct_samples = 0;
    
    while(fgets(prediction, MAX_LINE, testing_result_file) != NULL && fgets(truth, MAX_LINE, ground_truth_file) != NULL){    
        if(strcmp(prediction, truth) == 0) correct_samples++;
        num_samples++;
    }

    return correct_samples/(double) num_samples;

}

double viterbi(HMM *hmm, int s){
    
    int num_symbols = hmm->observ_num; // Number of Observed Symbols
    int N = hmm->state_num; // Number of States
    int T = sample_seq_len[s];
    //Initialize
    int ot = sample_seq[s][0], i = 0;
    for(i = 0; i < N; ++i){
        delta[s][0][i] = hmm->pi[i] * hmm->b[ot][i];
    }
    //Recursion
    int t = 0, j = 0;
    for(t = 0; t < T-1; ++t){
        int ot = sample_seq[s][t+1];
        for(j = 0; j < N; ++j){
            int arg_max = -1;
            double max = 0.0;
            for(i = 0; i < N; ++i){
                if(delta[s][t][i] * hmm->a[i][j] > max){
                    max = delta[s][t][i] * hmm->a[i][j];
                    arg_max = i;
                }
            }
            delta[s][t+1][j] = max * hmm->b[ot][j];
            psi[s][t+1][j] = arg_max;
        }
    }
    //Termination
    double P_star = 0;
    int qT_star = -1;
    //for(i = 0; i < N; ++i) printf("jizz %f\n", delta[s][3][i]);
    for(i = 0; i < N; ++i){
        if(delta[s][T-1][i] > P_star){
            P_star = delta[s][T-1][i];
            qT_star = i;
        }
    }
    return P_star;
}
void testHMM(FILE *modellist_file, FILE *testing_data_file, FILE *testing_result_file){
    //Load models
    int num_models = 0;
    char model_file_name[MAX_LINE];
    while(fgets(model_file_name, MAX_SEQ, modellist_file) != NULL){    
        model_file_name[strlen(model_file_name)-1] = '\0';
        loadHMM(&hmms[num_models++], model_file_name);
    }
    //Load samples
    int num_samples = 0;
    int s = 0;
    char raw_sample[MAX_SEQ];
    while(fgets(raw_sample, MAX_SEQ, testing_data_file) != NULL){
        int i = 0;
        sample_seq_len[s] = strlen(raw_sample);
        for(; i < sample_seq_len[s]; ++i){
            if(raw_sample[i] == '\n' || raw_sample[i] == '\0') continue;
            sample_seq[s][i] = raw_sample[i] - 'A';
        }
        s++;
    }
    num_samples = s;
    for(s = 0; s < num_samples; ++s){
        int m = 0, arg_max = -1;
        double max = 0.0;
        for(m = 0; m < num_models; ++m){
            double prob_m_O = viterbi(&hmms[m], s);
            if(prob_m_O > max){
                max = prob_m_O;
                arg_max = m;
            }
        }
        //printf("sample %d is class %s with probability %.100f\n", s, hmms[arg_max].model_name, max);
        fprintf(testing_result_file, "%s\n", hmms[arg_max].model_name);
    }
}

int main(int argc, char *argv[]){

    if(argc < 4){
        perror("Too few arguments, expected 4");
    }else if(argc > 4){
        perror("Too many arguments, expected 4");
    }

    //Open files
    FILE *modellist_file = open_or_die(argv[1], "r");
    FILE *testing_data_file = open_or_die(argv[2], "r");
    FILE *testing_result_file = open_or_die(argv[3], "a+");
    //Run prediction
    testHMM(modellist_file, testing_data_file, testing_result_file);
    fclose(modellist_file);
    fclose(testing_data_file);
    fclose(testing_result_file);
    
    //Accuracy
    FILE *ground_truth_file = open_or_die("testing_answer.txt", "r");
    testing_result_file = open_or_die(argv[3], "r");
    double testing_accuracy = accuracy(testing_result_file, ground_truth_file);
    printf("Testing complete, accuracy = %f", testing_accuracy);
    
    return 0;
}