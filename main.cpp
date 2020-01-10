#include <iostream>
#include "constants.h"
#include "helper_functions.h"
#include "Pair_HMM.h"
#include <math.h>
int main() {

    int i = 0;

    // reading emission and transitions_al matrices from files.

    double *emissions_al = read_matrix_file(PATH_EMISSIONS_FILE_ALIIGN, N_STATES, N_OBSERAVBLES);
    double *transitions_al = read_matrix_file(PATH_TRANSITIONS_FILE_ALIGN, N_STATES, N_STATES);

    double *emissions_rand = read_matrix_file(PATH_EMISSIONS_FILE_RANDOM, N_STATES, N_OBSERAVBLES);
    double *transitions_rand = read_matrix_file(PATH_TRANSITIONS_FILE_RANDOM, N_STATES, N_STATES);

    std::ofstream alignment_file_by_task(ALIGNMENT_OUTPUT_FILE.c_str());
    std::ofstream alignment_file(ALIGNMENT_OUTPUT_FILE_FULL.c_str());
    std::ofstream paths_file(ALIGNMENT_OUTPUT_FILE_PATH.c_str());

    std::ofstream random_prob_file(PROBAB_RANDOM_FILE.c_str());
    std::ofstream alignment_prob_file(PROBAB_ALIGNMENT_FILE.c_str());
    std::ofstream log_ratios_file(LOG_RATIOS_FILE.c_str());


    //print_matrix(emissions_al, 5, 25);
    //print_matrix(transitions_al, 5, 5);

    std::cout<<"Alignment pair HMM \n";
    Pair_HMM align_hmm(N_STATES, N_OBSERAVBLES, transitions_al, emissions_al);
    std::cout<<"Random pair HMM \n";
    Pair_HMM random_hmm(N_STATES, N_OBSERAVBLES, transitions_rand, emissions_rand);

    align_hmm.set_model_name("Alignment model");
    random_hmm.set_model_name("Random model");

    std::string s1 = "";
    std::string s2 = " " ;
    double alignment_prob,random_prob,log_odds;

    // Reading sequences in_pairs and calcualting in Alignment model
    int n_sequence_to_read =20;
    std::string* sequences = sequences_reader(PATH_INPUT_FILE,n_sequence_to_read);

    for (i=0;i<n_sequence_to_read;i+=2){
        //std::cout<<i<<std::endl;
        s1 = sequences[i];
        s2 = sequences[i+1];

        std::cout<<"X: "<<s1<<" "<<s1.length()<<'\n'<<"Y: "<<s2<<" "<<s2.length()<<std::endl;

        align_hmm.set_observations_x(s1);
        align_hmm.set_observations_y(s2);

        random_hmm.set_observations_x(s1);
        random_hmm.set_observations_y(s2);

        align_hmm.calculate_viterbi_alignment ();

        alignment_prob =  align_hmm.calculate_forward_alignment_prob();
        random_prob = random_hmm.calculate_forward_alignment_prob();

        alignment_prob_file << alignment_prob<<std::endl;
        random_prob_file << random_prob<<std::endl;



        alignment_file_by_task<< align_hmm.get_annotated_x()<<std::endl;
        alignment_file_by_task<< align_hmm.get_annotated_y()<<std::endl<<std::endl;



        paths_file << align_hmm.get_annotated_state_path()<<std::endl;


        log_odds = log2((alignment_prob/random_prob));
        if (log_odds>0){
            log_ratios_file<< log_odds<<"\trelated "<<std::endl;
            alignment_file<< log_odds<<"\trelated "<<std::endl;
        }
        else{
            log_ratios_file<< log_odds<<"\tunrelated "<<std::endl;
            alignment_file<< log_odds<<"\tunrelated "<<std::endl;
        }
        alignment_file<< align_hmm.get_annotated_state_path()<<std::endl;
        alignment_file<< align_hmm.get_annotated_x()<<std::endl;
        alignment_file<< align_hmm.get_annotated_y()<<std::endl<<std::endl<<std::endl;


        //std::cout<<"X: "<<s1<<" "<<s1.length()<<'\n'<<"Y: "<<s2<<" "<<s2.length()<<std::endl;
    }

    alignment_file.close();
    paths_file.close();
    alignment_prob_file.close();
    random_prob_file.close();
//
//
//    s1 ="AAAAAAATTTCC";
//    s2 = "AATTTTCCGG";
//
//    align_hmm.set_observations_x(s1);
//    align_hmm.set_observations_y(s2);
//    //align_hmm.test_public_call();
//
//    //align_hmm.set_model_name("Alignment");
//    align_hmm.set_model_name("Random");
//    double a = align_hmm.calculate_forward_alignment_prob();
//    double f = align_hmm.calculate_viterbi_alignment ();
//
//
//
//    std::cout<<"X: "<<s1<<" "<<s1.length()<<'\n'<<"Y: "<<s2<<" "<<s2.length()<<std::endl;
    return 0;
}
