//
// Created by esemenc on 1/6/20.
//

#ifndef ALIGNMENT_HMM_H
#define ALIGNMENT_HMM_H

#endif //ALIGNMENT_HMM_H
#include <iostream>
#include <cmath>
#include <limits>
#include <string>
#include <cstring>
// Class HMM takes transition and emission matrices and accept the int array of observations
// Viterbi state path reconstruction implementeted as a member function of this class.
// Other helper functions (array to string conversions, file readers) declared and implemented below int* calculate_viterbi_state_path()

// We assume that emission observable reflects its X-> pos in emission matrix.
// Observations should be already encoded to int and correspond to the index in the emission matrix



class Pair_HMM {

public:
    Pair_HMM(int n_states, int n_observables, float *transitions, float *emissions);

    void set_observations(std::string observation);
    void set_observations_x(const std::string& observations);
    void set_observations_y(const std::string& observations);
    void set_model_name(std::string model_name);

    int *calculate_viterbi_state_path();
    float calculate_forward_alignment_prob();

    ~Pair_HMM();

    void test_public_call();
    float calculate_viterbi_alignment();

private:

    float *_transitions;
    float *_emissions;

    int *_observations;
    int _n_observations;
    int _n_states, _n_observables; // number of states, amount of possible emissions.

    char *_sequence_test;

    int _n_x;
    int _n_y;
    char *_sequence_x;
    char *_sequence_y;

    std::string aligned_x;
    std::string aligned_y;
    std::string state_path;
    std::string _model_name;


    int (*_state_readings)[2];

    int delta_x(int state);
    int delta_y(int state);

    int m_index(int x, int y, int z);
    static int get_character_index(char character);

    float get_emission_proba(int state,int i,int j);
    void calculate_states_readings();

};