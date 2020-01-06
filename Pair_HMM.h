//
// Created by esemenc on 1/6/20.
//

#ifndef ALIGNMENT_HMM_H
#define ALIGNMENT_HMM_H

#endif //ALIGNMENT_HMM_H
#include <iostream>
#include <cmath>
#include <limits>


// Class HMM takes transition and emission matrices and accept the int array of observations
// Viterbi state path reconstruction implementeted as a member function of this class.
// Other helper functions (array to string conversions, file readers) declared and implemented below int* calculate_viterbi_state_path()

// We assume that emission observable reflects its X-> pos in emission matrix.
// Observations should be already encoded to int and correspond to the index in the emission matrix



class Pair_HMM {
public:
    Pair_HMM(int n_states, int n_observables, float *transitions, float *emissions);

    void set_observations(int *observations, int n_observations);

    void set_observations_x(std::string observations);
    void set_observations_y(std::string observations);

    int *calculate_viterbi_state_path();


    int *calculate_forward_alignment();

    ~Pair_HMM();

    int *_observations;
    int _n_observations;

    char *_sequence_x = NULL;
    char *_sequence_y = NULL;

    int _n_x;
    int _n_y;
    void public_call();
private:

    float *_transitions;
    float *_emissions;
    int (*_state_readings)[2];
    void calculate_states_readings();

    int delta_x(int state);
    int delta_y(int state);
    float get_emission_proba(int state,int i,int j);
    int _n_states, _n_observables; // number of states, amount of possible emissions.


    int get_character_index(char character);


};
