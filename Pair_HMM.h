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
    Pair_HMM(int n_states, int n_observables, double *transitions, double *emissions);
    ~Pair_HMM();

    void set_observations_x(const std::string& observations);
    void set_observations_y(const std::string& observations);
    void set_model_name(const std::string& model_name);

    double calculate_viterbi_alignment();
    double calculate_forward_alignment_prob();

    std::string get_annotated_x();
    std::string get_annotated_y();
    std::string get_annotated_state_path();


private:

    double *_transitions;
    double *_emissions;
    int _n_states, _n_observables; // number of states, amount of possible emissions.

    int _n_x;
    int _n_y;
    char *_sequence_x;
    char *_sequence_y;

    std::string _aligned_x;
    std::string _aligned_y;
    std::string _state_path;
    std::string _model_name;

    // here the properties of states saved (Exy or Ex or Ey or silent). init in calculate_state readings()
    int (*_state_readings)[2];

    int delta_x(int state);
    int delta_y(int state);

    int m_index(int x, int y, int z);

    //rewrite if alphabet and emission matrix structure have changed.
    int get_character_index(char character);
    double get_emission_proba(int state,int i,int j);
    void calculate_states_readings(); // save what each state reads

};