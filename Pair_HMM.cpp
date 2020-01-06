//
// Created by esemenc on 1/6/20.
//

//
// Created by esemenc on 12/12/19.
//


#include "Pair_HMM.h"


Pair_HMM::Pair_HMM(int n_states, int n_observables,
         float *transitions, float *emissions) {
    _transitions = transitions;
    _emissions = emissions;
    _n_states = n_states;
    _n_observables = n_observables;


    _observations = NULL;

    _sequence_x = NULL;
    _sequence_y = NULL;

    _n_x = 0;
    _n_y = 0;

    _state_readings = new int[_n_states][2];
    calculate_states_readings();




}
void Pair_HMM::public_call()
{
    float a = get_emission_proba(3,2,2) ;
    std::cout<<"test putput "<<a<<std::endl;
}
Pair_HMM::~Pair_HMM() {
    delete[] _transitions;
    delete[] _emissions;
    if (_observations) { delete[] _observations;}
}

// this function looks at emissions matrix and save what each state reads from sequence to
//  _state_readings [state index] [1-x seq 2-y ] = 0/1 reads or not from x/y

void Pair_HMM::calculate_states_readings(){
    for (int k = 0; k < _n_states; k++) {
        int i = 0;
        bool counted_two = false;
        bool counted_x = false;
        bool counted_y = false;
        while (i < 16 and not counted_two) { // check 0-15 emission matrix options
            if (_emissions[k * _n_observables + i] > 0)
                counted_two = true;
            i++;
        }
        if (counted_two){
            this->_state_readings[k][0] = 1;
            this->_state_readings[k][1] = 1;
            std::cout<<k<<" is pair-read state "<< std::endl;

        }
        else { //16 17 18 19
            while (i < 20 and not(counted_x)) {
                if (_emissions[k * _n_observables + i] > 0)
                    counted_x = true;
                i++;
            }
            if (counted_x) {
                this->_state_readings[k][0] = 1;
                this->_state_readings[k][1] = 0;
                std::cout<<k<<" is x- emit state "<< std::endl;

            }
            else{ // 20 21 22 23
                while (i < 24 and not(counted_y)) {
                    if (_emissions[k * _n_observables + i] > 0)
                        counted_y = true;
                    i++;
                }
                if (counted_y) {
                    this->_state_readings[k][0] = 0;
                    this->_state_readings[k][1] = 1;
                    std::cout<<k<<" is y- emit state "<< std::endl;

                }
                else{
                    this->_state_readings[k][0] = 0;
                    this->_state_readings[k][1] = 0;
                    std::cout<<k<<" is silent state "<< std::endl;

                }
            }
        }

    }


}

void Pair_HMM::set_observations(int *observations, int n_observations) {


    _observations = NULL;
    _observations = observations;
    _n_observations = n_observations;

}

void Pair_HMM::set_observations_x(std::string observations) {

    delete[] _sequence_x;
    _sequence_x = NULL;

    this->_n_x = observations.length();
    _sequence_x = new char[_n_x];
    for (int i = 0; i < _n_x; i++) {
        _sequence_x[i] = observations[i];
    }
}

void Pair_HMM::set_observations_y(std::string observations) {

    delete[] _sequence_y;
    _sequence_y = NULL;

    this->_n_y = observations.length();
    _sequence_y = new char[_n_y];
    for (int i = 0; i < _n_y; i++) {
        _sequence_y[i] = observations[i];
    }

}

int Pair_HMM::delta_x(int state){

    if (_state_readings[state][0] == 1){
        return 1;
    }
    else
        return 0;


}
int Pair_HMM::delta_y(int state){

    if (_state_readings[state][1] == 1){
        return 1;
    }
    else
        return 0;

}
int Pair_HMM::get_character_index(char character){

    int index;
    if (character == 'A') { index = 0; }
    if (character == 'T') { index = 1; }
    if (character == 'G') { index = 2; }
    if (character == 'C') { index = 3; }

    return index;
    }
float Pair_HMM::get_emission_proba(int state,int i,int j){

    // 0:silent
    // 1:emit x
    // 2:emit y
    int x_ind = get_character_index(_sequence_x[i]);
    int y_ind = get_character_index(_sequence_y[j]);
    std::cout<<_sequence_x<<std::endl;
    std::cout<<_sequence_y<<std::endl;

    std::cout<<_sequence_x[i]<<std::endl;
    std::cout<<_sequence_y[j]<<std::endl;

    if (_state_readings[state][0] == 1 and  _state_readings[state][1] == 1){
        //pair-read  state

        int emission_matrix_ind = x_ind*4+y_ind;
        std::cout<<"pair read "<<emission_matrix_ind<<std::endl;
        return _emissions[state * _n_observables + emission_matrix_ind];

    }
    if (_state_readings[state][0] == 1 and  _state_readings[state][1] == 0){
        //emit x
        int emission_matrix_ind = 16 + x_ind;
        std::cout<<"X read "<<emission_matrix_ind<<std::endl;
        return _emissions[state * _n_observables + emission_matrix_ind];

    }
    if (_state_readings[state][0] == 0 and  _state_readings[state][1] == 1){
        //emit y
        int emission_matrix_ind = 20 + y_ind;
        std::cout<<"Y read "<<emission_matrix_ind<<std::endl;
        return _emissions[state * _n_observables + emission_matrix_ind];

    }
    if (_state_readings[state][0] == 0 and  _state_readings[state][1] == 0){
        //emit y
        int emission_matrix_ind = 24;
        return _emissions[state * _n_observables + emission_matrix_ind];
    }

}


int *Pair_HMM::calculate_forward_alignment() {

    int k, m, t,i,j,z;

    //Forward Aln

    float *viterbi = new float[_n_x * _n_y * _n_states];
    int *pointers = new int[_n_x * _n_y * _n_states];

    // Initialization step

    for (k = 0; k < _n_states; k++) {
        for (i = 0; i < _n_x; i++) {
            for (j = 0; j < _n_y; j++) {
                viterbi[k * _n_x * _n_y + _n_x * j + i] = 0;
                pointers[k * _n_x * _n_y + _n_x * j + i] = -1;
            }
        }
    }
    viterbi[0 *_n_x*_n_y + _n_x* 0 + 0] = 1;

    //   Calculation

    float path_maximum_likelihood, previous_state_likelihood;
    int path_previous_state_pointer, observation;


    for (i = 0; i < _n_x; i++) {
        for (j = 0; j < _n_y; j++) {
            if (not(i==0 and j==0 )){
                for (k = 1; k < _n_states; k++){
                    viterbi[k * _n_x * _n_y + _n_x * j + i] = 0;
                }
            }

            pointers[k * _n_x * _n_y + _n_x * j + i] = -1;
        }
    }


    for (t = 1; t < _n_observations; t++) {
        observation = _observations[t]; //will serve as index in emissions matrix.
        for (k = 0; k < _n_states; k++) {

            //we calculate Maximum log L using all transitions from previous states
            path_maximum_likelihood = 0;
            path_previous_state_pointer = -2; //set to -2 to trace errors in filling pointer matrix

            for (m = 0; m < _n_states; m++) {
                previous_state_likelihood =
                        viterbi[m * _n_observations + (t - 1)] + log10(_transitions[m * _n_states + k]);
                if (previous_state_likelihood > path_maximum_likelihood) {
                    path_maximum_likelihood = previous_state_likelihood;
                    path_previous_state_pointer = m;
                }
            }
            viterbi[k * _n_observations + t] =
                    log10(_emissions[k * _n_observables + (observation - 1)]) + path_maximum_likelihood;
            pointers[k * _n_observations + t] = path_previous_state_pointer;
        }
    }

    // Termination

    float best_path_likelihood = 0;
    float currentEndPathLikelihood;
    int end_state = _n_states - 1; //will serve as index in transitions matrix. for readability
    int last_predicted_state = -1;

    //find max-l state and calculate transition to the end
    for (k = 0; k < _n_states; k++) {
        currentEndPathLikelihood = viterbi[k * _n_observations + (_n_observations - 1)] +
                                   log10(_transitions[k * _n_states + end_state]);
        if (currentEndPathLikelihood > best_path_likelihood) {
            best_path_likelihood = currentEndPathLikelihood;
            path_previous_state_pointer = pointers[k * _n_observations + (_n_observations - 1)];
            last_predicted_state = k;
        }
    }

    //Backtracking the states from the last to the "to start state" pointer

    int *reconstructed_states = new int[_n_observations + 1];
    reconstructed_states[_n_observations] = 3; //set the last element to end state
    for (t = _n_observations - 1; t >= 0; --t) {
        reconstructed_states[t] = path_previous_state_pointer;
        path_previous_state_pointer = pointers[path_previous_state_pointer * _n_observations + t];
    }

    std::cout << "best likelihood " << best_path_likelihood << " last_predicted_state " << last_predicted_state << std::endl;
    std::cout << path_previous_state_pointer << " <- the last backtrack pointer " << std::endl;
    std::cout << "Viterbi alg calculation completed\n";

    delete[] viterbi;
    delete[] pointers;

    return reconstructed_states;
}



int *Pair_HMM::calculate_viterbi_state_path() {

    int k, m, t;
    //Viterbi Alg
    float *viterbi = new float[_n_observations * _n_states];
    int *pointers = new int[_n_observations * _n_states];
    float negative_infinity = -1 * std::numeric_limits<float>::infinity();

    // Initialization step

    viterbi[0 * _n_observations + 0] = log10(1);
    pointers[0 * _n_observations + 0] = -1;
    for (k = 1; k < _n_states; k++) {
        viterbi[k * _n_observations + 0] = negative_infinity;
        pointers[k * _n_observations + 0] = -1;
    }

    //   Calculation

    float path_maximum_likelihood, previous_state_likelihood;
    int path_previous_state_pointer, observation;

    for (t = 1; t < _n_observations; t++) {
        observation = _observations[t]; //will serve as index in emissions matrix.
        for (k = 0; k < _n_states; k++) {

            //we calculate Maximum log L using all transitions from previous states
            path_maximum_likelihood = negative_infinity;
            path_previous_state_pointer = -2; //set to -2 to trace errors in filling pointer matrix

            for (m = 0; m < _n_states; m++) {
                previous_state_likelihood =
                        viterbi[m * _n_observations + (t - 1)] + log10(_transitions[m * _n_states + k]);
                if (previous_state_likelihood > path_maximum_likelihood) {
                    path_maximum_likelihood = previous_state_likelihood;
                    path_previous_state_pointer = m;
                }
            }
            viterbi[k * _n_observations + t] =
                    log10(_emissions[k * _n_observables + (observation - 1)]) + path_maximum_likelihood;
            pointers[k * _n_observations + t] = path_previous_state_pointer;
        }
    }

    // Termination

    float best_path_likelihood = negative_infinity;
    float currentEndPathLikelihood;
    int end_state = _n_states - 1; //will serve as index in transitions matrix. for readability
    int last_predicted_state = -1;

    //find max-l state and calculate transition to the end
    for (k = 0; k < _n_states; k++) {
        currentEndPathLikelihood = viterbi[k * _n_observations + (_n_observations - 1)] +
                                   log10(_transitions[k * _n_states + end_state]);
        if (currentEndPathLikelihood > best_path_likelihood) {
            best_path_likelihood = currentEndPathLikelihood;
            path_previous_state_pointer = pointers[k * _n_observations + (_n_observations - 1)];
            last_predicted_state = k;
        }
    }

    //Backtracking the states from the last to the "to start state" pointer

    int *reconstructed_states = new int[_n_observations + 1];
    reconstructed_states[_n_observations] = 3; //set the last element to end state
    for (t = _n_observations - 1; t >= 0; --t) {
        reconstructed_states[t] = path_previous_state_pointer;
        path_previous_state_pointer = pointers[path_previous_state_pointer * _n_observations + t];
    }

    std::cout << "best likelihood " << best_path_likelihood << " last_predicted_state " << last_predicted_state << std::endl;
    std::cout << path_previous_state_pointer << " <- the last backtrack pointer " << std::endl;
    std::cout << "Viterbi alg calculation completed\n";

    delete[] viterbi;
    delete[] pointers;

    return reconstructed_states;
}
