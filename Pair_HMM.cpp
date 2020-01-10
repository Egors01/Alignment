//
// Created by esemenc on 1/6/20.
//


#include "Pair_HMM.h"
#include "helper_functions.h"
#define  PRINT

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

    _model_name = "Alignment model (default)";
    _state_readings = new int[_n_states][2];
    calculate_states_readings();




}
void Pair_HMM::test_public_call()
{
    float a = get_emission_proba(3,2,2) ;
    std::cout<<"test putput "<<a<<std::endl;
}
Pair_HMM::~Pair_HMM() {
    delete[] _transitions;
    delete[] _emissions;

}



void Pair_HMM::calculate_states_readings(){

    // this function called once the emission matrix is read.
    // it looks in emissions matrix and saves if each state reads from X or Y sequence.
    // It saves reads from sequence to _state_readings
    // Acess: _state_readings[STATE][ 0|1  ( X or Y) ] = 0|1 (state reads| doesnt read)

    // Indexes in emission matrix X- > 0..25
    // 0 ('A', 'A') //9 ('G', 'T') //18 ('G', '-')
    //1 ('A', 'T')  //10 ('G', 'G') //19 ('C', '-')
    //2 ('A', 'G') //11 ('G', 'C') //20 ('-', 'A')
    //3 ('A', 'C') //12 ('C', 'A') //21 ('-', 'T')
    //4 ('T', 'A') //13 ('C', 'T') //22 ('-', 'G')
    //5 ('T', 'T') //14 ('C', 'G') //23 ('-', 'C')
    //6 ('T', 'G') //15 ('C', 'C') //24 ('-', '-')
    //7 ('T', 'C') //16 ('A', '-')
    //8 ('G', 'A') //17 ('T', '-')
    // States Y 0..N_STATES

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

void Pair_HMM::set_observations(std::string observations) {

    const char* seq_array = observations.c_str();
    for (int i = 0; i < observations.length(); i++) {
        //std::cout<<seq_array[i]<<std::endl;
    }

    _sequence_test = new char [observations.length()+1];
    std::strcpy (_sequence_test, observations.c_str());
    for (int i = 0; i < observations.length(); i++) {
        std::cout<<_sequence_test[i];
    }
    std::cout<<std::endl;


}

void Pair_HMM::set_observations_x(const std::string& observations) {

    _n_x = observations.length();
    _sequence_x = new char [observations.length()+1];
    std::strcpy (_sequence_x, observations.c_str());

}

void Pair_HMM::set_observations_y(const std::string& observations) {

    _n_y = observations.length();
    _sequence_y = new char [observations.length()+1];
    std::strcpy (_sequence_y, observations.c_str());

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

    int index = -777;

    if (character == 'A') {  index = 0; }
    if (character == 'T') {  index = 1; }
    if (character == 'G') {  index = 2; }
    if (character == 'C') {  index = 3; }
    if (index == -777){
        std::cout<<"Attempt to read not-specified character :"<<character<<std::endl;
        throw std::invalid_argument( "received wrong index" );
    }
    else
        return index;
    }
float Pair_HMM::get_emission_proba(int state,int i,int j){

    // Takes emission prob from the emission matrix by index of the i,j characters

    //  emission function should ignore  negative indexes in cases
    // "first element emission" SeqX(i)|SeqY(j) (i<0;j>0;state = Ey) and (i>0;j<0;state = Ex)
    //  where this  index is irrelevant ot the state emission

    // readings X, readings Y: state_annotation
    // 0 0 :silent
    // 1 1 :emit x
    // 2 2 :emit y

    //if pair-read  state
    if (_state_readings[state][0] == 1 and  _state_readings[state][1] == 1){

        int x_ind = get_character_index(_sequence_x[i]);
        int y_ind = get_character_index(_sequence_y[j]);
        if (i<0 or j<0){
            std::cout<<"Attempt to get emission for unexisting pos i,j: "<<i<<", "<<j<<std::endl;
            throw std::invalid_argument( "received wrong position" );
        }
        int emission_matrix_ind = x_ind*4 + y_ind;
        return _emissions[state * _n_observables + emission_matrix_ind];

    }
    //if emit x
    if (_state_readings[state][0] == 1 and  _state_readings[state][1] == 0){

        int x_ind = get_character_index(_sequence_x[i]);
        //extra check for errors
        if (i<0){
            std::cout<<"Attempt to get emission for unexisting pos Ex i,j: "<<i<<", "<<j<<std::endl;
            throw std::invalid_argument( "received wrong position" );
        }

        int emission_matrix_ind = 16 + x_ind;
        return _emissions[state * _n_observables + emission_matrix_ind];

    }
    //if emit y
    if (_state_readings[state][0] == 0 and  _state_readings[state][1] == 1){

        int y_ind = get_character_index(_sequence_y[j]);
        int emission_matrix_ind = 20 + y_ind;
        //extra check for errors
        if (j<0){
            std::cout<<"Attempt to get emission for unexisting pos Ey i,j: "<<i<<", "<<j<<std::endl;
            throw std::invalid_argument( "received wrong position" );
        }
        return _emissions[state * _n_observables + emission_matrix_ind];

    }
    //if silent state
    if (_state_readings[state][0] == 0 and  _state_readings[state][1] == 0){
        if (not(state == 0 or state ==_n_states-1)){ // not and and not start
            return 1; // dont count for emission proba
        }
        else
            return  0; // zero if it is start or end

        //int emission_matrix_ind = 24;
        //return _emissions[state * _n_observables + emission_matrix_ind];
    }
}

int Pair_HMM::m_index(int x, int y, int z){
    //return ((x) + (y) *_n_x + (z) *_n_x*_n_y);
    return ((x) + (y) *(_n_x+1) + (z) *(_n_x+1)*(_n_y+1));
}

float Pair_HMM::calculate_forward_alignment_prob() {

    int k, m, t,i,j,z;

    float *mforward = new float[(_n_x + 1) * (_n_y + 1) * _n_states];

    // Initialization step
    std::cout << "Forward prob initialization  ...  "<< std::endl;

    for (k = 0; k <  _n_states; k++) {
        for (i = 0; i <=  _n_x; i++) {
            for (j = 0; j <= _n_y; j++) {
                mforward[m_index(i, j, k)] = 0;
            }
        }
    }
    mforward[m_index(0, 0, 0)] = 1;

    //   Calculation
    std::cout << "Forward calculation ...  "<< std::endl;
    float path_elem = 0 ,ml_path_proba;
    int i_prev,j_prev;


    // Matrices are 0- based where i|j =0 has meaning of "not started reading" index
    // Strings are 0- based, therefore we -1 to get emission for symbol

    for (i = 0; i <= _n_x; i++) {
        for (j = 0; j <= _n_y; j++) {
            if (not(i==0 and j==0 ))
            {
                for (k = 1; k < _n_states; k++){
                    // sum k -> m calculation
                    for (m = 0;m < _n_states-1; m++){
                        i_prev = i-delta_x(k);
                        j_prev = j-delta_y(k);
                        if (i_prev>=0 and j_prev>=0){ // check if we previous path exists
                            mforward[m_index(i, j, k)]+= mforward[ m_index(i_prev, j_prev, m)] * _transitions[m * _n_states + k ]* get_emission_proba(k, i-1, j-1);
                            std::cout << " i,j " << i << " " << j << " for " << m << " ->  " << k << "\t trans "
                                      << _transitions[m * _n_states + k ] << " f[" << i_prev << "][" << j_prev << "][" << m << "] \t" << mforward[ m_index(i_prev, j_prev, m)]
                                      << std::endl << std::setprecision(3);
                        }
                        else {
                            std::cout << "Set zero emisssion i,j " << i << " " << j << " for " << m << " ->  " << k
                                      << "\t trans " << _transitions[m * _n_states + k ] << "\t i_pr, j_pr " << i_prev << " " << j_prev << " "
                                      << std::endl << std::setprecision(3);
                        }

                    }
                    if (mforward[m_index(i, j, k)]>0){
                        std::cout << "------\n Set i,j  f[" << i << "][" << j << "][" << k << "] \t" << mforward[ m_index(i, j, k)]
                                 << std::endl << std::endl<< std::setprecision(3);
                    }
                }
            }
        }
    }

    // termination. Calculating the max(m) -> end state
    std::cout << "Forward termination  ...  "<< std::endl;
    for (m = 1;m < _n_states-1; m++) {
        mforward[m_index(_n_x, _n_y, _n_states - 1)]+=
                mforward[m_index(_n_x, _n_y, m)] * _transitions[m * _n_states + (_n_states - 1)];
         //std::cout << " path elem termination " <<path_elem<<" trans "<< _transitions[m * _n_states + (_n_states - 1)]<<std::endl;
    }

#ifdef PRINT
    for (k = 0; k < _n_states; k++) {
        std::cout << "--------" << k << "---------" << std::endl;
        for (i = 0; i <= _n_x; i++) {
            for (j = 0; j <= _n_y; j++) {
                std::cout << mforward[m_index(i, j, k)] << " \t " << std::setprecision(3);
            }
            std::cout << std::endl;
        }

    }

#endif


    float res_ml  = mforward[m_index(_n_x, _n_y, _n_states - 1)] ;
    std::cout << "Likelihood in forward alg P(X, Y |" << _model_name << " )" << res_ml << std::endl;

    delete[] mforward;


    return res_ml;
}

float Pair_HMM::calculate_viterbi_alignment() {

    int k, m, t,i,j,z;

    float *mviterbi = new float[(_n_x+1) * (_n_y+1) * _n_states];
    int *pointers = new int[(_n_x+1) * (_n_y+1) * _n_states];

    // Initialization step
    std::cout << "Viterbi initialization  ...  "<< std::endl;

    for (k = 0; k <  _n_states; k++) {
        for (i = 0; i <=  _n_x; i++) {
            for (j = 0; j <= _n_y; j++) {
                mviterbi[m_index(i, j, k)] = 0;
                pointers[m_index(i, j, k)] = -1;
            }
        }
    }
    mviterbi[m_index(0, 0, 0)] = 1;

    //   Calculation
    std::cout << "Viterbi calculation ...  "<< std::endl;
    float path_elem = 0 ,ml_path_proba;
    int i_prev,j_prev,arg_max_previous_state;


    // Matrices are 0- based where i|j =0 has meaning of "not started reading" index
    // Strings are 0- based, therefore we -1 to get emission for symbol
    for (i = 0; i <= _n_x; i++) {
        for (j = 0; j <= _n_y; j++) {
            if (not(i == 0 and j == 0)) {
                for (k = 1; k < _n_states; k++) {
                    ml_path_proba = 0;
                    arg_max_previous_state = -1;

                    // max k -> m calculation
                    for (m = 0; m < _n_states - 1; m++) {
                        i_prev = i - delta_x(k);
                        j_prev = j - delta_y(k);
                        if (i_prev >= 0 and j_prev >= 0) { // check if we previous path exists

                            float trans = _transitions[m * _n_states + k];

                            if (i_prev >= 0 and j_prev >= 0) {

                                float vit = mviterbi[m_index(i_prev, j_prev, m)];
                                path_elem = mviterbi[m_index(i_prev, j_prev, m)] * _transitions[m * _n_states + k];
                                if (path_elem > ml_path_proba) {

                                    ml_path_proba = path_elem * get_emission_proba(k, i - 1, j - 1);;
                                    arg_max_previous_state = m;
                                }
                                std::cout << " i,j " << i << " " << j << " for " << m << " ->  " << k << "\t trans "
                                          << trans << " vit[" << i_prev << "][" << j_prev << "][" << m << "] \t" << vit
                                          << std::endl << std::setprecision(3);


                            } else { //if path does not exists - setting zero emission / ignore
                                path_elem = 0;
                                std::cout << "Set zero emisssion i,j " << i << " " << j << " for " << m << " ->  " << k
                                          << "\t trans " << trans << "\t i_pr, j_pr " << i_prev << " " << j_prev << " "
                                          << std::endl << std::setprecision(3);
                            }

                        }
                        mviterbi[m_index(i, j, k)] = ml_path_proba;
                        pointers[m_index(i, j, k)] = arg_max_previous_state;
                    }

                }
            }
        }
    }
    // termination. Calculating the max(m) -> end state
    std::cout << "Viterbi termination  ...  "<< std::endl;

    ml_path_proba = 0;
    arg_max_previous_state = -1;
    for (m = 1;m < _n_states-1; m++) {
        path_elem = mviterbi[m_index(_n_x, _n_y, m)] * _transitions[m * _n_states + (_n_states - 1)];
        if (path_elem > ml_path_proba){
            ml_path_proba = path_elem;
            arg_max_previous_state = m;
        }
    }
    mviterbi[m_index(_n_x, _n_y, _n_states-1)] = ml_path_proba;
    pointers[m_index(_n_x, _n_y, _n_states-1)] = arg_max_previous_state;

    if (ml_path_proba==0) {
        std::cout << "did not set ml end  " << std::endl;
    }




    // Traceback
    std::cout << "Viterbi traceback  ...  "<< std::endl;
    int next_state_pointer = arg_max_previous_state;

    int s=0;
    char *annotated_x = new char[_n_x+_n_y];
    char *annotated_y = new char[_n_x+_n_y];

    int  *states_path = new int[_n_x+_n_y];
    states_path[s] = next_state_pointer;

    int current_state;
    i = _n_x;
    j = _n_y;
    while(next_state_pointer != -1 ){
        current_state = next_state_pointer;
        states_path[s] = current_state ;
        if (_state_readings[current_state][0] == 1 and _state_readings[current_state][1] == 1){
            // pair read
            annotated_x[s] = _sequence_x[i-1];
            annotated_y[s] = _sequence_y[j-1];
        }
        if (_state_readings[current_state][0] == 1 and _state_readings[current_state][1] == 0){
            // X read
            annotated_x[s] = _sequence_x[i-1];
            annotated_y[s] = '-';
        }
        if (_state_readings[current_state][0] == 0 and _state_readings[current_state][1] == 1){
            // Y read
            annotated_x[s] = '-';
            annotated_y[s] = _sequence_y[j-1];
        }
        //silent if we have it i.e in random model
        if (_state_readings[current_state][0] == 0 and _state_readings[current_state][1] == 0){
            // Y read
            annotated_x[s] = '-';
            annotated_y[s] = '-';
        }

        //read next pointer
        next_state_pointer = pointers[m_index(i, j, current_state)];
        // get new i, j and
        i = i - delta_x(current_state);
        j = j - delta_y(current_state);
        s++;
    }
//    aligned_x = "";
//    aligned_y = "";
//    state_path = "";
    for (i = s-2;i >=0;i--){
        aligned_x+= annotated_x[i];
        aligned_y+= annotated_y[i];
        state_path+= (char)(states_path[i]+'0') ;
    }

    std::cout<<state_path <<'\n'<< aligned_x<<'\n'<<aligned_y<<'\n';
    float res_ml = mviterbi[m_index(_n_x, _n_y, _n_states - 1)];
    std::cout << "Likelihood in viterbi alg P(X, Y |" << _model_name << " )" << res_ml << std::endl;

    delete[] annotated_x;
    delete[] annotated_y;
    delete[] states_path;
    delete[] mviterbi;
    delete[] pointers;

    return res_ml;
}



void Pair_HMM::set_model_name(std::string model_name) {
    this->_model_name = model_name;
}
