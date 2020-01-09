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

    int index = 777;

    if (character == 'A') {  index = 0; }
    if (character == 'T') {  index = 1; }
    if (character == 'G') {  index = 2; }
    if (character == 'C') {  index = 3; }
    if (index == 777){
        std::cout<<"Attempt to read not-specified character :"<<character<<std::endl;
        throw std::invalid_argument( "received wrong index" );
    }
    else
        return index;
    }
float Pair_HMM::get_emission_proba(int state,int i,int j){

//   std::cout<<_sequence_x<<std::endl;
//    std::cout<<_sequence_y<<std::endl;

//    std::cout<<_sequence_x[i]<<std::endl;
//    std::cout<<_sequence_y[j]<<std::endl;

    // 0 0 :silent
    // 1 1 :emit x
    // 2 2 :emit y


    if (_state_readings[state][0] == 1 and  _state_readings[state][1] == 1){
        //pair-read  state

        int x_ind = get_character_index(_sequence_x[i]);
        int y_ind = get_character_index(_sequence_y[j]);
        if (i<0 or j<0){
            std::cout<<"Attempt to get emission for unexisting pos i,j: "<<i<<", "<<j<<std::endl;
            throw std::invalid_argument( "received wrong position" );
        }

        int emission_matrix_ind = x_ind*4 + y_ind;

        //std::cout<<"pair read "<<emission_matrix_ind<<" "<<_emissions[state * _n_observables + emission_matrix_ind]<<" "<<state<<std::endl;
        //std::cout<<"x: "<<_sequence_x[i]<<" "<<i<<std::endl;
        //std::cout<<"y: "<<_sequence_y[j]<<" "<<j<<std::endl;

        return _emissions[state * _n_observables + emission_matrix_ind];

    }
    if (_state_readings[state][0] == 1 and  _state_readings[state][1] == 0){
        //emit x

        int x_ind = get_character_index(_sequence_x[i]);

        if (i<0){
            std::cout<<"Attempt to get emission for unexisting pos Ex i,j: "<<i<<", "<<j<<std::endl;
            throw std::invalid_argument( "received wrong position" );
        }


        int emission_matrix_ind = 16 + x_ind;
       // std::cout<<"X read "<<emission_matrix_ind<<" "<<_emissions[state * _n_observables + emission_matrix_ind]<<" "<<state<<std::endl;
        //std::cout<<"x: "<<_sequence_x[i]<<" "<<i<<std::endl;
       // std::cout<<"y: "<<_sequence_y[j]<<" "<<j<<std::endl;
        return _emissions[state * _n_observables + emission_matrix_ind];

    }
    if (_state_readings[state][0] == 0 and  _state_readings[state][1] == 1){
        //emit y
        int y_ind = get_character_index(_sequence_y[j]);
        int emission_matrix_ind = 20 + y_ind;

        if (j<0){
            std::cout<<"Attempt to get emission for unexisting pos Ey i,j: "<<i<<", "<<j<<std::endl;
            throw std::invalid_argument( "received wrong position" );
        }
       // std::cout<<"Y read "<<emission_matrix_ind<<" "<<_emissions[state * _n_observables + emission_matrix_ind]<<" "<<state<<std::endl;
      //  std::cout<<"x: "<<_sequence_x[i]<<" "<<i<<std::endl;
      //  std::cout<<"y: "<<_sequence_y[j]<<" "<<j<<std::endl;
        return _emissions[state * _n_observables + emission_matrix_ind];

    }
    if (_state_readings[state][0] == 0 and  _state_readings[state][1] == 0){
       // std::cout<<"silent state "<<state<<std::endl;
       // std::cout<<"x: "<<_sequence_x[i]<<" "<<i<<std::endl;
       // std::cout<<"y: "<<_sequence_y[j]<<" "<<j<<std::endl;
        int emission_matrix_ind = 24;
        return _emissions[state * _n_observables + emission_matrix_ind];
    }

}

int Pair_HMM::m_index(int x, int y, int z){
    //return ((x) + (y) *_n_x + (z) *_n_x*_n_y);
    return ((x) + (y) *(_n_x+1) + (z) *(_n_x+1)*(_n_y+1));
}

float Pair_HMM::calculate_forward_alignment() {

    int k, m, t,i,j,z;

    //Forward Aln

    float *mforward = new float[_n_x * _n_y * _n_states];

    std::cout<<_sequence_x<<std::endl;
    std::cout<<_sequence_y<<std::endl;

    // Initialization step
    std::cout << "Forward initialization  ...  "<< std::endl;

    for (k = 0; k < _n_states; k++) {
        for (i = 0; i < _n_x; i++) {
            for (j = 0; j < _n_y; j++) {
                mforward[m_index(i,j,k)] = 0;
                //pointers[k * _n_x * _n_y + _n_x * j + i] = -1;
            }
        }
    }
    mforward[0 * _n_x * _n_y + _n_x * 0 + 0] = 1;

    //   Calculation
    std::cout << "Forward calculation ...  "<< std::endl;
    float add,mf_prev,trans;
    int i_prev,j_prev;
    for (i = 0; i < _n_x; i++) {
        for (j = 0; j < _n_y; j++) {
            if (not(i==0 and j==0 )){
                for (k = 1; k < _n_states-1; k++){
                    for (m = 0;m < _n_states-1; m++){

                        i_prev = i-delta_x(k);
                        j_prev = j-delta_y(k);

                        if (not (i_prev<0 or j_prev<0)){
                            mforward[m_index(i,j,k)]+= mforward[ m_index(i_prev,j_prev,m)] *_transitions[m*_n_states + k ]*get_emission_proba(k,i,j);
                        }
                    }

                }
            }
        }
    }

    // termination
    std::cout << "Forward termination  ...  "<< std::endl;

    for (m = 1;m < _n_states; m++) {
        mforward[m_index(_n_x - 1, _n_y - 1, _n_states - 1)] += mforward[m_index(_n_x - 1, _n_y - 1, m)] *_transitions[m * _n_states + (_n_states - 1)];

    }

    float res_ml  = mforward[m_index(_n_x-1,_n_y-1,_n_states-1)];
    std::cout << "Likelihood in forward alg P(X, Y |alignment model) "<< res_ml << std::endl;


    for (k = 0; k < _n_states ; k++) {
        std::cout <<"--------"<<k<<"---------" <<std::endl;
        for (i = 0; i < _n_x; i++) {
            for (j = 0; j < _n_y; j++) {
                std::cout << mforward[m_index(i,j,k)]<< " \t " << std::setprecision(3);
                }
            std::cout << std::endl;
            }

        }
    delete[] mforward;


    return res_ml;
}

float Pair_HMM::calculate_viterbi_alignment() {

    int k, m, t,i,j,z;

    //Forward Aln

    float *mviterbi = new float[(_n_x+1) * (_n_y+1) * _n_states];
    int *pointers = new int[(_n_x+1) * (_n_y+1) * _n_states];

    std::cout<<_sequence_x<<std::endl;
    std::cout<<_sequence_y<<std::endl;

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
    float path_elem = 0 ,ml_path_proba,trans,vit;
    int i_prev,j_prev,arg_max_previous_state;
    for (i = 0; i <= _n_x; i++) {
        for (j = 0; j <= _n_y; j++) {
            if (not(i==0 and j==0 ))
            {
                for (k = 1; k < _n_states; k++){
                    ml_path_proba = 0;
                    arg_max_previous_state = 1000*i+100*j+k;

                    for (m = 0;m < _n_states-1; m++){

                        i_prev = i-delta_x(k);
                        j_prev = j-delta_y(k);


                        trans = _transitions[m * _n_states + k ];

                        if (i_prev>=0 and j_prev>=0){

                            vit = mviterbi[ m_index(i_prev, j_prev, m)];
                            path_elem = mviterbi[ m_index(i_prev, j_prev, m)] * _transitions[m * _n_states + k ];


                            if ( ((i-1<0) and (delta_x(k)!=0)) or ((j-1<0) and (delta_y(k)!=0))){
                                std::cout<<"Error index\n";
                            }
                            if (path_elem > ml_path_proba){
                                ml_path_proba = path_elem * get_emission_proba(k, i-1, j-1);;
                                arg_max_previous_state = m;
                            }
                            std::cout<<" i,j "<<i<<" "<<j<<" for "<<m<<" ->  "<<k<<"\t trans "<<trans<<" vit["<<i_prev<<"]["<<j_prev<<"]["<<m<<"] \t"<<vit <<std::endl<< std::setprecision(3);


                        }
                        else{
                            path_elem = 0;
                            std::cout<<"Set zero emisssion i,j "<<i<<" "<<j<<" for "<<m<<" ->  "<<k<<"\t trans "<<trans<<"\t i_pr, j_pr "<< i_prev<<" "<< j_prev<<" " <<std::endl<< std::setprecision(3);
                        }

                    }
                    mviterbi[m_index(i, j, k)] = ml_path_proba;
                    pointers[m_index(i, j, k)] = arg_max_previous_state;


                    if (arg_max_previous_state ==  1000*i+100*j+k){
                        std::cout<<"-----\nSet vit["<<i<<"]["<<j<<"]["<<k<<"]:=  "<<ml_path_proba<<"p:=  "<<arg_max_previous_state<<" ->  "<<k<<"\t trans "<<trans<<"\t; pointer "<<arg_max_previous_state<<"i_pr, j_pr "<< i_prev<<" "<< j_prev<<" NOT CALC \n\n" <<std::endl<< std::setprecision(3);

                        //std::cout<<" not calculated i,j "<<i<<" "<<j<<" for -> "<<k<<" "<<std::endl;
                    }else{
                        std::cout<<"-----\nSet vit["<<i<<"]["<<j<<"]["<<k<<"]:= "<<ml_path_proba<<" p:= "<<arg_max_previous_state<<" ->  "<<k<<"\t trans "<<trans<<"\t; pointer "<<arg_max_previous_state<<"i_pr, j_pr "<< i_prev<<" "<< j_prev<<" \n\n" <<std::endl<< std::setprecision(3);

                    }

                }
            }
        }
    }

    // termination
    std::cout << "Viterbi termination  ...  "<< std::endl;
    ml_path_proba = 0;
    arg_max_previous_state = -1;
    for (m = 1;m < _n_states-1; m++) {
        path_elem = mviterbi[m_index(_n_x, _n_y, m)] * _transitions[m * _n_states + (_n_states - 1)];
        std::cout<<" TERMINATION "<<m<<" ->  "<<_n_states-1<<"\t trans "<< _transitions[m * _n_states + (_n_states - 1)]<<" vit["<<_n_x<<"]["<<_n_y<<"]["<<m<<"] \t"<<mviterbi[m_index(_n_x, _n_y, m)] <<std::endl<< std::setprecision(3);

        if (path_elem > ml_path_proba){
            ml_path_proba = path_elem;
            arg_max_previous_state = m;
            std::cout<<"^^^^ set this max\n";
        }
    }
    mviterbi[m_index(_n_x, _n_y, _n_states-1)] = ml_path_proba;
    pointers[m_index(_n_x, _n_y, _n_states-1)] = arg_max_previous_state;

    if (ml_path_proba==0) {
        std::cout << "did not set ml end  " << std::endl;
    }

#ifdef PRINT
    for (k = 0; k < _n_states; k++) {
        std::cout << "--------" << k << "---------" << std::endl;
        for (i = 0; i <= _n_x; i++) {
            for (j = 0; j <= _n_y; j++) {
                std::cout << mviterbi[m_index(i, j, k)] << " \t " << std::setprecision(3);
            }
            std::cout << std::endl;
        }

    }

    std::cout << "--------pointers--------" << std::endl;

    for (k = 0; k < _n_states; k++) {
        std::cout << "--------" << k << "---------" << std::endl;
        for (i = 0; i <= _n_x; i++) {
            for (j = 0; j <= _n_y; j++) {
                std::cout << pointers[m_index(i, j, k)] << " \t ";
            }
            std::cout << std::endl;
        }

    }
#endif



    // Traceback
    std::cout << "Viterbi traceback  ...  "<< std::endl;
    int next_state_pointer = arg_max_previous_state;

    int ix=0,iy = 0;
    char *annotated_x = new char[_n_x+_n_y];
    char *annotated_y = new char[_n_x+_n_y];

    int  *states_path = new int[_n_x+_n_y];
    int s = 0;
    states_path[s] = next_state_pointer;
    int current_state = next_state_pointer;
    i = _n_x;
    j = _n_y;

    while(next_state_pointer != -1 ){

        current_state = next_state_pointer;
        states_path[s] = current_state ;


        if (_state_readings[current_state][0] == 1 and _state_readings[current_state][1] == 1){
            // pair read
            std::cout << "E XY: [" << i << "][" << j << "] " << current_state << "\t " << _sequence_x[i-1] << " " << _sequence_y[j-1] << "  s(" << i - 1 << ")(" << j - 1 << ") ";
            annotated_x[ix] = _sequence_x[i-1];
            annotated_y[iy] = _sequence_y[j-1];
        }
        if (_state_readings[current_state][0] == 1 and _state_readings[current_state][1] == 0){
            // X read
            std::cout << "E  X: [" << i << "][" << j << "] " << current_state << "\t " << _sequence_x[i - 1] << " " << '-' << "  s(" << i - 1 << ")(" << " " << ") ";
            annotated_x[ix] = _sequence_x[i-1];
            annotated_y[iy] = '-';
        }
        if (_state_readings[current_state][0] == 0 and _state_readings[current_state][1] == 1){
            // Y read
            std::cout << "E  Y: [" << i << "][" << j << "] " << current_state << "\t " << '-' << " " << _sequence_y[j - 1] << "  s(" << " " << ")(" << j - 1 << ") ";
            annotated_x[ix] = '-';
            annotated_y[iy] = _sequence_y[j-1];
        }
        //silent without change
        // --

        // get new pointer
        next_state_pointer = pointers[m_index(i, j, current_state)];
        i = i - delta_x(current_state);
        j = j - delta_y(current_state);
        ix++; iy++; s++;
        std::cout << "\t > New: [" << i << "][" << j << "]["<< current_state <<"] :=" << next_state_pointer << std::endl;

    }
    std::cout << "\t > Final : [" << i << "][" << j << "] " << next_state_pointer << std::endl;
    int n_states_al = s;
    for (j = n_states_al-2; j>=0 ; j--){
        std::cout << states_path[j];
    }
    std::cout<< std::endl;


    int nx_a = ix;
    int ny_a = iy;

    for (i = nx_a-2;i >=0;i--){

        std::cout << annotated_x[i];
    }
    std::cout <<std::endl;

    for (j = ny_a-2; j>=0 ; j--){
        std::cout << annotated_y[j];
    }
    std::cout <<std::endl;

    for (j = 0;j< n_states_al; j++){
        std::cout << states_path[j];
    }
    std::cout<< std::endl;

    float res_ml  = mviterbi[m_index(_n_x, _n_y, _n_states - 1)];
    std::cout << "Likelihood in viterbi alg P(X, Y |alignment model) "<< res_ml << std::endl;




    delete[] mviterbi;

    return res_ml;
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
