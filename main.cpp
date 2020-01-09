#include <iostream>
#include "constants.h"
#include "helper_functions.h"
#include "Pair_HMM.h"
int main() {
    std::cout << "Hello, World!" << std::endl;

    int i = 0, j = 0;

    // reading emission and transitions matrices from files.

    float *emissions = read_matrix_file(PATH_EMISSIONS_FILE, 5, 25);
    float *transitions = read_matrix_file(PATH_TRANSITIONS_FILE,5, 5);

    print_matrix(emissions,5,25);
    print_matrix(transitions,5,5);

    Pair_HMM A(5,25,transitions,emissions);


    //A.set_observations("ACTGTAAAA");
    //A.set_observations("AAAABBBBAA");



    std::string s1 = "ATATATATATA";
    std::string s2 = "TATATATATAT" ;
    A.set_observations_x(s1);
    A.set_observations_y(s2);


    //A.test_public_call();
    float f = A.calculate_viterbi_alignment ();
    std::cout<<"X: "<<s1<<" "<<s1.length()<<'\n'<<"Y: "<<s2<<" "<<s2.length()<<std::endl;
    return 0;
}
