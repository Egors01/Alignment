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


    A.set_observations_x("AAGGGT");
    A.set_observations_y("AATAAGCGCGCGCGA");

    A.public_call();

    return 0;
}
