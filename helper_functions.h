//
// Created by esemenc on 1/6/20.
//

#ifndef ALIGNMENT_HELPER_FUNCTIONS_H
#define ALIGNMENT_HELPER_FUNCTIONS_H


#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#ifndef VITERBI_HOMEWORK_HELPER_FUNCTIONS_H
#define VITERBI_HOMEWORK_HELPER_FUNCTIONS_H

#endif //VITERBI_HOMEWORK_HELPER_FUNCTIONS_H



float *read_matrix_file(std::string filename, int N, int M);

void print_matrix(float *matrix, int N, int M);

std::string *sequences_reader(std::string filename, int number_of_sequences);

std::string array_to_string(int *arr, int n_elem);

int *string_to_arr(std::string s_sequence);

std::string annotate_path(std::string reconstructed_states) ;




#endif //ALIGNMENT_HELPER_FUNCTIONS_H
