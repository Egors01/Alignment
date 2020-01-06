//
// Created by esemenc on 1/6/20.
//

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include"helper_functions.h"

float *read_matrix_file(std::string filename, int N, int M) {
    std::fstream matrix_file;
    matrix_file.open(filename.c_str(), std::ios::in);
    float *matrix = NULL;

    if (matrix_file.is_open()) {
        std::cout << "File " << filename << " is open \n";
        matrix = new float[M * N];
    } else {
        std::cout << "cant open " << filename << " exiting \n";
        return matrix;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            matrix_file >> matrix[i * M + j];
            //cout<< matrix[i*M+j]<<"\t";
        }
    }
    return matrix;
}


void print_matrix(float *matrix, int N, int M) {

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            std::cout << matrix[i * M + j] << " \t " << std::setprecision(5);
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------" << std::endl;
}


std::string *sequences_reader(std::string filename, int number_of_sequences) {

    int i = 0;
    std::string line;
    std::ifstream f(filename.c_str());
    std::string *sequence_set = NULL;
    if (f.is_open()) {
        sequence_set = new std::string[number_of_sequences];
        while (getline(f, line)) {
            sequence_set[i] = line;
            //cout<< sequence_set[i];
            i++;
        }
    } else { std::cout << "cannot open " << filename << std::endl; };

    f.close();
    return sequence_set;
}

std::string array_to_string(int *arr, int n_elem) {
    std::string reconstructed_states_string,t_string;
    std::stringstream ss;
////
    for (int t = 0; t <= n_elem; t++) {
        ss << arr[t];
        t_string = ss.str();
        ss.str(std::string());
        reconstructed_states_string.append(t_string);
        //cout<<to_string(arr[t]);
    }

    return reconstructed_states_string;
}

int *string_to_arr(std::string s_sequence) {

    int n = s_sequence.length();
    int *arr = new int[n];

    for (int i = 0; i < n; i++) {
        arr[i] = s_sequence[i] - '0';
    }
    return arr;
}

std::string annotate_path(std::string reconstructed_states) {

    std::string annoteted_states;
    int i, n = reconstructed_states.length();
    for (i = 0; i < n; i++) {
        if (reconstructed_states[i] == '1') { annoteted_states += 'F'; }
        if (reconstructed_states[i] == '2') { annoteted_states += 'L'; }
    }
    return annoteted_states;
}