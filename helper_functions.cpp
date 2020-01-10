//
// Created by esemenc on 1/6/20.
//

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include"helper_functions.h"

double *read_matrix_file(std::string filename, int N, int M) {
    std::fstream matrix_file;
    matrix_file.open(filename.c_str(), std::ios::in);
    double *matrix = NULL;

    if (matrix_file.is_open()) {
        std::cout << "File " << filename << " is open \n";
        matrix = new double [M * N];
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


void print_matrix(double *matrix, int N, int M) {

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
        sequence_set = new std::string[number_of_sequences+1];
        while (i<number_of_sequences) {
            getline(f, line);
            sequence_set[i] = line;
            //cout<< sequence_set[i];
            std::cout << "read " << sequence_set[i] <<" " <<i<<std::endl;
            i++;
        }
    } else { std::cout << "cannot open " << filename << std::endl; };

    f.close();
    return sequence_set;
}

