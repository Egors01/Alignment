//
// Created by esemenc on 1/6/20.
//

#include <string>
#include <iostream>
#include <iomanip>
#include"helper_functions.h"
#include <sstream>

#define to_string( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()

double *read_matrix_file( std::string filename, const int N, const int M) {
    std::fstream matrix_file;
    std::string line;
    int line_count = 0;

    matrix_file.open(filename.c_str(), std::ios::in);

    while (getline(matrix_file, line))
        line_count++;
    matrix_file.close();

    if (line_count!=N){
        std::string error_msg =  "Read matrix from file(): wrong line count provided. Matrix will be misread.\nRecieved number of lines: "+to_string(N)+" Actual number of lines: "+to_string(line_count);
        throw std::invalid_argument(error_msg);
    }

    matrix_file.open(filename.c_str(), std::ios::in);
    double *matrix = NULL;

    if (matrix_file.is_open()) {
        std::cout << "File " << filename << " is open \n";
        matrix = new double [M * N];
    } else {
        std::cout << "cannot open " << filename << " exiting \n";
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


void print_matrix(double *matrix, const int N, const int M) {


    if (matrix[(N-1)*M + (M-1) +1]<=1e10 or matrix[(N-1)*M + (M-1) +1]>=1e-10){
        std::cout<<"Warning in print matrix(): Suspicios values in matrix. M,N provided:"<<M<<N<<std::endl;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            std::cout << matrix[i * M + j] << " \t " << std::setprecision(5);
        }
        std::cout << std::endl;
    }
    std::cout << "--------------------" << std::endl;
}


std::string *sequences_reader(const std::string filename, const int number_of_sequences) {

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
            //std::cout << "read " << sequence_set[i] <<" " <<i<<std::endl;
            i++;
        }
    } else { std::cout << "cannot open " << filename << std::endl; };

    f.close();
    return sequence_set;
}

