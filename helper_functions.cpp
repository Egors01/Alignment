//
// Created by esemenc on 1/6/20.
//

#include <string>
#include <iostream>
#include <iomanip>
#include"helper_functions.h"

// #define to_string( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()

double *create_2dim_matrix(const int N_lines, const int M_columns){

    // allocates memory for matrix and adds the size N,M to the last two elements  N lines M columns
    // access matrix [M * i + j] where i = {0..N-1}, j = {0..M-1}
    // size stored N := matrix[M*N]; M := matrix[N*M + 1]

    double * m_pointer;
    m_pointer = new double [N_lines*M_columns+2];
    m_pointer[N_lines*M_columns] = N_lines;
    m_pointer[N_lines*M_columns+1] = M_columns;
    return m_pointer;
}


double *read_matrix_file(std::string filename, const int N, const int M) {

    // creates *double 2 dim matrix and fills it with values from file
    // create_2dim_matrix is used to allocate memory

    std::fstream matrix_file;
    std::string line;
    int line_count = 0;
    matrix_file.open(filename.c_str(), std::ios::in);
    if (matrix_file.is_open()){
        while (getline(matrix_file, line))
            line_count++;
        matrix_file.close();

        if (line_count!=N){
            std::string error_msg =  "Read matrix from file(): wrong line count provided. Matrix will be misread.\n"
                                     "Recieved number of lines: "+to_string(N)+" Actual number of lines: "+to_string(line_count);
            throw std::invalid_argument(error_msg);
        }
    }
    else {
        std::cout << "cannot open " << filename << " exiting... \n";
        return NULL;
    }

    matrix_file.open(filename.c_str(), std::ios::in);
    if (matrix_file.is_open()) {
        std::cout << "File " << filename << " is open \n";
        double* matrix = create_2dim_matrix(N,M);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                matrix_file >> matrix[i * M + j];
                //cout<< matrix[i*M+j]<<"\t";
            }
        }
        return matrix;
    }
    else{
        std::cout << "failed to read from file " << filename << " exiting... \n";
        return NULL;
    }
}


void print_matrix(double *matrix, const int N, const int M) {

    //  Expect the last two elements to store the size of matrix
    //  double *read_matrix_file(std::string filename, const int N, const int M)

    if (matrix[N*M] !=N or matrix[N*M +1] !=M ){
        std::cout<<"Warning in print matrix(): Size stored in last two additional elements :"<<matrix[N*M]<<", "<<matrix[N*M+1]<<std::endl;
        throw std::invalid_argument("Size stored in last two elements does not match with the M,N provided: "+
        to_string(N)+", "+to_string(M) +"\nAdd the two extra elements to allocated memory. Store the actual N,M size there ");
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
            i++;
        }
    } else { std::cout << "cannot open " << filename << std::endl; };

    f.close();
    return sequence_set;
}

