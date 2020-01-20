//
// Created by esemenc on 1/6/20.
//

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#define to_string( x ) static_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()

double*read_matrix_file( std::string filename,  int N,  int M);

void print_matrix(double*matrix,  int N,  int M);

std::string *sequences_reader(std::string filename,  int number_of_sequences);


