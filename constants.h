//
// Created by esemenc on 12/9/19.
//

#ifndef UNTITLED_CONSTANTS_H
#define UNTITLED_CONSTANTS_H
static const int N_OBSERAVBLES = 6;
static const int N_STATES = 4;

//static const std::string PATH_EMISSIONS_FILE = "../Alignment_emissions.tsv";
//static const std::string PATH_TRANSITIONS_FILE = "../Alignment_transitions.tsv";

static const std::string PATH_EMISSIONS_FILE = "../Random_emissions.tsv";
static const std::string PATH_TRANSITIONS_FILE = "../Random_transitions.tsv";


static const std::string PATH_INPUT_FILE = "../in_pairs.txt";



static const std::string PATH_OUTPUT_FILE = "./predicted_state_path.txt";
static const std::string PATH_OUTPUT_ANNOTATED_FILE = "./predicted_annotation.txt";
static const std::string PATH_KNOWN_ANNOTATION_FILE = "./known_annotation.txt";
static const std::string PATH_OUTPUT_PERFORMANCE_FILE = "./performance.txt";
#endif //UNTITLED_CONSTANTS_H
