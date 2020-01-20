//
// Created by esemenc on 12/9/19.
//

#ifndef UNTITLED_CONSTANTS_H
#define UNTITLED_CONSTANTS_H
static const int N_OBSERAVBLES = 25;
static const int N_STATES = 5;

static const std::string PATH_PREFIX = "./";
static const std::string PATH_EMISSIONS_FILE_ALIIGN =  PATH_PREFIX+"Alignment_emissions.tsv";
static const std::string PATH_TRANSITIONS_FILE_ALIGN = PATH_PREFIX+"Alignment_transitions.tsv";

static const std::string PATH_EMISSIONS_FILE_RANDOM =   PATH_PREFIX+"Random_emissions.tsv";
static const std::string PATH_TRANSITIONS_FILE_RANDOM = PATH_PREFIX+"Random_transitions.tsv";

static const std::string PATH_INPUT_FILE = PATH_PREFIX+"in_pairs.txt";


static const std::string ALIGNMENT_OUTPUT_FILE_FULL = PATH_PREFIX+"alignments_with_states_and_ratios.txt";

static const std::string ALIGNMENT_OUTPUT_FILE =      PATH_PREFIX+"predicted_alignments.txt";
static const std::string ALIGNMENT_OUTPUT_FILE_PATH = PATH_PREFIX+"paths_alignment_model.txt";
static const std::string PROBAB_ALIGNMENT_FILE =      PATH_PREFIX+"probs_alignment_model.txt";
static const std::string PROBAB_RANDOM_FILE =         PATH_PREFIX+"probs_random_model.txt";
static const std::string LOG_RATIOS_FILE = PATH_PREFIX+"log_odds_ratios.txt";
#endif //UNTITLED_CONSTANTS_H
