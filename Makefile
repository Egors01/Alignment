CC =g++

CFLAGS= -Wall -pedantic -ansi
all:
	$(CC) $(CFLAGS) main.cpp Pair_HMM.cpp helper_functions.cpp  -o run.o