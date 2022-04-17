/*
 * PRIMEX 1.1 Last update Dec 12, 2003
 * Read the README.txt file associated with this code for instructions
 *
 * Authors: Matej Lexa, Giorgio Valle, 2002.
 *
 * LICENSE: You are free to use this program for academic use. This 
 * includes development of new programs for academic use based on the code.
 * No guarantees are given as to the performance of the program.
 * Commercial use allowed with authors' permission only.
 */

#ifndef I_MEX_H
#define I_MEX_H

#define MAX_MATCHES 10000
#define MISMATCH 1
#define NUM_CL_OPTIONS 22 
#define WORD_LENGTH 10
#define MAX_CHAR 160
#define MAX_OLIGO 45
#define MAX_HITS 5000
#define MAX_LTABLE 1048576
#define MAX_SEQUENCE_LENGTH 2700000
#define MAX_CLONES 10
#define MAX_RESULTS 5000
#define PACKET_LENGTH 1024
#define MESSAGE_LENGTH 131077
#define PORT 30000

//Includes for the whole program
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// #include <iostream>
// #include <fstream>
// #include <math.h>
// #include <unistd.h>
// #include <ctype.h>
//Includes for the server part
// #include <sys/types.h>
// #include <sys/socket.h>
// #include <netinet/in.h>
// #include <netdb.h>
// #include <sys/wait.h>
// #include <signal.h>
#include "input.h"
// extern char filein[128], queryin[128], xfilein[128], xfileout[128], mfileout[128], primer_label[128];
// extern int word_length, ltable_step, max_hits, mismatches, mismatches2, ltable_only, primer_num, label_start, label_length, boundary_off, indels, server_port, genome_size;

#endif // I_MEX_H
