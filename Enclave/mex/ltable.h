#ifndef I_LTABLE_H 
#define I_LTABLE_H 

#include "mex.h"

#if defined(__cplusplus)
extern "C" {
#endif 
// #include "sequence.h"

  // SequenceSet* genome;
  // char oligo[MAX_CHAR];
  // int word_length,ltable_step;
  
  // // int ltable_count[MAX_LTABLE];
  // // int *ltable[MAX_LTABLE];

  // int global_pos;
  // int c;
  // char match_sense[MAX_MATCHES];
  // int match_shift[MAX_MATCHES], matches[MAX_MATCHES];
  // int query_len, value, pos;
  // int query_full_len;
  // int j,sense_change,match_count;
  // // Parameters set from command lines
  // int mismatches, mismatches2, indels, max_hits, server_port, genome_size;
  // char primer_label[128];//,xfilein[128],xfileout[128];





  // int ecall_matchOligo(char* oligo,size_t size);
  int oligo_to_binary(char *string);  
  void reverse(char* string,char* new_string,size_t length);
  void find_seed_matches(char *string);
  // void justCall(char *string);
  
#if defined(__cplusplus)
}
#endif 


#endif // I_LTABLE_H
