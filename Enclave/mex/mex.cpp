#include "mex.h"
#include "ltable.h"

#define MAX_RAND_ACCESS 100


int main(int argc, char *argv[]){


  char *data = argv[1];
  char result[MESSAGE_LENGTH];  
  size_t size = 10;
  ecall_matchOligo(data,result,10);


  return 0;
}

