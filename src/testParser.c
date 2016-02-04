#include <stdio.h>

#include "parser.h"

int main(int argc, char *argv[]) {

  const char * input_file = argv[1];
  struct swift_params param_file;
 
  parseFile(&param_file,input_file);

  printParameters(&param_file);
}
