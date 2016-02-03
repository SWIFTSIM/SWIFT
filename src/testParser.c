#include <stdio.h>

#include "parser.h"

int main(int argc, char *argv[]) {

  struct model_params params;
  const char * input_file = argv[1];

  //parseFile(&params,"testInput.dat");
  parseFile(&params,input_file);

}
