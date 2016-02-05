#include "parser.h"

int main(int argc, char *argv[]) {

  const char * input_file = argv[1];
  struct swift_params param_file;
  int no_of_threads = 0;
  int no_of_time_steps = 0;
  float max_h = 0.0f;
  char ic_file [MAX_LINE_SIZE];

  parseFile(input_file,&param_file);

  printParameters(&param_file);
  
  getParamInt(&param_file,"no_of_threads",&no_of_threads);
  getParamInt(&param_file,"no_of_time_steps",&no_of_time_steps);
  getParamFloat(&param_file,"max_h",&max_h);
  getParamString(&param_file,"ic_file",ic_file);
  
  printf("no_of_threads: %d, no_of_time_steps: %d, max_h: %f, ic_file: %s\n",no_of_threads, no_of_time_steps, max_h, ic_file);

}
