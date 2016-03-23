#include "parser.h"

int main(int argc, char *argv[]) {

  const char * input_file = argv[1];
  
  /* Create a structure to read file into. */ 
  struct swift_params param_file;
  
  /* Create variables that will be set from the parameter file. */
  int no_of_threads = 0;
  int no_of_time_steps = 0;
  float max_h = 0.0f;
  char ic_file [MAX_LINE_SIZE];

  /* Read the parameter file. */
  parser_read_file(input_file,&param_file);

  /* Print the contents of the structure. */
  parser_print_params(&param_file);
  
  /* Retrieve parameters and store them in variables defined above. 
   * Have to specify the name of the parameter as it appears in the 
   * input file: testParserInput.yaml.*/
  parser_get_param_int(&param_file,"no_of_threads",&no_of_threads);
  parser_get_param_int(&param_file,"no_of_time_steps",&no_of_time_steps);
  parser_get_param_float(&param_file,"max_h",&max_h);
  parser_get_param_string(&param_file,"ic_file",ic_file);
  
  /* Print the variables to check their values are correct. */
  printf("no_of_threads: %d, no_of_time_steps: %d, max_h: %f, ic_file: %s\n",no_of_threads, no_of_time_steps, max_h, ic_file);
}
