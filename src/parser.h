/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_PARSER_H
#define SWIFT_PARSER_H

#define MAX_LINE_SIZE 128

struct model_params {
    //char *names[] = {"no_of_threads", "no_of_time_steps", "max_h", "ic_file"};

    int no_of_threads;
    int no_of_time_steps;
    float max_h;
    char *ic_file;
};

void parseFile(struct model_params *params, const char * file_name) {
 
    FILE *fp;
    char line[MAX_LINE_SIZE];

    /* Open file for reading */
    fp = fopen(file_name, "r");
    
    if(fp == NULL) {
      error("Error opening parameter file: %s",file_name);
    }
  
    /* Read until the end of the file is reached.*/
    while(!feof(fp)) {
      
      /* Read a line of the file */
      if(fgets(line,MAX_LINE_SIZE,fp)!=NULL) 
      {
        /* Check if the line contains a value */
        if(strchr(line,':')) {
         
          //int size = atoi(line);
          //int i;
          //char *p;
          //char *elements[size];
          //for(i=0; i<8; i++) elements[i] = (char *)malloc(MAX_LINE_SIZE);
          //

          ////for(i=0; i<size; i++) {
          //  fgets(line,MAX_LINE_SIZE,fp);
          //  strtok(line,":");
          //  strcpy(elements[0],line);

          //  p = (char *)strtok(NULL,":");
          //  strcpy(elements[1],p);
          ////}

          int i;
          char * token;
          char *elements[8];
          for(i=0; i<8; i++) elements[i] = (char *)malloc(MAX_LINE_SIZE);
          int element_count = 0;


          printf ("Splitting string: '%s' into tokens:\n",line);
          
          token = (char *)strtok(line,":");
         
          printf("Token:%s\n",token); 
          
          while (token != NULL) {
            strcpy(elements[element_count],token);  
            printf ("Element:%s\n",elements[element_count]);
            token = (char *)strtok (NULL, ":");
            element_count++;
          }
          
         //if(strcmp("no_of_threads",elements[0])) {
         //   params->no_of_threads = atoi(elements[1]);
         //} 
         //else if(strcmp("no_of_time_steps",elements[0])) {
         //   params->no_of_time_steps = atoi(elements[1]);
         //}
         //else if(strcmp("max_h",elements[0])) {
         //   params->max_h = atof(elements[1]);
         //}
         //else if(strcmp("ic_file",elements[0])) {
         //   params->ic_file = elements[1];
         //}
 
         //for(i=0; i<element_count; i++)
         //   printf("Element[%d]: %s\n",i,elements[i]);  
          
         //float fp_number;
          //char * str;
          //fscanf(fp, "%f %s", &fp_number,str);
          //
          ///* writing content to stdout */
          //printf("%s",line);
          ////printf("Float: %f, comment: %s\n",fp_number,str);
        }
      }
    
    }

    printf("no_of_threads:%d, no_of_time_steps:%d, max_h:%f, ic_file:%s\n",params->no_of_threads,params->no_of_time_steps,params->max_h,params->ic_file);

    fclose(fp);
}

void storeParam() {

}

#endif /* SWIFT_VECTOR_H */
