/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

#include "parser.h"

/* Needs to be included so that strtok returns char * instead of a int *. */
#include <string.h>
#include <stdlib.h>
#include "error.h"

void parseFile(const char *file_name, struct swift_params *params) {
 
    FILE *fp;
    
    params->count = 0;

    /* Open file for reading */
    fp = fopen(file_name, "r");
    
    if(fp == NULL) {
      error("Error opening parameter file: %s",file_name);
    }
  
    /* Read until the end of the file is reached.*/
    while(!feof(fp)) {
      readParameter(fp,params); 
    }

    fclose(fp);
}

static int countChar(char *str, char val) {

    int i, count = 0;

    for(i=0; i<strlen(str); i++) {
      
      /* Check if the line contains the character */
      if(str[i] == val) count++; 

    }

    return count;
}

static void readParameter(FILE *fp, struct swift_params *params) {

    char line [MAX_LINE_SIZE];
    char trim_line [MAX_LINE_SIZE];

    /* Read a line of the file */
    if(fgets(line,MAX_LINE_SIZE,fp) != NULL) {
        
        char *token;
        /* Remove comments */
        token = strtok(line,COMMENT);
        strcpy(trim_line,token);
        
        /* Check if the line contains a value */
        if(strchr(trim_line,VALUE)) {
          /* Check for more than one parameter on the same line. */
          if(countChar(trim_line,VALUE) > 1) { 
            error("Found more than one parameter in '%s', only one allowed.",line);
          }
          else { 
            /*Take first token as the parameter name. */
            token = strtok(trim_line,VALUE);
            strcpy(params->data[params->count].name,token);
            
            /*Take second token as the parameter value. */
            token = strtok (NULL, " #\n");
            strcpy(params->data[params->count++].value,token);
          }
        }
    }
}

void getParamInt(struct swift_params *params, char * name, int * retParam) {
  
   int i; 
   
   for(i=0; i<MAX_NO_OF_PARAMS; i++) {

     /*strcmp returns 0 if both strings are the same.*/
     if(!strcmp(name,params->data[i].name)) {
       
       /* Check if value is 0, to avoid integers of zero causing an error. */
       if(!strcmp("0",params->data[i].value)) {
           *retParam = 0;
       }
       /* Check if the value is a string. */
       else if(!atoi(params->data[i].value)) {
         error("Tried to parse '%s', when expecting an integer.",params->data[i].value);       
       }
       /* Check if the value is a float.*/
       else if( strchr(params->data[i].value,'.') || strchr(params->data[i].value,'e')|| 
                strchr(params->data[i].value,'E') ) {
         error("Tried to parse '%s', when expecting an integer.",params->data[i].value);       
       }
       else {
         *retParam = atoi(params->data[i].value);       
       }
       return;
       
     }
   }
}

void getParamFloat(struct swift_params *params, char * name, float * retParam) {
  
   int i; 
   
   for(i=0; i<MAX_NO_OF_PARAMS; i++) {

     /*strcmp returns 0 if both strings are the same.*/
     if(!strcmp(name,params->data[i].name)) {
       *retParam = atof(params->data[i].value);       
       return;
     }
   }
}

void getParamString(struct swift_params *params, char * name, char * retParam) {
  
   int i; 
   
   for(i=0; i<MAX_NO_OF_PARAMS; i++) {

     /*strcmp returns 0 if both strings are the same.*/
     if(!strcmp(name,params->data[i].name)) {
       strcpy(retParam, params->data[i].value);       
       return;
     }
   }
}

void printParameters(struct swift_params *params) {

    int i;

    printf("\n--------------------\n");
    printf("SWIFT Parameter File\n");
    printf("--------------------\n");
    
    for(i=0; i<params->count; i++) {
        printf("Name: %s\n",params->data[i].name);
        printf("Value: %s\n",params->data[i].value);
    }

}
