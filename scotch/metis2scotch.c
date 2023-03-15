#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MAX_LINE_LENGTH 1024
#include <scotch.h>


/* function declaration */
void read_ncells(const char* filename, int *ncells);
void read_metis_edges(const char* filename, int *adjncy);
void read_metis_weights(const char* filename, int *weights_e, int *weights_v);

int main () {
    // Read in the metis simple file dump
    char* edgesname = "";
    // Read in the metis weights file dump
    char* weightsname = "";
    /* local variable definition */
    int ncells = 0;
    read_ncells(edgesname, &ncells);
    int *xadj;
    if ((xadj = (int *)malloc(sizeof(int) * (ncells + 1))) == NULL)
        printf("Failed to allocate xadj buffer.");
    int *adjncy;
    if ((adjncy = (int *)malloc(sizeof(int) * 26 * ncells)) == NULL)
        printf("Failed to allocate adjncy array.");
    int *weights_v = NULL;
    if ((weights_v = (int *)malloc(sizeof(int) * ncells)) == NULL)
        printf("Failed to allocate vertex weights array");
    int *weights_e = NULL;
    if ((weights_e = (int *)malloc(26 * sizeof(int) * ncells)) == NULL)
        printf("Failed to allocate edge weights array");
    int *regionid;
    if ((regionid = (int *)malloc(sizeof(int) * ncells)) == NULL)
        printf("Failed to allocate regionid array");

    read_metis_edges(edgesname, adjncy);
    read_metis_weights(weightsname, weights_e, weights_v);
    // Setting up the Scotch graph
    SCOTCH_Graph graph;
    SCOTCH_Num baseval = 0;
    SCOTCH_Num vertnbr = ncells;
    SCOTCH_Num *verttab;   /* Vertex array [vertnbr+1] */
    SCOTCH_Num *vendtab = NULL;   /* Vertex array [vertnbr]   */
    SCOTCH_Num *velotab;   /* Vertex load array        */
    SCOTCH_Num *vlbltab = NULL;   /* Vertex label array       */ 
    SCOTCH_Num edgenbr = (26 * vertnbr);       /* Number of edges (arcs)   */    
    SCOTCH_Num *edgetab;   /* Edge array [edgenbr]     */
    SCOTCH_Num *edlotab;

    verttab = (SCOTCH_Num*) malloc((vertnbr+1) * sizeof(SCOTCH_Num));
    velotab = (SCOTCH_Num*) malloc((vertnbr) * sizeof(SCOTCH_Num));
    edgetab = (SCOTCH_Num*) malloc(edgenbr * sizeof(SCOTCH_Num));
    edlotab = (SCOTCH_Num*) malloc(edgenbr * sizeof(SCOTCH_Num));

    printf("Done the set up \n");
    int i;
    for (i = 0; i <= vertnbr; i++) {
        verttab[i] = i*26;
        velotab[i] = weights_v[i];
    }

    for (i = 0; i < edgenbr; i++) {
        edgetab[i] = adjncy[i];
        edlotab[i] = weights_e[i];
    }

    printf("Initialise graph \n");
    SCOTCH_graphInit(&graph);

    if (SCOTCH_graphBuild(&graph, baseval, vertnbr, verttab, vendtab, velotab, NULL, edgenbr, edgetab, edlotab) != 0) {
        printf("Error: Cannot build Scotch Graph.\n");
        exit(EXIT_FAILURE);
    }

    printf("Scotch Graph built successfully.\n");

    FILE *file = fopen("", "w");
    if (file == NULL) {
        printf("Error: Cannot open output file.\n");
        exit(EXIT_FAILURE);
    }

    if (SCOTCH_graphSave(&graph, file) != 0) {
        printf("Error: Cannot save Scotch Graph.\n");
        exit(EXIT_FAILURE);
    }

    printf("Scotch Graph saved to file.\n");

    fclose(file);
    SCOTCH_graphExit(&graph);

    // Free memory
    free(verttab);
    free(velotab);
    free(edgetab);
    free(edlotab);
    free(xadj);
    free(adjncy);
    free(weights_v);
    free(weights_e);
    free(regionid);
   return 0;
}


void read_ncells(const char* filename, int *ncells) {
    // Read in the number of cells/vertices
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("printf opening file %s\n", filename);
        return;
    }
    char line[MAX_LINE_LENGTH];
    int line_num = 0;
    if (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        char * pch;
        pch = strtok (line," ");
        *ncells = atoi(pch);
    }
    fclose(fp);
}



void read_metis_edges(const char* filename, int *adjncy) {
    // Read in the vertex neighbours
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("printf opening file %s\n", filename);
        return;
    }
    char line[MAX_LINE_LENGTH];
    int index = 0;
    int line_num = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        if (line_num > 0) {
            char * pch;
            pch = strtok(line," ,.-");
            while (pch != NULL){
                adjncy[index] = atoi(pch);
                pch = strtok (NULL, " ,.-");
                index +=1;
            }
        }
    line_num += 1;
    }
    fclose(fp);
}

void read_metis_weights(const char* filename, int *weights_e, int *weights_v) {
    // Read in the vertex and edge weights
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("printf opening file %s\n", filename);
        return;
    }

    char line[MAX_LINE_LENGTH];
    int v_index = 0;
    int e_index = 0;
    int line_num = 0;
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        int vertex_ind = 0;
        if (line_num > 0) {
            printf ("Line number is %i\n",line_num);
            char * pch;
            pch = strtok(line," ,.-");
            weights_v[v_index] = atoi(pch);
            printf ("Vertex weight is %s\n",pch);
            while (pch != NULL){
                if (vertex_ind>0){
                    printf ("Edge Weight is %s\n",pch);
                    weights_e[e_index] = atoi(pch);
                    e_index +=1;
                }
            vertex_ind += 1;
            pch = strtok (NULL, " ,.-");
            }
        v_index+=1;
        }
        line_num += 1;
    }
    fclose(fp);
}