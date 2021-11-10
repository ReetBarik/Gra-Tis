#include "input_output.h"
void writeGraphGorderEdgeList(graph* G, char *filename) {
    //Get the iterators for the graph:
    long NVer     = G->numVertices;
    long NEdge    = G->numEdges;       //Returns the correct number of edges (not twice)
    long *verPtr  = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
    edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
    printf("NVer= %ld --  NE=%ld\n", NVer, NEdge);
    
    printf("Writing graph in Gorder format - Zero indexed and each edge represented ONLY ONCE!\n");
    printf("Graph will be stored in file: %s\n", filename);
    
    FILE *fout;
    fout = fopen(filename, "w");
    if (!fout) {
        printf("Could not open the file \n");
        exit(1);
    }
    printf("Num vertices: %ld    Num edges: %ld\n", NVer, NEdge);
    
    //Write the edges (lower triangle only):
    for (long v=0; v<NVer; v++) {
        long adj1 = verPtr[v];
        long adj2 = verPtr[v+1];
        for(long k = adj1; k < adj2; k++ ) {
            if (v < verInd[k].tail ) { //Print only once
                fprintf(fout, "%ld %ld\n", v, (verInd[k].tail));
            }
        }
    }
    fclose(fout);
    printf("Graph has been stored in file: %s\n",filename);
}//End of writeGraphGorderEdgeList()

void writeGraphGorderEdgeListReordered(graph* G, char *filename, long* old2New) {
    //Get the iterators for the graph:
    long NVer     = G->numVertices;
    long NEdge    = G->numEdges;       //Returns the correct number of edges (not twice)
    long *verPtr  = G->edgeListPtrs;   //Vertex Pointer: pointers to endV
    edge *verInd = G->edgeList;       //Vertex Index: destination id of an edge (src -> dest)
    printf("NVer= %ld --  NE=%ld\n", NVer, NEdge);
    
    printf("Writing graph in Gorder format - Zero indexed and each edge represented ONLY ONCE!\n");
    printf("Graph will be stored in file: %s\n", filename);
    
    FILE *fout;
    fout = fopen(filename, "w");
    if (!fout) {
        printf("Could not open the file \n");
        exit(1);
    }
    printf("Num vertices: %ld    Num edges: %ld\n", NVer, NEdge);
    //Write the edges (lower triangle only):
    for (long v=0; v<NVer; v++) {
        long adj1 = verPtr[v];
        long adj2 = verPtr[v+1];
        for(long k = adj1; k < adj2; k++ ) {
            if (v < verInd[k].tail) { //Print only once
                fprintf(fout, "%ld %ld\n", old2New[v], old2New[verInd[k].tail]);
            }
        }
    }
    fclose(fout);
    printf("Graph has been stored in file: %s\n",filename);
}//End of writeGraphGorderEdgeList()
