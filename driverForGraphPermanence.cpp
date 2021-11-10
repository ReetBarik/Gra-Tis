//
// Created by sriramsrinivas@unomaha.edu on 6/1/20.
//
#include "defs.h"
#include "input_output.h"
#include "basic_comm.h"
#include "basic_util.h"
#include "utilityClusteringFunctions.h"
#include "color_comm.h"
#include "sync_comm.h"
#include "permanence.h"

using namespace std;

int main(int argc, char** argv)
{
    
    printf("I am here");
    clustering_parameters opts;
    if (!opts.parse(argc, argv)) {
        return -1;
    }
    int nT = 1; //Default is one thread
#pragma omp parallel
    {
        nT = omp_get_num_threads();
    }
    if (nT < 1) {
        printf("The number of threads should be greater than one.\n");
        return 0;
    }
    
    double time1, time2;
    graph* G = (graph *) malloc (sizeof(graph));
    
    
    int fType = opts.ftype; //File type
    char *inFile = (char*) opts.inFile;
    switch (fType) {
        case 1: parse_MatrixMarket_Sym_AsGraph(G, inFile); break;
        case 2: parse_Dimacs9FormatDirectedNewD(G, inFile); break;
        case 3: parse_PajekFormat(G, inFile); break;
        case 4: parse_PajekFormatUndirected(G, inFile); break;
        case 5: loadMetisFileFormat(G, inFile); break;
        case 6: parse_UndirectedEdgeList(G, inFile); break;
            //parse_UndirectedEdgeListDarpaHive(G, inFile); break;
        case 7: /*parse_DirectedEdgeList(G, inFile);*/ printf("This routine is under development.\n"); exit(1); break;
        case 8: parse_SNAP(G, inFile); break;
        case 9: parse_EdgeListBinaryNew(G,inFile); break;
        case 10:
#ifdef USEHDF5
            //parse_EdgeListCompressedHDF5(G,inFile);
            parse_EdgeListCompressedHDF5NoDuplicates(G,inFile);
#endif
            break;
        case 11: parse_UndirectedEdgeListFromJason(G, inFile); break;
        case 12: parse_UndirectedEdgeListWeighted(G, inFile); break; // for John F's graphs
        case 13: parse_UndirectedEdgeListDarpaHive(G, inFile); break;
        case 14: parse_EdgeListFromGorder(G, inFile); break;
        default:  cout<<"A valid file type has not been specified"<<endl; exit(1);
    }
    
    displayGraphCharacteristics(G);
    
    int threadsOpt = 0;
    if(opts.threadsOpt)
        threadsOpt =1;
    threadsOpt =1;
    
    int replaceMap = 0;
    if(  opts.basicOpt == 1 )
        replaceMap = 1;
    if( opts.VF ) {
        printf("Vertex following is enabled.\n");
        time1 = omp_get_wtime();
        long numVtxToFix = 0; //Default zero
        long *C = (long *) malloc (G->numVertices * sizeof(long)); assert(C != 0);
        numVtxToFix = vertexFollowing(G,C); //Find vertices that follow other vertices
        if( numVtxToFix > 0) {  //Need to fix things: build a new graph
            printf("Graph will be modified -- %ld vertices need to be fixed.\n", numVtxToFix);
            graph *Gnew = (graph *) malloc (sizeof(graph));
            long numClusters = renumberClustersContiguously(C, G->numVertices);
            buildNewGraphVF(G, Gnew, C, numClusters);
            //Get rid of the old graph and store the new graph
            free(G->edgeListPtrs);
            free(G->edgeList);
            free(G);
            G = Gnew;
        }
        free(C); //Free up memory
        printf("Graph after modifications:\n");
        displayGraphCharacteristics(G);
        
    }//End of if( VF == 1 )
    
    
    
    runPermanence(G, nT);
    
    
    
    return 0;
}
