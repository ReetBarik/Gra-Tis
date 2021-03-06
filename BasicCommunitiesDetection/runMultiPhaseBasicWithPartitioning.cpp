// ***********************************************************************
//
//            Grappolo: A C++ library for graph clustering
//               Mahantesh Halappanavar (hala@pnnl.gov)
//               Pacific Northwest National Laboratory
//
// ***********************************************************************
//
//       Copyright (2014) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************

#include "defs.h"
#include "basic_comm.h"
#include "basic_util.h"
#include "utilityGraphPartitioner.h"

using namespace std;

//WARNING: This will overwrite the original graph data structure to
//         minimize memory footprint
// int nParts: Number of partitions to call Metis
// Return: C_orig will hold the cluster ids for vertices in the original graph
//         Assume C_orig is initialized appropriately
//WARNING: Graph G will be destroyed at the end of this routine
void runMultiPhaseBasicWithPartitioning(graph *G, long *C_orig, int basicOpt, long minGraphSize, long minPhaseSize,
                        double threshold, double C_threshold, int numThreads, int threadsOpt, int nParts)
{
    double totTimeClustering=0, totTimeBuildingPhase=0, totTimeColoring=0, tmpTime=0;
    int tmpItr=0, totItr = 0;
    long NV = G->numVertices;
    long NVOrig = NV;
    
    /* Step 1: Find communities */
    double prevMod = -1;
    double currMod = -1;
    long phase = 1;
    
    graph *Gnew; //To build new hierarchical graphs
    long numClusters;
    long *C = (long *) malloc (NV * sizeof(long));
    assert(C != 0);
#pragma omp parallel for
    for (long i=0; i<NV; i++) {
        C[i] = -1;
    }
    
    while(1){
        printf("===============================\n");
        printf("Phase %ld\n", phase);
        printf("===============================\n");
        prevMod = currMod;
        
        
        if(basicOpt == 1){
            currMod = parallelLouvianMethodNoMap(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        }else if(threadsOpt == 1){
            currMod = parallelLouvianMethod(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
            //currMod = parallelLouvianMethodApprox(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        }else{
            currMod = parallelLouvianMethodScale(G, C, numThreads, currMod, threshold, &tmpTime, &tmpItr);
        }
        
        totTimeClustering += tmpTime;
        totItr += tmpItr;
        
        //Renumber the clusters contiguiously
        numClusters = renumberClustersContiguously(C, G->numVertices);
        printf("Number of unique clusters: %ld\n", numClusters);
        
        //printf("About to update C_orig\n");
        //Keep track of clusters in C_orig
        if(phase == 1) {
#pragma omp parallel for
            for (long i=0; i<NV; i++) {
                C_orig[i] = C[i]; //After the first phase
            }
        } else {
#pragma omp parallel for
            for (long i=0; i<NV; i++) {
                assert(C_orig[i] < G->numVertices);
                if (C_orig[i] >=0)
                    C_orig[i] = C[C_orig[i]]; //Each cluster in a previous phase becomes a vertex
            }
        }
        printf("Done updating C_orig\n");

        if (numClusters <= nParts) {
            // shouldn't send a graph to METIS with n < #partitions
            break;
        }
        
        //build the graph for next phase
        //In case coloring is used, make sure the non-coloring routine is run at least once

        Gnew = (graph *) malloc (sizeof(graph)); assert(Gnew != 0);
        tmpTime =  buildNextLevelGraphOpt(G, Gnew, C, numClusters, numThreads);
        totTimeBuildingPhase += tmpTime;
        //Free up the previous graph
        //free(G->edgeListPtrs);
        //free(G->edgeList);
        //free(G);
        G = Gnew; //Swap the pointers
        G->edgeListPtrs = Gnew->edgeListPtrs;
        G->edgeList = Gnew->edgeList;
        
        //Free up the previous cluster & create new one of a different size
        //free(C);
        C = (long *) malloc (numClusters * sizeof(long)); assert(C != 0);
            
#pragma omp parallel for
        for (long i=0; i<numClusters; i++) {
            C[i] = -1;
        }
        phase++; //Increment phase number

        //Break if too many phases or iterations
        if((phase > 200)||(totItr > 100000)) {
            // printf("too many phases or iterations\n");
            break;
        }
        //Break if the graph will be less than a given parameter
        if(numClusters < minPhaseSize) {
            if(phase == 1)
                prevMod = currMod;
            // printf("graph will be less than a given parameter\n");
            break;
        }
        //Break if minGraphsize is set to -1
        if(minPhaseSize == -1) {
            if(phase == 1)
                prevMod = currMod;
            // printf("minGraphsize is set to -1\n");
            break;
        }

        //Check for modularity gain
        if( (currMod - prevMod) < threshold ) {
            break; //Modularity gain is not enough. Exit.
        }  //End of  if( (currMod - prevMod) > threshold )
    } //End of while(1)


    /* Step 2: Find Partitions */
    //free(C);
    long numV = G->numVertices;
    C = (long *) malloc (numV * sizeof(long));
    assert(C != 0);
    #pragma omp parallel for
    for (long i=0; i<numV; i++) {
        C[i] = -1;
    }
    //Call the K-way graph partitioner:
    double time1 = omp_get_wtime();
    MetisGraphPartitioner(G, C, nParts, C_orig, NVOrig);
    double time2 = omp_get_wtime();
    //Renumber the clusters:
    #pragma omp parallel for
    for (long i=0; i<NV; i++) {
        assert(C_orig[i] < G->numVertices);
        if (C_orig[i] >=0)
            C_orig[i] = C[C_orig[i]]; //Each cluster in a previous phase becomes a vertex
    }
    
    printf("********************************************\n");
    printf("*********    Compact Summary   *************\n");
    printf("********************************************\n");
    printf("Number of threads              : %ld\n", numThreads);
    printf("Total number of phases         : %ld\n", phase);
    printf("Total number of iterations     : %ld\n", totItr);
    printf("Final number of clusters       : %ld\n", numClusters);
    printf("Final modularity               : %lf\n", prevMod);
    printf("Total time for clustering      : %lf\n", totTimeClustering);
    printf("Total time for building phases : %lf\n", totTimeBuildingPhase);
    printf("********************************************\n");
    printf("Number of partitions (Metis)   : %ld\n", nParts);
    printf("Time for graph partitioning    : %lf\n", time2 - time1);
    printf("********************************************\n");
    printf("TOTAL TIME                     : %lf\n", (totTimeClustering+totTimeBuildingPhase+totTimeColoring+(time2 - time1)) );
    printf("********************************************\n");
    
    //Clean up:
    //free(C);
    if(G != 0) {
    //    free(G->edgeListPtrs);
    //    free(G->edgeList);
    //    free(G);
    }
}//End of runMultiPhaseBasicWithPartitioning()
