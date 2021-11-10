//
// Created by sriramsrinivas@unomaha.edu on 5/30/20.
//

#include "defs.h"
#include "utilityClusteringFunctions.h"
# include "permanence.h"

using namespace std;

void printEdges(long NV, long  *vtxPtr,  edge * vtxInd)
{

    for(int i=0;i<NV;i++) {
        long adj1 = vtxPtr[i];	    //Begin
        long adj2 = vtxPtr[i+1];	//End
        double totalWt = 0;
        for(long j=adj1; j<adj2; j++) {
            cout <<vtxInd[j].tail<<" ";
        }
        cout<<"\n";
    }

    ////    for(int i=0;i<vector_info.size();i++)
////    {
////
////        cout<<i<<"\t";
////        cout<<vector_info.at(i).degree<<"\t";
////
////        cout<<vector_info.at(i).Comm<<"\t";
////
////
////        cout<<vector_info.at(i).perm<<"\t";
////
////
////        for(int j=0;j<vector_info.at(i).comms.size();j++)
////        {
////            cout<<vector_info.at(i).comms.at(j)<<"\t";
////        }
//
////        for(int j=0;j<vector_info.at(i).in_degree.size();j++)
////        {
////            cout<<vector_info.at(i).in_degree.at(j)<<"\t";
////        }
////
////
////        for(int j=0;j<vector_info.at(i).cc.size();j++)
////        {
////            cout<<vector_info.at(i).cc.at(j)<<"\t";
////        }
//
//
//        cout<<"\n";
//
//
//    }


}



void runPermanence(graph *G, int numThreads)
{

    cout<<"Inside Permanence Code"<<"\n";
    cout<<"Number of Threads"<<numThreads<<"\n";
    long    NV        = G->numVertices;
    long    NS        = G->sVertices;
    long    NE        = G->numEdges;
    long    *vtxPtr   = G->edgeListPtrs;
    edge    *vtxInd   = G->edgeList;

    //printEdges(NV,vtxPtr,vtxInd);
    long* pastCommAss = (long *) malloc (NV * sizeof(long)); assert(pastCommAss != 0);
    long* currCommAss = (long *) malloc (NV * sizeof(long)); assert(currCommAss != 0);

    /* Mahantesh's implementation of initial Community
     *  We will be using Seeding initially past and current community assigniment will be based on degree seeding even for hybrid implementation
     * */
    initCommAss(pastCommAss, currCommAss, NV);
    vector<Perm_Info> vector_info;
    int max_comms;
    Perm_Info dummy_perm;
    vector_info.resize(NV, dummy_perm);
    degreeMin_seed(&NV,vtxPtr,vtxInd,&vector_info);
    cout <<"Seeding Done \n";
    initialize_perminfo(&NV,vtxPtr,vtxInd,&max_comms,&vector_info);
    cout <<"Initialization Done \n";

    //Compute Permanence
    //Set this to false if you do not want to have singleton vertices
    bool allow_singleton=true;
    cluster_by_permanence_old(&NV,vtxPtr,vtxInd, max_comms, &vector_info,allow_singleton);
    double permanence=0.0;
    for(int z=0;z< vector_info.size();z++)
    {   permanence=permanence+vector_info[z].perm;
        cout << z << "::"<< vector_info[z].Comm << "::" << vector_info[z].perm <<"\n";}
    cout << "Total Permanence : ==" << permanence  << "\n";

}//End of runPermanence()

