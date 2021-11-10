

#ifndef GRAPPOLO_PERMANENCE_H
#define GRAPPOLO_PERMANENCE_H

#include <algorithm>
#include <iostream>
#include <list>

typedef pair<int,double>int_double;
struct Perm_Info
{
    double degree; //vertex degree
    vector<list<int>> comm_neighs;//neighbors in each community
    vector<int>  comms;//list of neighboring communities including own
    vector<double> in_degree; //weighed degree for each community
    vector<double> cc; //clustering coefficient for each community
    std::map<int,int> comm_map; //Map between community id and their position in comm_neighs

    int Comm;  //community id
    double perm ; //permanence

    //Constructor
    Perm_Info()
    {
        degree=0.0;
        comm_neighs.clear();
        in_degree.clear();
        comms.clear();
        comm_map.clear();
        cc.clear();

        Comm=-1;
        perm=-1.0;
    }
};

void get_neighbors(int node,long  *vtxPtr ,  edge  *vtxInd, vector<int> *neighbors);

void compute_CC(long *NV,long  *vtxPtr , edge  *vtxInd , vector<int> node_set, double *cc);

//Seeding Algorithm
void degreeMin_seed(long *NV,long  *vtxPtr ,  edge  *vtxInd , vector <Perm_Info> *vector_info);

//Permanance calculations:
void initialize_perminfo(long *NV,long  *vtxPtr ,  edge  *vtxInd , int *max_comm, vector<Perm_Info> *vector_info);

double get_permanence_old(int i, Perm_Info myvector, int mycomm);

void cluster_by_permanence_old(long *NV,long  *vtxPtr ,  edge  *vtxInd , int max_comm,vector<Perm_Info> *vector_info, bool allow_singleton);

void printEdges(long NV, long  *vtxPtr,  edge * vtxInd);

void runPermanence(graph *G, int numThreads);

#endif  //GRAPPOLO_PERMANENCE_H
