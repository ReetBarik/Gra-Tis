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

using namespace std;

clustering_parameters::clustering_parameters()
: ftype(7), strongScaling(false), output(false), VF(false), coloring(0), numColors(16), percentage(100), syncType(0),
threadsOpt(false), basicOpt(0), C_thresh(0.01), minGraphSize(100000), minPhaseSize(1000), threshold(0.000001), compute_metrics(false), vtxOrderingAlgorithm(0), nParts(16)
{}

void clustering_parameters::usage() {
    cout << "***************************************************************************************"<< endl;
    cout << "Basic usage: Driver <Options> FileName\n";
    cout << "***************************************************************************************"<< endl;
    cout << "Input Options: \n";
    cout << "***************************************************************************************"<< endl;
    cout << "File-type  : -f <1-8>   -- default=7" << endl;
    cout << "File-Type  : (1) Matrix-Market format     (2) DIMACS#9 format \n";
    cout << "           : (3) Pajek (each edge once)   (4) Pajek (twice)   (5) Metis (DIMACS#10) \n";
    cout << "           : (6) Simple edgelist (stored once; unweighted) (7) Simple edgelist (stored once; weighted) \n";
    cout << "           : (8) SNAP format (9) Binary format  (10) HDF5 (11) Data from Jason   \n";
    cout << "           : (12) Data from Jason (13) Graph Challenge format (HIVE) \n";
    cout << "           : (14) Gorder format (zero based edge list) \n";
    cout << "--------------------------------------------------------------------------------------" << endl;
    cout << "Threshold       : -t <value> [default=0.000001]" << endl;
    cout << "Percentage      : -p <value> [default=100]"      << endl;
    cout << "Min-size (Phase): -h <value> [default=1000] (-1) Run only one phase"   << endl;
    cout << "BasicOpt        : -b <value>  [default=0]  (0) basic (1) replaceMap    " << endl;
    cout << "syncType        : -y <value>  [default=0]  (1) FullSync (2) NeighborSync (3) EarlyTerm (4) 1+3   " << endl;
    cout << "--------------------------------------------------------------------------------------" << endl;
    cout << "Coloring        : -c <value>  [default=0]  (0) None (1) d1-coloring (2) balanced coloring " <<endl;
    cout << "                                           (3) incomplete coloring -- need to set n (#colors)" << endl;
    cout << "Min-size (color): -m <value> [default=100000] (-1) Run only one phase"   << endl;
    cout << "C-threshold     : -d <value> [default=0.01]"     << endl;
    cout << "# of colors     : -n <value> [default=16]"       << endl;
    cout << "--------------------------------------------------------------------------------------" << endl;
    cout << "Vtx Ordering    : -z <value> [default=0]  (0) Print Statistics (1) RCM (2) Grappolo    " << endl;
    cout << "                                          (3) Print graph for Gorder                   " << endl;
    cout << "                                          (4) RCM-Strict (5) Grappolo-ML (6) Degree    " << endl;
    cout << "--------------------------------------------------------------------------------------" << endl;
    cout << "Strong scaling  : -s   [default=false] " << endl;
    cout << "VF              : -v   [default=false] " << endl;
    cout << "Output          : -o   [default=false] " << endl;
    cout << "Metrics         : -r   [default=false] " << endl;
    cout << "--------------------------------------------------------------------------------------" << endl;
    cout << "Num partitions  : -k   [default=16   ] " << endl;
    cout << "***************************************************************************************"<< endl;
}//end of usage()

bool clustering_parameters::parse(int argc, char *argv[]) {
    static const char *opt_string = "c:z:n:p:b:y:svof:t:d:m:h:r:k:";
    int opt = getopt(argc, argv, opt_string);
    while (opt != -1) {
        switch (opt) {
            case 'c': coloring = atoi(optarg);
                if((coloring <0) || ((coloring >3))) {
                    cout << "Coloring: " << coloring << endl;
                    cout << "Coloring needs to be an integer between 0 and 2. Options:" << endl;
                    cout << "0) None (1) d1-coloring (2) balanced coloring (3) incomplete coloring -- need to set n (#colors)" << endl;
                    coloring = 0;
                }
                break;
            case 'y': syncType = atol(optarg); break;
            case 'b': basicOpt = atol(optarg); break;
            case 's': strongScaling = true; break;
            case 'v': VF = true; break;
            case 'o': output = true; break;
            case 'r': compute_metrics = true; break;
            case 'p': percentage = atol(optarg);
                if((percentage <0) || ((percentage >100))) {
                    cout << "Percentage is set to default of 100 percent" << endl;
                    percentage = 100;
                }
                break;

            case 'f': ftype = atoi(optarg);
                if((ftype >14)||(ftype<0)) {
                    cout << "ftype must be an integer between 1 to 13" << endl;
                    return false;
                }
                break;

            case 't': threshold = atof(optarg);
                if (threshold < 0.0) {
                    cout << "Threshold must be non-negative" << endl;
                    return false;
                }
                break;

            case 'd': C_thresh = atof(optarg);
                if (C_thresh < 0.0) {
                    cout << "Threshold must be non-negative" << endl;
                    return false;
                }
                break;

            case 'm': minGraphSize = atol(optarg);
                if(minGraphSize < 0) {
                    cout << "minGraphSize must be non-negative" << endl;
                    return false;
                }
                break;

            case 'h': minPhaseSize = atol(optarg);
            if(minPhaseSize < -1) {
                cout << "minGraphSize must be -1 or higher" << endl;
                return false;
            }
            break;
                
            case 'n': numColors = atol(optarg);
                if((numColors <1) || ((numColors >1024))) {
                    cout << "Number of colors should be between 1 and 1024 (default=16)" << endl;
                    numColors = 16;
                }
                break;
                
            case 'z': vtxOrderingAlgorithm = atoi(optarg);
                if((vtxOrderingAlgorithm > 6)||(vtxOrderingAlgorithm<1)) {
                    cout << "Vertex ordering parameter should be between 0 and 4 (default=0)" << endl;
                    cout << "You entered: " <<vtxOrderingAlgorithm << endl;
                    return false;
                }
                break;
            case 'k': nParts = atoi(optarg);
                if((nParts > 100000)||(nParts<1)) {
                    cout << "Number of partitions should be between 1 and 100,000 (default=16)" << endl;
                    cout << "You entered: " <<nParts << endl;
                    return false;
                }
                break;
            default:
                cerr << "Unknown argument" << endl;
                return false;
        }
        opt = getopt(argc, argv, opt_string);
    }

    if (argc - optind != 1) {
        cout << "Problem name not specified.  Exiting." << endl;
        usage();
        return false;
    } else {
        inFile = argv[optind];
    }

#ifdef PRINT_DETAILED_STATS_
    cout << "********************************************"<< endl;
    cout << "Input Parameters: \n";
    cout << "********************************************"<< endl;
    cout << "Input File      : " << inFile << endl;
    cout << "File type       : " << ftype  << endl;
    cout << "Threshold       : " << threshold << endl;
    cout << "C-threshold     : " << C_thresh << endl;
    cout << "Min-size        : " << minGraphSize << endl;
    cout << "Min-phase       : " << minPhaseSize << endl;
    cout << "basicOpt        : " << basicOpt << endl;
    cout << "SyncType        : " << syncType << endl;
    cout << "Percentage      : " << percentage << endl;
    cout << "# of colors     : " << numColors << endl;
    cout << "Compute_metrics : " << compute_metrics << endl;
    cout << "Vertex ordering : " << vtxOrderingAlgorithm << endl;
    cout << "--------------------------------------------" << endl;
    if (coloring) {
        cout << "Coloring   : TRUE" << endl;
        if (coloring == 1)
            cout << "Coloring Type  : distance-1 coloring" << endl;
        if (coloring == 2)
            cout << "Coloring Type  : distance-1 coloring with relabalancing" << endl;
        if (coloring == 3)
            cout << "Coloring Type  : incomplete distance-1 coloring with "<< numColors <<" colors"<< endl;
    }
    else
        cout << "Coloring   : FALSE" << endl;
    if (basicOpt)
        cout << "Replace map : TRUE" << endl;
    else
        cout << "Replace map : FALSE" << endl;
    if (strongScaling)
        cout << "Strong scaling : TRUE" << endl;
    else
        cout << "Strong scaling : FALSE" << endl;
    if(VF)
        cout << "VF         : TRUE" << endl;
    else
        cout << "VF         : FLASE" << endl;
    if(output)
        cout << "Output     : TRUE"  << endl;
    else
        cout << "Output     : FALSE"  << endl;
    if (compute_metrics)
        cout << "Compute metrics : TRUE" << endl;
    else
        cout << "Compute metrics : FALSE" << endl;
    cout << "********************************************"<< endl;
#endif

    return true;
}
