/**
This C++ source code is written by Ninh Pham (pham@di.ku.dk) as a part of the DABAI project
We used the random generator from "http://www.cs.rutgers.edu/~muthu/massdal-code-index.html"
Feel free to re-use and re-distribute this source code for any purpose,
and cite our work when you re-use this code.
**/

//#include "Header.h"
//#include "Test.h"
#include "InputParser.h"
#include "Utilities.h"

#include "L1Depth.h"
#include "VOA.h"
#include "ABOF.h"
#include "Distance.h"
#include "RP_VOA.h"
#include "SamplingL1Depth.h"
#include "SamplingKNN.h"

#include <limits> // for max and min of integer
#include <iostream> // cin, cout
#include <time.h> // for time(0) to generate different random number
#include <stdlib.h>

/**
Parameters:
- PARAM_DATA_N  : number of points (N)
- PARAM_DATA_D  : number of dimensions (D)
- FILE_NAME     : filename of dataset with the matrix format N x D
- METHOD_NAME   : used method

"kNN": SIGMOD 00
    - PARAM_KNN_K: number k in kNN
"Sam1NN" : NIPS 13
    - PARAM_1NN_S: number of samples

"kNNW": TKDE 05
    - PARAM_KNNW_K: number k in kNN
"LOF": SIGMOD 00
    - PARAM_LOF_K: minPts

"ABOF": KDD 08
    - No parameter
"approxABOF": KDD 08
    - PARAM_ABOF_K: number k in kNN

"VOA": KDD 12
    - No parameter
"FastVOA": KDD 12
    - PARAM_RP_VOA_R: number of random projections
    - PARAM_AMS_S1: AMS Sketch size
    - PARAM_AMS_S2: AMS Sketch size

"L1D": PKDD 18
    - No parameter
"BasicSamL1D": PKDD 18
    - PARAM_L1D_S: number of sample pairs
"SamDepth": PKDD 18
    - PARAM_SAMDEPTH_S: number of samples
**/
int main(int nargs, char** args)
{

    srand(time(NULL)); // should only be called once

	/************************************************************************/
	/* Load input file                                                      */
	/************************************************************************/
	int iType = loadInput(nargs, args);

	// Set up some internal parameters
	PARAM_INTERNAL_SAVE = true; // saving results

	/************************************************************************/
	/* Approaches                                             */
	/************************************************************************/

	switch (iType)
	{
    /** All Exact solutions: kNN, kNNW, LOF, ABOF, VOA, MOC **/
    case 0:
        kNN();
        kNNW();
        LOF();
        ABOF();
        VOA();
        L1_Depth();
        break;

    /** Global outliers **/
    case 11:
        kNN();
        break;
    case 12:
        kNNW();
        break;
    case 13:
        Sam1NN();
        break;

    /** Local outliers **/
    case 21:
        LOF();
        break;

    /** ABOF **/
    case 31:
        ABOF();
        break;
    case 32:
        approxABOF();
        break;

    /** VOA **/
    case 41:
        VOA();
        break;
    case 42:
        FastVOA();
        break;

    /** L1 Depth **/
    case 51:
        L1_Depth();
        break;
    case 52:
        BasicSampling_L1D();
        break;
    case 53:
        SamDepth();
        break;

    }
	// system("PAUSE");
	return 0;
}
