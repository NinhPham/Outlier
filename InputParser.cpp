#include "InputParser.h"
#include "Header.h"
#include "Utilities.h"

#include "stdlib.h"     /* atoi */
#include <iostream> // cin, cout
#include <fstream> // fscanf, fopen, ofstream

#include "string.h" // strcmp is C lib
#include <math.h> // ceil


using namespace std;

// declare PARAM for data and algorithm

/*
PARAM of data set
- PARAM_DATA_N: number of points
- PARAM_DATA_D: number of dimensions
*/
size_t PARAM_DATA_N; // Number of data points
size_t PARAM_DATA_D; // Number of dimensions
vector<double> POINTS; // vector of data points

/*
PARAM of distance-based algorithms, including kNN, kNNW, LOF
- PARAM_kNN_K: kNN
- PARAM_kNNW_K: kNNW
- PARAM_LOF_K: LOF
- PARAM_1NN_S; One time sampling 1NN
*/
size_t PARAM_KNN_K; // kNN search
size_t PARAM_KNNW_K; // kNN search
size_t PARAM_LOF_K; // kNN search
size_t PARAM_1NN_S; // One time sampling 1NN

/*
PARAM of VOA algorithms, including VOA, RP_VOA, SamplingVOA
- PARAM_VOA_RP_R: number of random projections
- PARAM_AMS_S1: AMS Sketches size S1
- PARAM_AMS_S2: AMS Sketches size S2
- PARAM_VOA_S: sample size of SamplingVOA
- PARAM_VOA_C2N: (PARAM_DATA_N - 1) * (PARAM_DATA_N - 2) / 2
- PARAM_VOA_C2R: PARAM_VOA_RP_R * (PARAM_VOA_RP_R - 1) / 2
- PARAM_VOA_NORMAL: Normal distribution of size PARAM_VOA_RP_R * PARAM_DATA_D
*/
size_t PARAM_VOA_RP_R;
size_t PARAM_AMS_S1;
size_t PARAM_AMS_S2;
size_t PARAM_VOA_S;
double PARAM_VOA_C2N;
double PARAM_VOA_C2R;
vector<double> PARAM_VOA_NORMAL;

/*
PARAM of ABOF algorithms, including ABOF, SamplingABOF
- PARAM_ABOF_K: K in approxABOF
- PARAM_ABOF_S: samples size of sampling ABOF
*/
size_t PARAM_ABOF_K;
size_t PARAM_ABOF_S;

/*
PARAM of L1D algorithms, including L1Depth, SamplingL1Depth
- PARAM_L1D_S: sample size of basic sampling
- PARAM_SAMDEP_S: sample size of SamDepth
- PARAM_FAST_SAMDEPTH_S: sample size of Fast SamDepth
*/
size_t PARAM_L1D_S;
size_t PARAM_SAMDEPTH_S;
size_t PARAM_FAST_SAMDEPTH_S;

int loadInput(int nargs, char** args)
{
    if (nargs < 4)
        exit(1);

    // Parse arguments: Note that don't know what args[0] represents for !!!
    PARAM_DATA_N = atoi(args[1]);
    cout << "Number of points: " << PARAM_DATA_N << endl;

    PARAM_VOA_C2N = ((double)PARAM_DATA_N - 1) * ((double)PARAM_DATA_N - 2) / 2.0;

    PARAM_DATA_D = atoi(args[2]);
    cout << "Number of dimensions: " << PARAM_DATA_D << endl;

    // Read data
    size_t idxPoint, idxDim;
    cout << "Read data..." << endl;
    if (args[3])
    {
        FILE *f = fopen(args[3], "r+");
        if (!f)
        {
            printf("Data file does not exist");
            exit(1);
        }

        POINTS = vector<double>(PARAM_DATA_N * PARAM_DATA_D, 0.0);
        for (idxPoint = 0; idxPoint < PARAM_DATA_N; idxPoint++)
        {
            for (idxDim = 0; idxDim < PARAM_DATA_D; idxDim++)
            {
                fscanf(f, "%lf", &POINTS[get1DIndex(idxPoint, idxDim)]);
                //cout << POINTS[get1DIndex(idxPoint, idxDim)] << " ";
            }
            //cout << endl;
        }

        fclose(f);
        delete f;
    }

    // Algorithm
    int iType = 0;

    /**
    Internal function: All deterministic solutions
    **/
    if (strcmp(args[4], "exact") == 0)
    {
        PARAM_KNN_K = atoi(args[5]);
        PARAM_KNNW_K = atoi(args[6]);
        PARAM_LOF_K = atoi(args[7]);

        iType = 0;
        cout << "KNN with K = " << PARAM_KNN_K << "... " << endl;
        cout << "KNNW with K = " << PARAM_KNNW_K << "... " << endl;
        cout << "LOF with K = " << PARAM_LOF_K << "... " << endl;
    }

    /**
    kNN distance
    **/
    else if (strcmp(args[4], "kNN") == 0)
    {
        PARAM_KNN_K = atoi(args[5]);

        iType = 11;
        cout << "Using kNN with K = " << PARAM_KNN_K << "... " << endl;
    }

    /**
    Accumulate kNN distances
    **/
    else if (strcmp(args[4], "kNNW") == 0)
    {
        PARAM_KNNW_K = atoi(args[5]);

        iType = 12;
        cout << "Using kNNW with K = " << PARAM_KNNW_K << "... " << endl;
    }

    /**
    Accumulate kNN distances
    **/
    else if (strcmp(args[4], "Sam1NN") == 0)
    {
        PARAM_1NN_S = atoi(args[5]);

        iType = 13;
        cout << "Using one time sampling for 1NN with S = " << PARAM_1NN_S << "... " << endl;
    }

    /**
    LOF
    **/
    else if (strcmp(args[4], "LOF") == 0)
    {
        PARAM_LOF_K = atoi(args[5]);

        iType = 21;
        cout << "Using LOF with K = " << PARAM_LOF_K << "... " << endl;
    }

    /**
    ABOF
    **/
    else if (strcmp(args[4], "ABOF") == 0)
    {
        iType = 31;
        cout << "Using naive ABOF... " << endl;
    }

    /**
    FastABOF
    **/
    else if (strcmp(args[4], "approxABOF") == 0)
    {
        PARAM_ABOF_K = atoi(args[5]);

        iType = 32;
        cout << "Using approximate ABOD with K = " << PARAM_ABOF_K << "... " << endl;
    }

    /**
    Important sampling approach for ABOF, using O(N) samples
    **/
    else if (strcmp(args[4], "SamABOF") == 0)
    {
        PARAM_ABOF_S = atoi(args[5]);

        iType = 33;
        cout << "Basic importance sampling approach for ABOF with " << PARAM_ABOF_S << " pairs" << endl;
    }

    /**
    VOA
    **/
    else if (strcmp(args[4], "VOA") == 0)
    {
        iType = 41;
        cout << "Using naive VOA... " << endl;
    }

    /**
    FastVOA - KDD 12
    **/
    else if (strcmp(args[4], "FastVOA") == 0)
    {
        PARAM_VOA_RP_R = atoi(args[5]);
        PARAM_VOA_C2R = double(PARAM_VOA_RP_R) * (double(PARAM_VOA_RP_R) - 1) / 2.0;

        PARAM_AMS_S1 = atoi(args[6]);
        PARAM_AMS_S2 = atoi(args[7]);

        iType = 42;
        cout << "Using Fast VOA with "<< PARAM_VOA_RP_R << " projections and  AMS size " \
             << PARAM_AMS_S1 << "x" << PARAM_AMS_S2 << "... " << endl;
    }

    /**
    Fast RP for VOA - Avoid AMS Sketch to obtain a quadratic time
    **/
    else if (strcmp(args[4], "FastRP_VOA") == 0)
    {
        PARAM_VOA_RP_R = atoi(args[5]);
        PARAM_VOA_C2R = double(PARAM_VOA_RP_R) * (double(PARAM_VOA_RP_R) - 1) / 2.0;

        iType = 43;
        cout << "Using Random Projections with "<< PARAM_VOA_RP_R << " projections and without AMS... " << endl;
    }

    /**
    Basic sampling approach for VOA, using unbiased variance (MOA2 - MOA1^2) * (S - 1) / S
    **/
    else if (strcmp(args[4], "BasicSamVOA") == 0)
    {
        PARAM_VOA_S = atoi(args[5]);

        iType = 44;
        cout << "Direct sampling and sample variance approaches for VOA with " << PARAM_VOA_S << " pairs" << endl;
    }

    /**
    Adaptive sampling approach for VOA, using epsilon = 0.1, delta = 0.1
    **/
    else if (strcmp(args[4], "AdaptSamVOA") == 0)
    {
        iType = 45;
        cout << "Adaptive sampling approach for VOA..." << endl;
    }

    /**
    Fast Adaptive sampling approach for VOA, using epsilon = 0.1, delta = 0.1
    **/
    else if (strcmp(args[4], "FastAdaptSamVOA") == 0)
    {
        iType = 46;
        cout << "Adaptive sampling approach on MOA2 for VOA..." << endl;
    }

    /**
    L1 Depth
    **/
    else if (strcmp(args[4], "L1D") == 0)
    {
        iType = 51;
        cout << "Using L1 Depth..." << endl;
    }

    /**
    Basic sampling for L1D, using O(log(1/delta) * (1 / epsilon^2))
    **/
    else if (strcmp(args[4], "BasicSamL1D") == 0)
    {
        PARAM_L1D_S = atoi(args[5]);

        iType = 52;
        cout << "Basic sampling for L1D with " << PARAM_L1D_S << " pairs..." << endl;
    }

    /**
    Batch sampling for L1D
    **/
    else if (strcmp(args[4], "SamDepth") == 0)
    {
        PARAM_SAMDEPTH_S = atoi(args[5]);

        iType = 53;
        cout << "Batch sampling for L1D with " << PARAM_SAMDEPTH_S << " points..." << endl;
    }

    /**
    Fast Batch sampling for L1D
    **/
    else if (strcmp(args[4], "FastSamDepth") == 0)
    {
        PARAM_FAST_SAMDEPTH_S = atoi(args[5]);

        iType = 54;
        cout << "Fast Batch sampling for L1D with " << PARAM_FAST_SAMDEPTH_S << " points..." << endl;
    }

    /**
    Internal functions: All sampling approaches for L1Depth
    **/
    else if (strcmp(args[4], "SamL1D") == 0)
    {
        PARAM_L1D_S = atoi(args[5]);
        PARAM_SAMDEPTH_S = atoi(args[6]);

        iType = 50;
        cout << "Basic sampling approaches for L1D with " << PARAM_L1D_S << " pairs..." << endl;
        cout << "Batch sampling approaches for L1D with " << PARAM_SAMDEPTH_S << " points..." << endl;
    }

    return iType;
}

