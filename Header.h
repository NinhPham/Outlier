#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

#include <vector>

using namespace std;

#define PI				3.14159265359
#define LAMBDA          0.71828182846
#define EPSILON         0.1
#define DELTA           0.1
#define MERSENNE_PRIME_31  2147483647 // 2^31 - 1
#define HL              31 // 2^31 - 1 is a Mersenne prime
/*
PARAM of data set
- PARAM_DATA_N: number of points
- PARAM_DATA_D: number of dimensions
*/
extern size_t PARAM_DATA_N; // Number of data points
extern size_t PARAM_DATA_D; // Number of dimensions
extern vector<double> POINTS; // vector of data points

/*
PARAM of distance-based algorithms, including kNN, kNNW, LOF
- PARAM_kNN_K: kNN
- PARAM_kNNW_K: kNNW
- PARAM_LOF_K: LOF
- PARAM_1NN_S: One time sampling
*/
extern size_t PARAM_KNN_K;      // kNN
extern size_t PARAM_KNNW_K;     // kNNW
extern size_t PARAM_LOF_K;      // LOF
extern size_t PARAM_1NN_S;      // One time sampling for 1NN
/*
PARAM of VOA algorithms, including VOA, RP_VOA, SamplingVOA
- PARAM_VOA_RP_R: number of random projections
- PARAM_AMS_S1: AMS Sketches size S1
- PARAM_AMS_S2: AMS Sketches size S2
- PARAM_VOA_C2N: (PARAM_DATA_N - 1) * (PARAM_DATA_N - 2) / 2
- PARAM_VOA_C2R: PARAM_VOA_RP_R * (PARAM_VOA_RP_R - 1) / 2
- PARAM_VOA_S: sample size of SamplingVOA
*/

extern size_t PARAM_VOA_RP_R; // Number of random projections
extern size_t PARAM_AMS_S1; // AMS params
extern size_t PARAM_AMS_S2; // AMS params
extern double PARAM_VOA_C2N; // (NUM_POINT - 1) * (NUM_POINT - 2) / 2;
extern double PARAM_VOA_C2R; // PARAM_VOA_RP_R * (PARAM_VOA_RP_R - 1) / 2;
extern size_t PARAM_VOA_S;


/*
PARAM of ABOF algorithms, including ABOF, SamplingABOF
- PARAM_ABOF_K: K in approxABOF
- PARAM_ABOF_S: samples size of sampling ABOF
*/
extern size_t PARAM_ABOF_K;
extern size_t PARAM_ABOF_S;

/*
PARAM of L1D algorithms, including L1Depth, SamplingL1Depth
- PARAM_L1D_S: sample size of basic sampling
- PARAM_SAMDEPTH_S: sample size of SamDepth
- PARAM_FAST_SAMDEPTH_S: sample size of Fast SamDepth
*/
extern size_t PARAM_L1D_S;
extern size_t PARAM_SAMDEPTH_S;
extern size_t PARAM_FAST_SAMDEPTH_S;

/* Internal param
- PARAM_INTERNAL_SAVE: save to disk
*/
extern bool PARAM_INTERNAL_SAVE; // save the results


typedef pair<size_t, double> IDPair;
typedef pair<size_t, size_t> IIPair;
typedef vector<double> DVector;
typedef vector<size_t> IVector;

bool ascentSortBySec(const IDPair &, const IDPair &);
bool descentSortBySec(const IDPair &, const IDPair &);

#endif // HEADER_H_INCLUDED
