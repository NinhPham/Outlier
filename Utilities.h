#ifndef UTILITIES_H_INCLUDED
#define UTILITIES_H_INCLUDED

#include "Header.h"
#include "prng.h"

//#include <string>
#include <sstream> // already included string
#include <utility> // for pair<int, int>

#include "time.h" // clock()
#include <algorithm> // sort()
#include <cstdio> // printf()
#include <math.h> // log()
#include <iostream> // cout

//const double Gamma = ceil(4 * LAMBDA * log(2 / DELTA) / pow(EPSILON, 2));

/**
Manipulate 2D vector by 1D vector
**/
size_t get1DIndex(size_t, size_t); // Get index from (pointIdx, dimIDx)
vector<double> getPoint(size_t); // Get vector of point

/**
KNN operations
**/
void kNearestNeighbors(vector<IDPair> &, const size_t, const size_t); // get kNN list of (pointIdx, distance)
double kNN_Dist(const size_t, const size_t); // Get kNN distance

/**
Generate random uniform variables
**/
void generateRndPair(const size_t, size_t &, size_t &); // generate of random pair (B, C) where A != B != C
vector<size_t> samplingSubset(vector<size_t>, size_t ); // subsampling a subset from [N]
void generateNormal(vector<double> &); // generate Normal distribution variables
void generate4Wises(vector<int> &, vector<int> &); // generate 2 4-wise independent vectors
long long hash31(long long, long long, long long); // generate 2-wise independent h(x) = ax + b mod 2^31

/**
Linear Algebra operations
**/
double dotProduct(const vector<double> &, const vector<double> &);
double l2Dist(const vector<double> &, const vector<double> &);
void computeAB_AC(const size_t, const size_t, const size_t , double &, double &, double &);

/**
Print vector in command line
**/
void printVector(const vector<double> &);
void printVector(const vector<int> &);
void printVector(const vector<size_t> &);

/**
Type conversion
**/
string int2str(int);

/**
Save outputs
**/
void saveRanking(const vector<IDPair> &, string);
void saveVOAInfo(const vector<double> &, const vector<double> &, const vector<double> &, string);
void saveVector(const vector<int> &, string );

/**
Arithmetic operations
**/
double truncateCosine(double);
double hackAcos(double);

#endif // UTILITIES_H_INCLUDED
