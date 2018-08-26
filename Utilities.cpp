#include "Utilities.h"
#include "Header.h"

#include <random>
#include <fstream> // fscanf, fopen, ofstream

/**
Get position from 1D vector

Input:
- p_iPoint : the point idx
- p_iDim : the dimension idx

Output: the position corresponding to 1D vector
**/
size_t get1DIndex(size_t p_iPoint, size_t p_iDim)
{
    return p_iPoint * PARAM_DATA_D + p_iDim;
}

/**
Get vector of point

Input:
- A : the point index

Output: A vector of point of size D
**/
vector<double> getPoint(size_t A)
{
    return vector<double>(POINTS.begin() + A * PARAM_DATA_D, POINTS.begin() + (A + 1) * PARAM_DATA_D);
}

/**
Generate a random pair (B, C) in [N] where B != C != A.

Input:
- An integer A

Output:
- A random pair of integer (B, C) in range [N]

**/
void generateRndPair(const size_t A, size_t &B, size_t &C)
{
    do
    {
        B = rand() % PARAM_DATA_N;
        C = rand() % PARAM_DATA_N;
    }
	while (A == B || A == C || B == C);
}

/**
Randomly sampling a subset of size K points
- Use Fisher-Yates shuffle algorithm
**/
vector<size_t> samplingSubset(vector<size_t> vectorIndex, size_t K)
{
    vector<size_t>::iterator iterFirst, iterRandom;
    iterFirst = vectorIndex.begin();
    size_t left = PARAM_DATA_N - 1;

    while (K--)
    {
        //cout << *iterFirst << endl;
        iterRandom = iterFirst;

        // Find a random position in [n]
        advance(iterRandom, rand() % left);
        //cout << *iterRandom << endl;
        // Swap value
        swap(*iterFirst, *iterRandom);
        //cout << *iterFirst << endl;
        //cout << *iterRandom << endl;

        // Increase the iterFirst
        ++iterFirst;

        // Decrease the size of vector
        --left;
    }

    return vector<size_t>(vectorIndex.begin(), iterFirst);
}

/**
Generate a vector of {+1, -1} 4-wise independent elements (SODA'04)

Note:
- We need to use *prng* object and initialize it with a very large random number.
- Then we generate large random seeds to generate 4-wise random variables.

Output:
- FOUR_WISE_LEFT, FOUR_WISE_RIGHT

**/
void generate4Wises(vector<int> &FOUR_WISE_LEFT, vector<int> &FOUR_WISE_RIGHT)
{
	// vectorRan contains random seeds to use in pairwise generator hash31()
	prng_type * prng;
	vector<int> vectorRan(8, 0);

    size_t i;
    long long x;

    // Generate a large random integer
    long long lLargeNumber = rand();
    lLargeNumber =  lLargeNumber ^ ((rand() & 1) << 31);

    prng = prng_Init(lLargeNumber, 2);

	for ( i = 0; i < 8; i++ )
        //vectorRan[i] = (long long)rand();
        vectorRan[i] = (long long) prng_int(prng);

    //printVector(vectorRan);
    prng_Destroy(prng);

	long long lValue;
	for ( i = 0; i < PARAM_DATA_N; i++)
	{
		x = (long long)(i + 1);

		// Left vector
		lValue	= hash31(hash31(hash31(x, vectorRan[0], vectorRan[1]), vectorRan[2], x), vectorRan[3], x);
		FOUR_WISE_LEFT[i] = 2 * (lValue & 1) - 1;

		// Right vector
		lValue	= hash31(hash31(hash31(x, vectorRan[4], vectorRan[5]), vectorRan[6], x), vectorRan[7], x);
		FOUR_WISE_RIGHT[i] = 2 * (lValue & 1) - 1;
	}

	//printVector(FOUR_WISE_LEFT);
	//printVector(FOUR_WISE_RIGHT);

}

/**
Hash function for 4-wise independent value generator (SODA 04)

Input:
- a, b, x to compute a*x + b

Output: A 4-wise long long value
**/
long long hash31(long long a, long long b, long long x)
{
  long long result;

// return a hash of x using a and b mod (2^31 - 1)
// may need to do another mod afterwards, or drop high bits
// depending on d, number of bad guys
// 2^31 - 1 = 2147483647

  //  result = ((long long) a)*((long long) x)+((long long) b);
  result	= (a * x) + b;
  result	= ((result >> HL) + result) & MERSENNE_PRIME_31;

  return result;
}

/**
Generate Normal random variables

Output: vectorNormal
**/
void generateNormal(vector<double> &vectorNormal)
{
	//prng_type * prng;
	//prng = prng_Init(generateRndInt(), 2);

	default_random_engine generator;
    normal_distribution<double> distribution(0.0, 1.0);

	for ( size_t iRan = 0; iRan < vectorNormal.size(); iRan++ )
	{
		//vectorNormal.push_back(prng_normal(prng));
		vectorNormal[iRan] = distribution(generator);
		// printVector(vectorNormal);
	}

	//prng_Destroy(prng);
}

/**
DOT PRODUCT: dot(AB, AC)

Input:
- AB, AC: vectors to compute their dot product

Output:
- The dot product dot(AB, AC)
**/

double dotProduct(const vector<double> &AB, const vector<double> &AC)
{
	double dInner = 0.0;

	for ( size_t d = 0; d < PARAM_DATA_D; d++ )
		dInner += AB[d] * AC[d];

	return dInner;
}

/**
Euclidean distance between 2 vectors

Input:
- AB, AC: vectors to compute their dot product

Output:
- The Eulidean distance between AB and AC
**/

double l2Dist(const vector<double> &AB, const vector<double> &AC)
{
	double dDist = 0.0;

	for ( size_t d = 0; d < PARAM_DATA_D; d++ )
		dDist += pow(AB[d] - AC[d], 2);

	return sqrt(dDist);
}

/**
Get k nearest neighbors for the point A

Input:
- A : the index of point A
- K : k value in kNN

Output:
- vector_kNN : the sorted vector containing k nearest neighbor of A

**/
void kNearestNeighbors(vector<IDPair> & vectorkNN, size_t A, const size_t p_iK)
{
	vector<IDPair>::iterator iLow;
	double dLength = 0.0;
	size_t B, iDim;

	vectorkNN.clear();

	for ( B = 0; B < PARAM_DATA_N; B++ )
	{
	    if ( B!= A)
        {
            // Compute vector AB and its length
            dLength = 0.0;
            for ( iDim = 0; iDim < PARAM_DATA_D; iDim++ )
                dLength	+= pow(POINTS[get1DIndex(B, iDim)] - POINTS[get1DIndex(A, iDim)], 2);

            dLength = sqrt(dLength);

            IDPair tempPair = make_pair(B, dLength);

            // vector_kNN is always sorted in ascent order
            iLow = lower_bound(vectorkNN.begin(), vectorkNN.end(), tempPair, ascentSortBySec);

            if ( vectorkNN.size() < p_iK )
                vectorkNN.insert(iLow, tempPair); // Insert when there are not enough k nearest neighbors
            else
            {
                // Remove the furthest point when there are enough k nearest neighbors
                if ( iLow != vectorkNN.end() )
                {
                    vectorkNN.insert(iLow, tempPair);
                    vectorkNN.pop_back();
                }
            }
        }

	}
}

/**
Get the k nearest neighbors distance for the point A

Input:
- A : the index of point A
- K : k value in kNN

Output: the kNN distance of A
**/
double kNN_Dist(size_t A, size_t p_iK)
{
	vector<double>::iterator iLow;
	double dLength = 0.0;
	size_t B, iDim;

	vector<double> vectorkNN;

	for ( B = 0; B < PARAM_DATA_N; B++ )
	{
		// Compute vector AB and its length
		if (B != A)
        {
            dLength = 0.0;
            for ( iDim = 0; iDim < PARAM_DATA_D; iDim++ )
                dLength	+= pow(POINTS[get1DIndex(B, iDim)] - POINTS[get1DIndex(A, iDim)], 2);

            dLength = sqrt(dLength);

            // vector_kNN is always sorted in ascent order
            iLow = lower_bound(vectorkNN.begin(), vectorkNN.end(), dLength);

            if ( vectorkNN.size() < p_iK )
                vectorkNN.insert(iLow, dLength); // Insert when there are not enough k nearest neighbors
            else
            {
                // Remove the furthest point when there are enough k nearest neighbors
                if ( iLow != vectorkNN.end() )
                {
                    vectorkNN.insert(iLow, dLength);
                    vectorkNN.pop_back();
                }
            }
        }

	}

	return vectorkNN[p_iK - 1];
}

/**
Print a vector
**/
void printVector(const vector<double> &vecPrint)
{
	cout << "Vector is: ";
	for (size_t i = 0; i < vecPrint.size(); i++)
		cout << vecPrint[i] << " ";

	cout << endl;
}

/**
Print a vector
**/
void printVector(const vector<int> &vecPrint)
{
	cout << "Vector is: ";
	for (size_t i = 0; i < vecPrint.size(); i++)
		cout << vecPrint[i] << " ";

	cout << endl;
}

/**
Print a vector
**/
void printVector(const vector<size_t> &vecPrint)
{
	cout << "Vector is: ";
	for (size_t i = 0; i < vecPrint.size(); i++)
		cout << vecPrint[i] << " ";

	cout << endl;
}

/**
Save a vector
**/
void saveVector(const vector<int> &vecResult, string fileName)
{
	ofstream myfile;
    myfile.open(fileName);
    for (size_t i = 0; i < vecResult.size(); i++)
        myfile << vecResult[i] << endl;

    myfile.close();
}

/**
Save ranking: [pointID outlier_score]
**/
void saveRanking(const vector<IDPair> &vecResult, string fileName)
{
    if (!PARAM_INTERNAL_SAVE)
        return;

    ofstream myfile;
    myfile.open(fileName);
    for (size_t i = 0; i < vecResult.size(); i++)
    {
        myfile << vecResult[i].first + 1;
        myfile << " ";
        myfile << vecResult[i].second << endl;
    }

    myfile.close();
}

/**
Save all info: [pointID F1 F2 VOA]
**/
void saveVOAInfo(const vector<double> &F1, const vector<double> &F2, const vector<double> &VOA, string fileName)
{
    if (!PARAM_INTERNAL_SAVE)
        return;

    ofstream myfile;
    myfile.open(fileName);
    for (size_t i = 0; i < PARAM_DATA_N; i++)
    {
        myfile << (i + 1);
        myfile << " ";
        myfile << F1[i];
        myfile << " ";
        myfile << F2[i];
        myfile << " ";
        myfile << VOA[i] << endl;
    }

    myfile.close();
}



/**
Truncate the cosine value to be in [-1, +1]
**/
double truncateCosine(double x)
{
    return max(-1.0, min(x, 1.0));
}

/**
Approximate acos() by Nvidia Toolkits
https://web.archive.org/web/20161223122122/http://http.developer.nvidia.com:80/Cg/acos.html
**/
double hackAcos(double x)
{
    double negate = double(x < 0);
    x = abs(x);
    double ret = -0.0187293;
    ret = ret * x;
    ret = ret + 0.0742610;
    ret = ret * x;
    ret = ret - 0.2121144;
    ret = ret * x;
    ret = ret + 1.5707288;
    ret = ret * sqrt(1.0-x);
    ret = ret - 2 * negate * ret;
    return negate * 3.14159265358979 + ret;
}

/**
Compute all information, including |AB|, |AC|, <AB, AC>

Input:
    - A, B, C: point indexes

Output:
    - p_dInner = <AB, AC>
    - p_dLengthAB = |AB|
    - p_dLengthAC = |AC|
**/
void computeAB_AC(const size_t A, const size_t B, const size_t C, double &p_dInner, double &p_dLengthAB, double &p_dLengthAC)
{
    vector<double> vectorA = getPoint(A);
    vector<double> vectorB = getPoint(B);
    vector<double> vectorC = getPoint(C);

    vector<double> vectorAB(PARAM_DATA_D, 0.0), vectorAC(PARAM_DATA_D, 0.0);

    p_dLengthAB = 0.0;
    p_dLengthAC = 0.0;
    p_dInner = 0.0;

    double dTemp1, dTemp2;

    for ( size_t iDim = 0; iDim < PARAM_DATA_D; iDim++ )
    {
        dTemp1 = vectorA[iDim] - vectorB[iDim];
        p_dLengthAB += dTemp1 * dTemp1;

        dTemp2 = vectorA[iDim] - vectorC[iDim];
        p_dLengthAC += dTemp2 * dTemp2;

        p_dInner += dTemp1 * dTemp2;
    }

    p_dLengthAB = sqrt(p_dLengthAB);
    p_dLengthAC = sqrt(p_dLengthAC);

}

/**
Convert an integer to string
**/
string int2str(int x)
{
    stringstream ss;
    ss << x;
    return ss.str();
}
