#include "SamplingKNN.h"
#include "Utilities.h"

/**
Ranking data based on batch sampling estimators - NIPS 13
We sample without replacement a random subset S of size 20
**/
void Sam1NN()
{
	double dStart = clock();

	/************************************************************************/
	/* Compute L1Depth                             */
	/************************************************************************/
	vector<IDPair> vectorRank(PARAM_DATA_N);
	vector<double> vectorB, vectorA;
	vector<size_t> vecSubset;

	size_t A, B;
	double d1NN_Dist;

	// Subsampling a set
    vector<size_t> vectorIndex(PARAM_DATA_N, 0);
    for ( A = 0; A < PARAM_DATA_N; A++)
        vectorIndex[A] = A;
    vecSubset = samplingSubset(vectorIndex, PARAM_1NN_S);

    //printVector(vecSubset);

    for ( A = 0; A < PARAM_DATA_N; A++ )
	{
	    vectorA = getPoint(A);
        d1NN_Dist = (double)MERSENNE_PRIME_31; // suppose 1NNDist is very large 1^31 - 1 at the beginning

		// Compute 1NN distance over the subset
        for ( B = 0; B < PARAM_1NN_S; B++)
        {
            vectorB = getPoint(vecSubset[B]);
            d1NN_Dist = min(d1NN_Dist, l2Dist(vectorA, vectorB));
        }

		vectorRank[A] = make_pair(A, d1NN_Dist);
	}

    printf("One Time Sampling 1NN - # Samples: %d \n", PARAM_1NN_S);
    printf("One Time Sampling 1NN - Time: %f \n\n", clock() - dStart);

	/** Data depth **/
    saveRanking(vectorRank, "Sam1NN_noSort_" + int2str(PARAM_1NN_S) + ".txt");

	sort(vectorRank.begin(), vectorRank.end(), descentSortBySec);

	saveRanking(vectorRank, "Sam1NN_" + int2str(PARAM_1NN_S) + ".txt");
}


