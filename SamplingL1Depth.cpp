#include "SamplingL1Depth.h"
#include "Utilities.h"

/**
Ranking data based on uniform sampling estimators using s = 4 * ln(2/delta) / epsilon^2 sample pairs
We sample with replacement a random pair of different points (B, C) then estimate avgCosine as outlier factors
**/
void BasicSampling_L1D()
{
	double dStart = clock();

	/************************************************************************/
	/* Compute L1Depth                             */
	/************************************************************************/
	vector<IDPair> vectorRank(PARAM_DATA_N);

	vector<double> vectorA, vectorB, vectorC;
	vector<double> vectorAB(PARAM_DATA_D, 0.0), vectorAC(PARAM_DATA_D, 0.0);

	size_t A, B, C;
	double dSumCosine; // sum of 1/(AB *AC)
	double dLengthAB, dLengthAC, dInner; // |AB|, |AC|, <AB, AC>
	double dCosine, dL1D;

    size_t s = 0;
	for ( A = 0; A < PARAM_DATA_N; A++ )
	{
		// Reset everything
        s = 0;
        dSumCosine = 0.0;

		while (s < PARAM_L1D_S)
		{
		    generateRndPair(A, B, C);
            computeAB_AC(A, B, C, dInner, dLengthAB, dLengthAC);

            // Sampling results
            if (dLengthAB == 0 || dLengthAC == 0)
                continue;
            else
            {
                s++;

                dCosine = dInner / (dLengthAB * dLengthAC);
                dCosine = truncateCosine(dCosine);

                // Accumulate cos(\theta_{apb})
                dSumCosine += dCosine;
            }
		}

		// Compute estimators as average value
        dSumCosine = dSumCosine / PARAM_L1D_S;

        // Compute L1D
        dL1D = (dSumCosine * (PARAM_DATA_N - 2) + 1)/ (PARAM_DATA_N - 1);
        dL1D = 1 - sqrt(dL1D);

        // Compute scores VOA, MOA1, MOA2 for each point
		vectorRank[A] = make_pair(A, dL1D);
	}

    printf("Basic Sampling L1 Depth - # Pairs: %d \n", PARAM_L1D_S);
    printf("Basic Sampling L1 Depth - Time: %f \n\n", clock() - dStart);

	/** Data depth **/
	saveRanking(vectorRank, "Basic_Sampling_L1D_noSort_" + int2str(PARAM_L1D_S) + ".txt");

	sort(vectorRank.begin(), vectorRank.end(), ascentSortBySec);

	saveRanking(vectorRank, "Basic_Sampling_L1D_" + int2str(PARAM_L1D_S) + ".txt");

}

/**
Ranking data based on batch sampling estimators - PKDD'18
We sample without replacement a random subset S of size sqrt(n)
**/
void SamDepth()
{
	double dStart = clock();

	/************************************************************************/
	/* Compute L1Depth                             */
	/************************************************************************/
	vector<IDPair> vectorRank(PARAM_DATA_N);

	vector<double> vectorX, vectorA;
	vector<double> vectorEDGE(PARAM_DATA_D, 0.0), vectorSum(PARAM_DATA_D, 0.0);

    vector<size_t> vecSubset;

	size_t A, B, X;
	double dLength, iDim, dTemp, dL1D;

	// Initialize vector index of points: 1 --> (N-1)
    vector<size_t> vectorIndex(PARAM_DATA_N - 1, 0);
    for ( A = 0; A < PARAM_DATA_N - 1; A++)
        vectorIndex[A] = A + 1;

    size_t P = 0;
    size_t iDuplicate, iRealNumSamples;

    for ( A = 0; A < PARAM_DATA_N; A++ )
	{
	    vectorA = getPoint(A);

		// Replace the A by the removed element from the previous iteration, hence the vectorIndex does not contain A
		if (A > 0)
            swap(vectorIndex[A - 1], P); // note here that P = A for each iteration

	    // Reset everything
        fill(vectorSum.begin(), vectorSum.end(), 0.0);
        iDuplicate = 0;

		// Generate a random permutation then pick top S position without A
        vecSubset = samplingSubset(vectorIndex, PARAM_SAMDEPTH_S);

		for ( B = 0; B < PARAM_SAMDEPTH_S; B++ )
		{
		    X = vecSubset[B];
            vectorX = getPoint(X);

            // Need to clear the length
            dLength = 0.0;
            fill(vectorEDGE.begin(), vectorEDGE.end(), 0.0);

            for ( iDim = 0; iDim < PARAM_DATA_D; iDim++ )
            {
                dTemp = vectorA[iDim] - vectorX[iDim];
                vectorEDGE[iDim] = dTemp; // Coordinates of vector AX

                dLength += dTemp * dTemp;
            }

            dLength = sqrt(dLength);

            if (dLength == 0)
                iDuplicate++;
            else
            {
                // Normalize vector AB and add up to the vector sum
                for ( iDim = 0; iDim < PARAM_DATA_D; iDim++ )
                    vectorSum[iDim] += vectorEDGE[iDim] / dLength;
            }

		}

		// Compute the L1_Depth
		dLength = dotProduct(vectorSum, vectorSum);
		iRealNumSamples = PARAM_SAMDEPTH_S - iDuplicate; // t in the paper formula

		dL1D = (dLength - iRealNumSamples) / ((iRealNumSamples - 1) * iRealNumSamples);
		dL1D = (dL1D * (PARAM_DATA_N - 2) + 1)/ (PARAM_DATA_N - 1);
		dL1D = 1 - sqrt(dL1D);

		vectorRank[A] = make_pair(A, dL1D);
	}

    printf("Batch Sampling L1 Depth - # Samples: %d \n", PARAM_SAMDEPTH_S);
    printf("Batch Sampling L1 Depth - Time: %f \n\n", clock() - dStart);

	/** Data depth **/
    saveRanking(vectorRank, "SamDepth_noSort_" + int2str(PARAM_SAMDEPTH_S) + ".txt");

	sort(vectorRank.begin(), vectorRank.end(), ascentSortBySec);

	saveRanking(vectorRank, "SamDepth_" + int2str(PARAM_SAMDEPTH_S) + ".txt");
}



