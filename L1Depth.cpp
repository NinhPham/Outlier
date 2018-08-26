#include "L1Depth.h"
#include "Utilities.h"


/**
Ranking data based on L1 Depth factors that requires O(dn^2) - PKDD'18
**/
void L1_Depth()
{
	double dStart = clock();

	vector<IDPair> vectorL1Depth(PARAM_DATA_N); // L1_Depth
    vector<IDPair> vectorMOA2(PARAM_DATA_N); // Second moment of angle

	vector<double> vectorSum(PARAM_DATA_D, 0.0); // length of |AB|, |AC|
	vector<double> vectorEDGE(PARAM_DATA_D, 0.0); // vector AB
	vector<double> vectorA; // point A
	vector<double> vectorB; // point B

	double dTemp, dLength;
	size_t A, B, iDim;
	size_t iCountZero; // counting number of duplicates since we need to avoid

	for ( A = 0; A < PARAM_DATA_N; A++ )
	{
		vectorA = getPoint(A);

		iCountZero = 0;
        fill(vectorSum.begin(), vectorSum.end(), 0.0);

		// EDGE stores the coordinates, LENGTH stores the norm_2
		for ( B = 0; B < PARAM_DATA_N; B++ )
		{
		    if (B == A)
                continue;

            vectorB = getPoint(B);

            // Need to clear the length
            dLength = 0.0;
            fill(vectorEDGE.begin(), vectorEDGE.end(), 0.0);

            for ( iDim = 0; iDim < PARAM_DATA_D; iDim++ )
            {
                dTemp = vectorB[iDim] - vectorA[iDim];
                vectorEDGE[iDim] = dTemp; // Coordinates of vector AB

                dLength += dTemp * dTemp;
            }

            if (dLength > 0.0)
            {
                dLength = sqrt(dLength);

                // Normalize vector AB and add up to the vector sum
                for ( iDim = 0; iDim < PARAM_DATA_D; iDim++ )
                    vectorSum[iDim] += vectorEDGE[iDim] / dLength;
            }
            else
                iCountZero++;

		}

		// Compute the L1_Depth
		dLength = dotProduct(vectorSum, vectorSum);
		vectorL1Depth[A] = make_pair(A, 1 - sqrt(dLength) / (PARAM_DATA_N - 1));

		// Compute second moment
		vectorMOA2[A] = make_pair(A, (pow(PARAM_DATA_N - 1, 2) - dLength) / PARAM_VOA_C2N);
	}

	printf("L1_Depth time is %f \n", clock() - dStart);

	// Saving
	saveRanking(vectorMOA2, "L1D_MOA2_noSort.txt");
	sort(vectorMOA2.begin(), vectorMOA2.end(), ascentSortBySec);
	saveRanking(vectorMOA2, "L1D_MOA2.txt");

	saveRanking(vectorL1Depth, "L1D_noSort.txt");
	sort(vectorL1Depth.begin(), vectorL1Depth.end(), ascentSortBySec);
	saveRanking(vectorL1Depth, "L1D.txt");
}


