#include "ABOF.h"
#include "Utilities.h"


/**
Compute the ABOD for any point A (KDD'08)

Input:
- vectorEDGES: AB, AC, ...
- vectorLENGTHS: |AB|, |AC|, ...

Output:
- ABOF: variance of cosine(AB, AC)/ (|AB|^2 * |AC|^2) with weight 1 / (|AB| * |AC|)
- ABOF1: the first moment
- ABOF2: the second moment
- Note: We need ABOF1, ABOF2 to measure the accuracy of sampling methods
**/
double computeABOF(const vector<DVector> &vectorEDGES, const vector<double> &vectorLENGTHS, double &ABOF1, double &ABOF2)
{
	double dABOF = 0.0, dWeight = 0.0;
	double dABOFSum = 0.0, dABOFSum2 = 0.0, dWeightSum = 0.0;

	vector<double> AB, AC;
	size_t B, C;

	for ( B = 0; B < PARAM_DATA_N; B++ )
	{
		if ( vectorLENGTHS[B] == 0.0 )
			continue;

		AB = vectorEDGES[B];
		for ( C = B + 1; C < PARAM_DATA_N; C++)
		{
			if ( vectorLENGTHS[C] == 0.0 )
				continue;

			AC = vectorEDGES[C];

			dABOF	= dotProduct(AB, AC) / pow(vectorLENGTHS[B] * vectorLENGTHS[C], 2);

            dWeight = 1 / (vectorLENGTHS[B] * vectorLENGTHS[C]);

			dWeightSum  += dWeight;
			dABOFSum	+= dWeight * dABOF;
			dABOFSum2	+= dWeight * dABOF * dABOF;
		}
	}

	// ABOD
	ABOF1 = dABOFSum / dWeightSum;
	ABOF2 = dABOFSum2 / dWeightSum;

	return ABOF2 - ABOF1 * ABOF1;
}

/**
Approximate ABOF for any point A (FastABOD - KDD'08)

Input:
- vector_kNN : the k nearest neighbors of that point A
- The point A

Output:
- The approximation of weighted ABOF which is based on kNN
**/
double computeApproxABOF(const vector<IDPair> & vector_kNN, const size_t A)
{
	double dABOF = 0.0, dWeight = 0.0;
	double dABOFSum1 = 0.0, dABOFSum2 = 0.0, dWeightSum = 0.0;
	size_t i, j, B, C, iDim;

	vector<double> AB(PARAM_DATA_D, 0.0), AC(PARAM_DATA_D, 0.0);
	double dLengthAB, dLengthAC;

	for ( i = 0; i < PARAM_ABOF_K; i++ )
	{
		B = vector_kNN[i].first;
        dLengthAB = vector_kNN[i].second;

        if ( dLengthAB == 0 )
            continue;

		for ( j = i + 1; j < PARAM_ABOF_K; j++ )
		{
			C = vector_kNN[j].first;
			dLengthAC = vector_kNN[j].second;

			if ( dLengthAC == 0 )
                continue;

			for (iDim = 0; iDim < PARAM_DATA_D; iDim++)
			{
				AB[iDim] = POINTS[get1DIndex(B, iDim)] - POINTS[get1DIndex(A, iDim)];
				AC[iDim] = POINTS[get1DIndex(C, iDim)] - POINTS[get1DIndex(A, iDim)];
			}

			dABOF	= dotProduct(AB, AC) / pow(dLengthAB * dLengthAC, 2);
			dWeight = 1 / (dLengthAB * dLengthAC);

			dWeightSum  += dWeight;
			dABOFSum1	+= dWeight * dABOF;
			dABOFSum2	+= dWeight * dABOF * dABOF;
		}
	}

	// Variance of weighted ABOD based on kNN
	return (dABOFSum2 / dWeightSum) - pow(dABOFSum1 / dWeightSum, 2);
}

/**
Ranking data based on outlier factors ABOF that requires O(dn^3).
**/
void ABOF()
{
	double dStart = clock();

	/************************************************************************/
	/* Compute ABOD for each point                                          */
	/************************************************************************/
	vector<IDPair> vectorABOF(PARAM_DATA_N);

	//These vector contains FirstMoment, SecondMoment, Variance corresponding to ABOF
	vector<double> vectorABOF0(PARAM_DATA_N, 0.0);
	vector<double> vectorABOF1(PARAM_DATA_N, 0.0);
	vector<double> vectorABOF2(PARAM_DATA_N, 0.0);

	vector<DVector> vectorEDGES(PARAM_DATA_N, vector<double>(PARAM_DATA_D, 0.0)); // list of AA, AB, AC, ...
	vector<double> vectorLENGTHS(PARAM_DATA_N, 0.0); // list of |AA|, |AB|, |AC|

	vector<double> vectorEDGE(PARAM_DATA_D, 0.0); // vector AB
	vector<double> vectorA; // point A
	vector<double> vectorB; // point B

	double dTemp = 0.0, dLength = 0.0, dABOF = 0.0;
	size_t A, B, iDim;

	for ( A = 0; A < PARAM_DATA_N; A++ )
	{
		vectorA = getPoint(A);

        // Each new point has new vectorEDGES and vectorLENGTHS
		fill(vectorEDGES.begin(), vectorEDGES.end(), vector<double>(PARAM_DATA_D, 0.0));
		fill(vectorLENGTHS.begin(), vectorLENGTHS.end(), 0.0);

		// EDGE stores the coordinates, LENGTH stores the norm_2
		for ( B = 0; B < PARAM_DATA_N; B++ )
		{
		    // If B == A, then everything is zero as initialization
		    if (B == A)
                continue;

            vectorB = getPoint(B);

            // Need to clear everything
            dLength = 0.0;
            fill(vectorEDGE.begin(), vectorEDGE.end(), 0.0);

            for ( iDim = 0; iDim < PARAM_DATA_D; iDim++ )
            {
                dTemp = vectorB[iDim] - vectorA[iDim];
                vectorEDGE[iDim] = dTemp;

                dLength += dTemp * dTemp;
            }

            vectorEDGES[B] = vectorEDGE; // Add vectors AA, AB, AC, ...
            vectorLENGTHS[B] = sqrt(dLength); // Add lengths |AA|, |AB|, |AC|, ...
		}

		// Compute ABOF for each point
		dABOF = computeABOF(vectorEDGES, vectorLENGTHS, vectorABOF1[A], vectorABOF2[A]);

		vectorABOF0[A] = dABOF;
        vectorABOF[A] = make_pair(A, dABOF);

	}

	printf("ABOF time is %f \n", clock() - dStart);

	saveVOAInfo(vectorABOF1, vectorABOF2, vectorABOF0, "ABOF_noSort.txt");

	sort(vectorABOF.begin(), vectorABOF.end(), ascentSortBySec);
    saveRanking(vectorABOF, "ABOF.txt");
}

/**
Ranking data based on approximate ABOF (approxABOF - KDD'08).

Input:
- We need to pre-compute kNN for each point.
**/
void approxABOF()
{
	double dStart = clock();

	vector<IDPair> vectorApproxABOF(PARAM_DATA_N);
	vector<IDPair> vectorkNN;
	for ( size_t A = 0; A < PARAM_DATA_N; A++ )
	{
		// cout << "The point: " << A << endl;
        kNearestNeighbors(vectorkNN, A, PARAM_ABOF_K);
		vectorApproxABOF[A] = make_pair(A, computeApproxABOF(vectorkNN, A));
	}

	printf("Approximate ABOD Time is %f \n", clock() - dStart);

	saveRanking(vectorApproxABOF, "approxABOF_noSort_" + int2str(PARAM_ABOF_K) + ".txt");
    //saveRanking(vectorApproxABOF, "_approxABOF_noSort.txt");

	sort(vectorApproxABOF.begin(), vectorApproxABOF.end(), ascentSortBySec);

    saveRanking(vectorApproxABOF, "approxABOF_" + int2str(PARAM_ABOF_K) + ".txt");
    //saveRanking(vectorApproxABOF, "_approxABOF.txt");
}

