#include "VOA.h"
#include "Utilities.h"

/**
Compute the VOA for any point A in KDD'12

Input:
- vectorEDGES: AB, AC, ...
- vectorLENGTHS: |AB|, |AC|, ...

Output:
- VOA = MOA2(Angle) - MOA1^2(Angle)
**/
double computeVOA(const vector<DVector> &vectorEDGES, const vector<double> &vectorLENGTHS, double &MOA1, double &MOA2)
{
	double dCosine = 0.0, dAngle = 0.0;
	double dAngleSum = 0.0, dAngleSum2 = 0.0;

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

			dCosine	= dotProduct(AB, AC) / (vectorLENGTHS[B] * vectorLENGTHS[C]);
			dCosine = truncateCosine(dCosine);

			dAngle  = hackAcos(dCosine);

			dAngleSum	+= dAngle;
			dAngleSum2	+= dAngle * dAngle;
		}
	}

	// VOA
	MOA1 = dAngleSum / PARAM_VOA_C2N;
	MOA2 = dAngleSum2 / PARAM_VOA_C2N;
	return MOA2 - MOA1 * MOA1;
}

/**
Ranking data based on outlier factors VOA that requires O(dn^3).
**/
void VOA()
{
	double dStart = clock();

	/************************************************************************/
	/* Compute VOA for each point                                           */
	/************************************************************************/
	vector<IDPair> vectorVOA(PARAM_DATA_N);

	// These vector contains First Moment, Second Moment and Variance
	vector<double> vectorMOA0(PARAM_DATA_N, 0.0);
	vector<double> vectorMOA1(PARAM_DATA_N, 0.0);
	vector<double> vectorMOA2(PARAM_DATA_N, 0.0);

	vector<DVector> vectorEDGES(PARAM_DATA_N, vector<double>(PARAM_DATA_D, 0.0)); // list of AA, AB, AC, ...
	vector<double> vectorLENGTHS(PARAM_DATA_N, 0.0); // length of |AA|, |AB|, |AC|

	vector<double> vectorEDGE(PARAM_DATA_D, 0.0); // vector AB
	vector<double> vectorA; // point A
	vector<double> vectorB; // point B

	double dTemp, dLength, dVOA;
	size_t A, B, iDim;

	for ( A = 0; A < PARAM_DATA_N; A++ )
	{
		vectorA = getPoint(A);

		// Need to clear everything for each point A
		fill(vectorEDGES.begin(), vectorEDGES.end(), vector<double>(PARAM_DATA_D, 0.0));
		fill(vectorLENGTHS.begin(), vectorLENGTHS.end(), 0.0);

		// EDGE stores the coordinates, LENGTH stores the norm_2
		for ( B = 0; B < PARAM_DATA_N; B++ )
		{
		    // IF B == A, everything is zero as initialization
		    if (B == A)
                continue;

			vectorB = getPoint(B);

            dLength = 0.0;
			fill(vectorEDGE.begin(), vectorEDGE.end(), 0.0);

			for ( iDim = 0; iDim < PARAM_DATA_D; iDim++ )
			{
				dTemp = vectorB[iDim] - vectorA[iDim];
				vectorEDGE[iDim] = dTemp; // Coordinates of vector AB

				dLength += dTemp * dTemp;
			}

			vectorEDGES[B] = vectorEDGE; // Add vectors AB, AC, ...
			vectorLENGTHS[B] = sqrt(dLength); // Add lengths |AB|, |AC|, ...
		}

		// Compute scores VOA, MOA1, MOA2 for each point
		dVOA = computeVOA(vectorEDGES, vectorLENGTHS, vectorMOA1[A], vectorMOA2[A]);

		vectorVOA[A] = make_pair(A, dVOA);
		vectorMOA0[A] = dVOA;
	}

	printf("VOA time is %f \n", clock() - dStart);

	saveVOAInfo(vectorMOA1, vectorMOA2, vectorMOA0, "VOA_noSort.txt");

	sort(vectorVOA.begin(), vectorVOA.end(), ascentSortBySec);
	saveRanking(vectorVOA, "VOA.txt");
}


