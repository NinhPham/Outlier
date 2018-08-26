#include "RP_VOA.h"
#include "Utilities.h"

/**
Compute AMS Sketch-based VOA and then rank points (FastVOA - KDD'12)

Algorithm:
- F1 is computed as the equation (1)
- F2 is computed as the equation (2)' in KDD'12 using AMS Sketches
**/
void FastVOA()
{
	double dStart = clock();

    vector<double> vectorNormal(PARAM_VOA_RP_R * PARAM_DATA_D); // vector of size R * N
	generateNormal(vectorNormal);

	size_t r = 0, n = 0, s1 = 0, s2 = 0;
	size_t iBaseIndex = 0, iIdx1 = 0, iIdx2 = 0;

	/************************************************************************/
	/* Sorting based on dot product                                                              */
	/************************************************************************/
	vector<IDPair> vectorOrders(PARAM_DATA_N * PARAM_VOA_RP_R);
	orderStats(vectorOrders, vectorNormal);

	/************************************************************************/
	/* Only for testing the accuracy - Not used when measure running time   */
	/************************************************************************/
	//NaiveRP(vectorOrders);

	/************************************************************************/
	/* Compute first moment                                                         */
	/************************************************************************/
    vector<double> FIRST_MOMENT(PARAM_DATA_N, 0.0);

    for (r = 0; r < PARAM_VOA_RP_R; r++ )
	{
	    // The random projection *iRan*
	    iBaseIndex = r * PARAM_DATA_N;

		for ( n = 0; n < PARAM_DATA_N; n++ )
		{
			iIdx1			        = vectorOrders[iBaseIndex + n].first; // The position n is the rank of the point iIdx1
			FIRST_MOMENT[iIdx1]     += n * (PARAM_DATA_N - 1 - n);
		}
	}

	/************************************************************************/
	/* Compute second moment using AMS Sketch                                                      */
	/************************************************************************/

	vector<int> AMS_LEFT(PARAM_DATA_N, 0), AMS_RIGHT(PARAM_DATA_N, 0);
	vector<int> AMS(PARAM_DATA_N, 0);

	vector<DVector> AMS_S2(PARAM_AMS_S2, vector<double> (PARAM_DATA_N, 0.0));

    vector<int> FOUR_WISE_LEFT(PARAM_DATA_N, 0);
    vector<int> FOUR_WISE_RIGHT(PARAM_DATA_N, 0);

	/************************************************************************/
	/* S2 Processing                                                        */
	/************************************************************************/

	for ( s2 = 0; s2 < PARAM_AMS_S2; s2++ )
	{

		/************************************************************************/
		/* S1 Processing                                                        */
		/************************************************************************/

		for ( s1 = 0; s1 < PARAM_AMS_S1; s1++ )
		{
			fill(AMS.begin(), AMS.end(), 0);

			// No need to be cleared since it will assign new random values
			generate4Wises(FOUR_WISE_LEFT, FOUR_WISE_RIGHT);

			/************************************************************************/
			/* NUM_RAN Processing                                                   */
			/************************************************************************/

			for ( r = 0; r < PARAM_VOA_RP_R; r++ )
			{
				fill(AMS_LEFT.begin(), AMS_LEFT.end(), 0);
				fill(AMS_RIGHT.begin(), AMS_RIGHT.end(), 0);

				iBaseIndex = r * PARAM_DATA_N;

				/************************************************************************/
				/* Left part                                                            */
				/************************************************************************/

				for ( n = 1; n < PARAM_DATA_N; n++ )
				{
					iIdx1 = vectorOrders[iBaseIndex + n].first;
                    iIdx2 = vectorOrders[iBaseIndex + n - 1].first;

					AMS_LEFT[iIdx1] = AMS_LEFT[iIdx2] + FOUR_WISE_LEFT[iIdx2];
				}

				/************************************************************************/
				/* Right part                                                           */
				/************************************************************************/

				for ( n = PARAM_DATA_N - 1; n >= 1; n-- ) // Note here since n-- and size_t n, we might get 0 then 2^31 -1
				{
					iIdx1 = vectorOrders[iBaseIndex + n - 1].first;
					iIdx2 = vectorOrders[iBaseIndex + n].first;

					AMS_RIGHT[iIdx1] = AMS_RIGHT[iIdx2] + FOUR_WISE_RIGHT[iIdx2];
				}

				/***************************************************************/
				/* SUM OF ( AMS_LEFT * AMS_RIGHT )                             */
				/***************************************************************/

				for ( n = 0; n < PARAM_DATA_N; n++ )
					AMS[n] += AMS_LEFT[n] * AMS_RIGHT[n];

			}

			/************************************************************************/
			/* SUM OF ( SECOND MOMENT / NUM_S1 ) for each S1                        */
			/************************************************************************/

			for ( n = 0; n < PARAM_DATA_N; n++ )
			{
				AMS_S2[s2][n] += pow((double)AMS[n], 2) / PARAM_AMS_S1;
			}
		}
	}


	/************************************************************************/
	/* Sort and Get Median                                                  */
	/************************************************************************/

	size_t iMedian = PARAM_AMS_S2 / 2;
	vector<double> vectorS2(PARAM_AMS_S2);
	vector<double> SECOND_MOMENT(PARAM_DATA_N, 0.0);

	for ( n = 0; n < PARAM_DATA_N; n++ )
	{
		//fill(vectorS2.begin(), vectorS2.end(), 0.0);

		for ( s2 = 0; s2 < PARAM_AMS_S2; s2++ )
			vectorS2[s2] = AMS_S2[s2][n];

		//sort(vectorS2.begin(), vectorS2.end());
		nth_element(vectorS2.begin(), vectorS2.begin() + iMedian, vectorS2.end()); // faster than sort
		SECOND_MOMENT[n] = vectorS2[iMedian];
	}


	/************************************************************************/
	/* Variance                                                             */
	/************************************************************************/

	vector<IDPair> vectorRank(PARAM_DATA_N);
	vector<double> VOA(PARAM_DATA_N, 0.0);

	for ( n = 0; n < PARAM_DATA_N; n++ )
	{
		FIRST_MOMENT[n]  = PI * FIRST_MOMENT[n] / (PARAM_VOA_C2N * PARAM_VOA_RP_R);
		SECOND_MOMENT[n] = PI * PI * SECOND_MOMENT[n] / (PARAM_VOA_C2N * PARAM_VOA_C2R);
		SECOND_MOMENT[n] = SECOND_MOMENT[n] - 2 * PI * FIRST_MOMENT[n] / (PARAM_VOA_RP_R - 1);


		VOA[n] = SECOND_MOMENT[n] - pow(FIRST_MOMENT[n], 2);
		vectorRank[n] = make_pair(n, VOA[n]);
	}

	printf("AMS-based VOA time is %f \n", clock() - dStart);


    saveVOAInfo(FIRST_MOMENT, SECOND_MOMENT, VOA, "FastVOA_noSort_" + int2str(PARAM_VOA_RP_R) + "_" + int2str(PARAM_AMS_S1) + "_" + int2str(PARAM_AMS_S2) + ".txt");

	sort(vectorRank.begin(), vectorRank.end(), ascentSortBySec);

	saveRanking(vectorRank, "FastVOA_" + int2str(PARAM_VOA_RP_R) + "_" + int2str(PARAM_AMS_S1) + "_" + int2str(PARAM_AMS_S2) + ".txt");
}

/**
Naive Random projection based F1, F2, and derived VAR for all points requires O(n^3)
This function is used for testing the accuracy of Random Project-based techniques

Algorithm:
- The first moment is computed as the equation (1) in KDD'12
- The second moment is computed as the equation (2) in KDD'12
**/
void NaiveRP(const vector<IDPair> &vectorOrders)
{
	double dStart = clock();

	size_t r = 0, n = 0, iIndex = 0;
	size_t iBaseIndex = 0;
	size_t iIdx1 = 0, iIdx2 = 0;
	size_t A, B;

	/************************************************************************/
	/* Sorting - We do not need it since vectorRank will be the input                                                              */
	/************************************************************************/
	//vector<double> vectorNormal(PARAM_VOA_RP_R * PARAM_DATA_D); // vector of size R * N
	//generateNormal(vectorNormal);
	//vector<IDPair> vectorOrders(PARAM_VOA_RP_R * PARAM_DATA_N);
	//orderStats(vectorOrders, vectorNormal);

	/************************************************************************/
	/* First moment                                                         */
	/************************************************************************/
	vector<double> FIRST_MOMENT(PARAM_DATA_N, 0.0);

	for ( r = 0; r < PARAM_VOA_RP_R; r++ )
	{
		iBaseIndex = r * PARAM_DATA_N;

		for ( n = 0; n < PARAM_DATA_N; n++ )
		{
			iIdx1			    = vectorOrders[iBaseIndex + n].first; // n is the rank of the point iIdx1
			FIRST_MOMENT[iIdx1] += n * (PARAM_DATA_N - 1 - n);
		}
	}

	/************************************************************************/
	/* Second moment                                                        */
	/************************************************************************/

	vector<double> SECOND_MOMENT(PARAM_DATA_N, 0.0);
	vector<DVector> MATRIX_P(PARAM_DATA_N, vector<double>(PARAM_DATA_N, 0.0));

	for ( n = 0; n < PARAM_DATA_N; n++ )
	{
		fill(MATRIX_P.begin(), MATRIX_P.end(), vector<double> (PARAM_DATA_N, 0.0));

		for ( r = 0; r < PARAM_VOA_RP_R; r++ )
		{
			iBaseIndex = r * PARAM_DATA_N;

			for ( iIdx1 = 0; iIdx1 < PARAM_DATA_N; iIdx1++ )
			{
				if ( vectorOrders[iBaseIndex + iIdx1].first == n)
                {
                    iIndex = iIdx1; // find the index of the point
                    break;
                }
			}

			for ( iIdx1 = 0; iIdx1 < iIndex; iIdx1++ ) // from 0 to the (index of the point - 1)
			{
			    A = vectorOrders[iBaseIndex + iIdx1].first;

				for ( iIdx2 = iIndex + 1; iIdx2 < PARAM_DATA_N; iIdx2++ ) // from the (index of point + 1) to (NUM_POINT - 1)
                {
                    B = vectorOrders[iBaseIndex + iIdx2].first;
                    MATRIX_P[A][B] = MATRIX_P[A][B] + 1.0;
                }
			}
		}

		// Compute squared Frobenius norm of P in Equation (2) - KDD 12
		for ( iIdx1 = 0; iIdx1 < PARAM_DATA_N; iIdx1++ )
		{
			for ( iIdx2 = 0; iIdx2 < PARAM_DATA_N; iIdx2++ )
				SECOND_MOMENT[n] += pow(MATRIX_P[iIdx1][iIdx2], 2);
		}
	}

	/************************************************************************/
	/* The First and Second Moments based on Random Projection              */
	/************************************************************************/
    vector<double> VOA(PARAM_DATA_N, 0.0);
	for ( n = 0; n < PARAM_DATA_N; n++ )
	{
		FIRST_MOMENT[n] = PI * FIRST_MOMENT[n] / ( PARAM_VOA_C2N * PARAM_VOA_RP_R );

        SECOND_MOMENT[n] = PI * PI * SECOND_MOMENT[n] / (PARAM_VOA_C2N * PARAM_VOA_C2R);
		SECOND_MOMENT[n] = SECOND_MOMENT[n] - 2 * PI * FIRST_MOMENT[n] / (PARAM_VOA_RP_R - 1);

        VOA[n] = SECOND_MOMENT[n] - pow(FIRST_MOMENT[n], 2);
	}

	printf("Naive_RP_F1_F2 time is %f \n", clock() - dStart);

	saveVOAInfo(FIRST_MOMENT, SECOND_MOMENT, VOA, "NaiveRP_noSort.txt");
}

/**
Compute dot product between each point and random vectors
Sort point ID based on their product
**/
void orderStats(vector<IDPair> &vectorOrders, const vector<double> &vectorNormal)
{
    vector<double> vectorA, vectorR;
    vector<IDPair> vectorInner(PARAM_DATA_N);
    size_t r, A;

	for ( r = 0; r < PARAM_VOA_RP_R; r++ )
	{
	    vectorR = vector<double>(vectorNormal.begin() + r * PARAM_DATA_D, vectorNormal.begin() + (r + 1) * PARAM_DATA_D);
		fill(vectorInner.begin(), vectorInner.end(), make_pair(0, 0.0));

		for ( A = 0; A < PARAM_DATA_N; A++)
		{
			vectorA = getPoint(A);
			vectorInner[A] = make_pair(A, dotProduct(vectorR, vectorA));
		}

		sort(vectorInner.begin(), vectorInner.end(), ascentSortBySec);
		copy(vectorInner.begin(), vectorInner.end(), vectorOrders.begin() + r * PARAM_DATA_N);
	}
}
