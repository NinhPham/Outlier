#include "Distance.h"
#include "Utilities.h"

/**
Ranking data based on its kNN distance (Distance-based Outlier - SIGMOD'00).

Input:
- We need to pre-compute kNN for each point.
**/
void kNN()
{
	double dStart = clock();
    vector<IDPair> vectorRank(PARAM_DATA_N);
    double kNN_dist;

	for ( size_t A = 0; A < PARAM_DATA_N; A++ )
	{
		// cout << "The point: " << A << endl;
		kNN_dist = kNN_Dist(A, PARAM_KNN_K);
		vectorRank[A] = make_pair(A, kNN_dist);
	}

	printf("kNN Time is %f \n", clock() - dStart);

    saveRanking(vectorRank, "kNN_noSort_" + int2str(PARAM_KNN_K) + ".txt");

	sort(vectorRank.begin(), vectorRank.end(), descentSortBySec);
    saveRanking(vectorRank, "kNN_" + int2str(PARAM_KNN_K) + ".txt");
}

/**
Ranking data based on its sum of kNN distances (kNN-weight - TKDE 05).

Input:
- We need to pre-compute kNN for each point.
**/
void kNNW()
{
	double dStart = clock();

	vector<IDPair> vectorkNN;
	vector<IDPair> vectorRank(PARAM_DATA_N);
	double kNNW_dist = 0.0;
	size_t A, k;

	for ( A = 0; A < PARAM_DATA_N; A++ )
	{
		// cout << "The point: " << A << endl;
		kNearestNeighbors(vectorkNN, A, PARAM_KNNW_K);

		kNNW_dist = 0.0; // Refresh the distance
		for ( k = 0; k < PARAM_KNNW_K; k++ )
            kNNW_dist += vectorkNN[k].second;

		vectorRank[A] = make_pair(A, kNNW_dist);
	}

	printf("kNNW Time is %f \n", clock() - dStart);

    saveRanking(vectorRank, "kNNW_noSort_" + int2str(PARAM_KNNW_K) + ".txt");

	sort(vectorRank.begin(), vectorRank.end(), descentSortBySec);
    saveRanking(vectorRank, "kNNW_" + int2str(PARAM_KNNW_K) + ".txt");
}

/**
Ranking data based on LOF (SIGMOD'00).

Input:
- We need to pre-compute kNN for each point.
**/
void LOF()
{
	double dStart = clock();

	vector<IDPair> vectorkNN, vector_tempkNN;
	vector<IDPair> vectorLOF_Rank(PARAM_DATA_N);

	vector<vector<IDPair> > vectorkNNs(PARAM_DATA_N); // contains kNN for each points
	vector<double> vectorLRD(PARAM_DATA_N, 0.0); // contain all local reachability density (lrd) for each points

	size_t P, O, k;
	double dSum = 0.0, dDensity = 0.0, dReachDist = 0.0;

	// Get kNN for each point
	for ( P = 0; P < PARAM_DATA_N; P++ )
	{
	    kNearestNeighbors(vectorkNN, P, PARAM_LOF_K);
		vectorkNNs[P] = vectorkNN; // Store kNN vectors for all points to compute the LOF
	}

	// Compute local reachability density of each point
	for ( P = 0; P < PARAM_DATA_N; P++ )
	{
		// cout << "The point: " << P << endl;

		// Compute local reachability density
		dSum = 0.0;
		vectorkNN = vectorkNNs[P];

		for ( k = 0; k < PARAM_LOF_K; k++ )
		{
			// Loop all points O in kNN(P)
			O = vectorkNN[k].first;

			// Get kNN(O)
			vector_tempkNN = vectorkNNs[O];

			// dReachDist = max(kNN_Dist(O), d(P, O))
			dReachDist = max(vector_tempkNN[PARAM_LOF_K - 1].second, vectorkNN[k].second);
			dSum += dReachDist;
		}

		vectorLRD[P] = PARAM_LOF_K / dSum;
	}

	// Compute LOF
	for ( P = 0; P < PARAM_DATA_N; P++ )
	{
		vectorkNN = vectorkNNs[P];
        dDensity = vectorLRD[P];

        dSum = 0.0;
		for ( k = 0; k < PARAM_LOF_K; k++ )
			dSum += vectorLRD[vectorkNN[k].first] / dDensity;

		vectorLOF_Rank[P] = make_pair(P, dSum / PARAM_LOF_K);
	}

	printf("LOF Time is %f \n", clock() - dStart);

    saveRanking(vectorLOF_Rank, "LOF_noSort_" + int2str(PARAM_LOF_K) + ".txt");

	sort(vectorLOF_Rank.begin(), vectorLOF_Rank.end(), descentSortBySec);
    saveRanking(vectorLOF_Rank, "LOF_" + int2str(PARAM_LOF_K) + ".txt");
}


