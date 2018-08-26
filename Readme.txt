The C++ source code (Code::Block IDE) for outlier detection in high dimensions.
It covers standard methods and its approximations to detect outliers in high-dimensional data sets, including
- kNN, kNNW, Sam1NN
- LOF
- ABOF, approxABOF
- VOA, FastVOA
- L1Depth, SamDepth


Parameters: <NUM_POINT> <NUM_DIM> <FILE_NAME> <METHOD> <ADDITIONAL_PARAMS>

- NUM_POINT  : number of points (N)
- NUM_DIM  : number of dimensions (D)
- FILE_NAME     : filename of dataset with the matrix format N x D

- <METHOD> with <ADDITIONAL_PARAMS>

"kNN": SIGMOD 00 - Efficient Algorithms for Mining Outliers from Large Data Sets
    - k (recommendation: 10)
"Sam1NN": NIPS 13 - Rapid Distance-Based Outlier Detection via Sampling
    - number of sample points (recommendation: 20)

"kNNW": TKDE 05 - Outlier Mining in Large High-Dimensional Data Sets
    - k (recommendation: 10)
"LOF": SIGMOD 00 - LOF: Identifying Density-Based Local Outliers
    - minPts (recommendation: 40)

"ABOF": KDD 08 - Angle-Based Outlier Detection in High-dimensional Data
"approxABOF": KDD 08 - Angle-Based Outlier Detection in High-dimensional Data
    - k (recommendation: 0.1 * N)

"VOA": KDD 12 - A near-linear time approximation algorithm for angle-based outlier detection in high-dimensional data
"FastVOA": KDD 12 - A near-linear time approximation algorithm for angle-based outlier detection in high-dimensional data
    - number of random projections (recommendation: 100)
    - AMS Sketch size S1 (recommendation: 3200)
    - AMS Sketch size S2 (recommendation: 5)

"L1D": PKDD 18 - L1-Depth Revisited - A Robust Angle-based Outlier Factor in High-dimensional Space
"BasicSamL1D": PKDD 18 - L1-Depth Revisited - A Robust Angle-based Outlier Factor in High-dimensional Space
    - number of sample pairs (recommendation: sqrt(N))
"SamDepth": PKDD 18 - L1-Depth Revisited - A Robust Angle-based Outlier Factor in High-dimensional Space
    - number of sample points (recommendation: sqrt(N))


Example: 
- 60839 41 "C:\_Data\Dataset\_L1D\Datasets\KDDCup99_norm_idf_60839_41_246.txt" "Sam1NN" 20
- 60839 41 "C:\_Data\Dataset\_L1D\Datasets\KDDCup99_norm_idf_60839_41_246.txt" "FastVOA" 100 3200 5

