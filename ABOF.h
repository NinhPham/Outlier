#ifndef ABOF_H_INCLUDED
#define ABOF_H_INCLUDED

#include "Header.h"


double computeABOF(const vector<DVector> &, const vector<double> &, double &, double &);
void ABOF();

double computeApproxABOF(const vector<IDPair> & , const size_t);
void approxABOF();

#endif // ABOD_H_VOA
