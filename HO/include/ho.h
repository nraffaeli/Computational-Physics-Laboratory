#ifndef HO_H
#define HO_H

#include "global.h"

extern double Action(void);

extern double DeltaAction(double varx, int n);

extern void Metropolis(void);

extern double ConfCorrelator(int delta);

extern void BinnedCorrelators(double Correlators[], int binSize);

extern double BinnedCorrelator(int delta, int binSize);

extern double Gamma(const double Correlator[], double MeanCorrelator[], int n, int t, int j);

extern double Energy(const double Correlator[], int j);

extern double MatrixElement(const double Correlator[], double Energy, int j);

extern void cluster_jackknife(double data[], int Ncluster);

extern double variance_jackknife(double data[], int Ncluster, double mean);

#endif
