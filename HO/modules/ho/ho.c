/*******************************************************************************
*
* File ho.c
*
* Harmonic Oscillator modules implemented
*
*******************************************************************************/
#define HO_C

#include <math.h>
#include <stdlib.h>
#include "ho.h"
#include "global.h"
#include "random.h"

double Action(void) {
    int i;
    double s = .0;

    for (i = 0; i < N; i++) {
        s += 0.5*M*(xx[(i+1)%N]-xx[i])*(xx[(i+1)%N]-xx[i]) + 0.5*M*W*W*xx[i]*xx[i];
    }

    return s;
}

double DeltaAction(double varx, int n) {
    return - M*(xx[(n+1)%N] + xx[(n+N-1)%N])*(varx - xx[n]) + M*(1 + W*W/2.)*(varx*varx - xx[n]*xx[n]);
}

void Metropolis(void) {
    int i;
    double dx = .0, dS = .0, r[2*N];
    ranlxd(r, 2*N);

    for (i = 0; i < N; i++) {
        dx = 2*DELTA*(r[i] - 0.5);
        dS = DeltaAction(xx[i] + dx, i);

        if (dS <= 0.0 || r[i + N] < exp(-dS)) {
            xx[i] += dx;
        }
    }
}

double ConfCorrelator(int delta) {
    int i;
    double res = .0;

    for (i = 0; i < N; i++) {
        res += xx[i]*xx[(i+delta)%N];
    }

    return res/N;
}

void BinnedCorrelators(double Correlators[], int binSize) {
    int i, delta;

    for (delta = 0; delta < N; delta++) {
        Correlators[delta] = .0;
    }

    for (i = 0; i < binSize; i++) {
        Metropolis();
        for (delta = 0; delta < N; delta++) {
            Correlators[delta] += ConfCorrelator(delta);
        }
    }

    for (delta = 0; delta < N; delta++) {
        Correlators[delta] /= binSize;
    }
}

double BinnedCorrelator(int delta, int binSize) {
    double correlators[N];

    BinnedCorrelators(correlators, binSize);
    return correlators[delta];
}

double Gamma(const double Correlator[], double MeanCorrelator[], int n, int t, int j) {
    int i;
    double res = 0.0;

    for (i = 0; i < n - t; i++) {
        res += Correlator[i] * Correlator[i + t];
    }
    return (res)/(n-t) - MeanCorrelator[j] * MeanCorrelator[j];
}

double Energy(const double Correlator[], int j) {
    double res = Correlator[(j + 1)%N] + Correlator[(j - 1 + N)%N];
    return acosh(res / (2 * Correlator[j]));
}

double MatrixElement(const double Correlator[], double Energy, int j) {
    double denominator = 2 * cosh((N/2 - j) * Energy);
    return Correlator[j] * exp(N/2 * Energy) / denominator;
}

void cluster_jackknife(double data[], int Ncluster)
{
    int i, t, cluster;
    double *totals;

    if (Ncluster < 2) {
        return;
    }

    totals = (double *)calloc(N, sizeof(double));
    if (totals == NULL) {
        return;
    }

    for (cluster = 0; cluster < Ncluster; cluster++) {
        for (t = 0; t < N; t++) {
            totals[t] += data[cluster * N + t];
        }
    }

    for (i = 0; i < Ncluster; i++) {
        for (t = 0; t < N; t++) {
            data[i * N + t] =
                (totals[t] - data[i * N + t]) / (Ncluster - 1);
        }
    }

    free(totals);
}

double variance_jackknife(double data[], int Ncluster, double mean)
{
    int j;
    double variance = 0.0;

    for (j = 0; j < Ncluster; j++) {
        variance += (data[j] - mean) * (data[j] - mean);
    }

    return variance * (Ncluster - 1) / Ncluster;
}