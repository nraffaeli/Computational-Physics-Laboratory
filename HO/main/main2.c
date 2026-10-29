/*
*******************************************************************************
*
* File main2.c
*
* Harmonic Oscillator.
* This program calculates the correlation function and the autocorrelation in
* Markovian time implementing binning procedure.
*
******************************************************************************/

#define MAIN_PROGRAM
#define MAIN2_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "random.h"
#include "ho.h"

int main(void)
{
	int i, maxLag; // Markov time
	int j, k; // Physical time
	int Nbin, binSize;
	double Gamma0;
	double value;
	double BinnedValues[N];
	double MeanCorrelator[N];
	double *Correlator;
	FILE *fCorr;
	FILE *fRead;
	FILE *fAuto;

	rlxd_init(1, SEED);
	ranlxd(xx, N);

	binSize = 10 * TAU;
	Nbin = NSTAT;

	maxLag = Nbin / 10;
	if (maxLag > 1000) {
		maxLag = 1000;
	}

	for (j = 0; j < N; j++) xx[j] = 0.0;
	for (i = 0; i < NTHERM; i++) Metropolis();

	Correlator = (double *)malloc(Nbin * sizeof(double));
	if (Correlator == NULL) {
		perror("Error allocating memory for correlator data");
		return EXIT_FAILURE;
	}

	fCorr = fopen("BinnedCorrelator.txt", "w");
	if (fCorr == NULL) {
		perror("Error opening BinnedCorrelator.txt");
		free(Correlator);
		return EXIT_FAILURE;
	}
	for (k = 0; k < Nbin; k++) {
		BinnedCorrelators(BinnedValues, binSize);
		for (j = 0; j < N; j++) {
			fprintf(fCorr, "%.15e%s", BinnedValues[j],
					j == N - 1 ? "\n" : " ");
		}
	}
	fclose(fCorr);

	fRead = fopen("BinnedCorrelator.txt", "r");
	if (fRead == NULL) {
		perror("Error opening BinnedCorrelator.txt for reading");
		free(Correlator);
		return EXIT_FAILURE;
	}

	fAuto = fopen("BinnedAutocorrelator.txt", "w");
	if (fAuto == NULL) {
		perror("Error opening BinnedAutocorrelator.txt");
		fclose(fRead);
		free(Correlator);
		return EXIT_FAILURE;
	}
	for (j = 0; j < N; j++) {
		rewind(fRead);
		MeanCorrelator[j] = 0.0;
		for (i = 0; i < Nbin; i++) {
			for (k = 0; k < N; k++) {
				if (fscanf(fRead, "%lf", &value) != 1) {
					perror("Error reading BinnedCorrelator.txt");
					fclose(fAuto);
					fclose(fRead);
					free(Correlator);
					return EXIT_FAILURE;
				}
				if (k == j) {
					Correlator[i] = value;
				}
			}
			MeanCorrelator[j] += Correlator[i];
		}
		MeanCorrelator[j] /= Nbin;

		Gamma0 = Gamma(Correlator, MeanCorrelator, Nbin, 0, j);

		for (i = 0; i < maxLag; i++) {
			fprintf(fAuto, "%.15e%s",
					Gamma(Correlator, MeanCorrelator, Nbin, i, j) / Gamma0,
					i == maxLag - 1 ? "\n" : " ");
		}
	}
	fclose(fAuto);
	fclose(fRead);

	free(Correlator);

	return EXIT_SUCCESS;
}
