/*
*******************************************************************************
*
* File main3.c
*
* Given the correlation function, we perform a jackknife analysis to estimate
* Energy and its error, togheter with the matrix element.
*
******************************************************************************/

#define MAIN_PROGRAM
#define MAIN3_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "random.h"
#include "ho.h"


int main(void)
{
    int j, k;
    int Nbin;
    double threshold = 0.2; // Threshold for trusting Energy and Matrix Element results
    double value;
    double MeanEnergy, MeanME, variance_energy, variance_me;
    double *bins;
    double *MeanCorrelator;
    double *JackEnergy;
    double *JackME;
    FILE *fRead;
    FILE *fEnergy;
    FILE *fME;

    Nbin = NSTAT;

    bins = (double *)malloc((size_t)Nbin * N * sizeof(double));
    MeanCorrelator = (double *)calloc(N, sizeof(double));
    JackEnergy = (double *)malloc((size_t)Nbin * sizeof(double));
    if (bins == NULL || MeanCorrelator == NULL || JackEnergy == NULL) {
        perror("Error allocating correlator data");
        free(bins);
        free(MeanCorrelator);
        free(JackEnergy);
        return EXIT_FAILURE;
    }
    JackME = (double *)malloc((size_t)Nbin * sizeof(double));
    if (JackME == NULL) {
        perror("Error allocating matrix element data");
        free(bins);
        free(MeanCorrelator);
        free(JackEnergy);
        return EXIT_FAILURE;
    }

    fRead = fopen("BinnedCorrelator.txt", "r");
	if (fRead == NULL) {
		perror("Error opening BinnedCorrelator.txt for reading");
		free(bins);
		free(MeanCorrelator);
        free(JackEnergy);
		return EXIT_FAILURE;
	}

    for (k = 0; k < Nbin; k++) {
        for (j = 0; j < N; j++) {
            if (fscanf(fRead, "%lf", &value) != 1) {
                perror("Error reading BinnedCorrelator.txt");
                fclose(fRead);
                free(bins);
                free(MeanCorrelator);
                free(JackEnergy);
                return EXIT_FAILURE;
            }
            bins[k * N + j] = value;
            MeanCorrelator[j] += value;
        }
    }
    for (j = 0; j < N; j++) {
        MeanCorrelator[j] = MeanCorrelator[j] / Nbin;
    }

    /* Run Clustering */
    cluster_jackknife(bins, Nbin);

    fEnergy = fopen("Energy.txt", "w");
    if (fEnergy == NULL) {
        perror("Error opening Energy.txt for writing");
        fclose(fRead);
        free(bins);
        free(MeanCorrelator);
        free(JackEnergy);
        return EXIT_FAILURE;
    }
    fME = fopen("MatrixElement.txt", "w");
    if (fME == NULL) {
        perror("Error opening MatrixElement.txt for writing");
        fclose(fRead);
        fclose(fEnergy);
        free(bins);
        free(MeanCorrelator);
        free(JackEnergy);
        return EXIT_FAILURE;
    }

    for (j = 2; j < N/2 - 1; j++) {
        MeanEnergy = Energy(MeanCorrelator, j);
        MeanME = MatrixElement(MeanCorrelator, MeanEnergy, j);

        for (k = 0; k < Nbin; k++) {
            JackEnergy[k] = Energy(&bins[k * N], j);
            JackME[k] = MatrixElement(&bins[k * N], JackEnergy[k], j);
        }
        variance_energy = variance_jackknife(JackEnergy, Nbin,MeanEnergy);
        fprintf(fEnergy, "%d\t%.15e\t%.15e\n", j, MeanEnergy, sqrt(variance_energy));

        variance_me = variance_jackknife(JackME, Nbin, MeanME);
        fprintf(fME, "%d\t%.15e\t%.15e\n", j, MeanME, sqrt(variance_me));

        if (sqrt(variance_energy) / fabs(MeanEnergy) > threshold) {
            printf("Warning: Energy error exceeds threshold at j = %d\n", j);
        }
    }

    free(bins);
    free(MeanCorrelator);
    free(JackEnergy);
    fclose(fRead);
    fclose(fEnergy);
    fclose(fME);

    return EXIT_SUCCESS;
}
