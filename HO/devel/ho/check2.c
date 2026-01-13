
/*******************************************************************************
*
* File check2.c
*
* Check for Mentropolis algorithm and Correlation behavior over sweeps
*
* Author: Nicolò Raffaeli
*
*******************************************************************************/
#define MAIN_PROGRAM
#define CHECK2_C

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "global.h"
#include "random.h"
#include "ho.h"

int main(void)
{
    int i, Nsweeps;

    rlxd_init(1, 12345);
    ranlxd(xx, N);


    printf("Thermalization time for Harmonic Oscillator:\n");

    FILE *file1 = fopen("ThermalizationTime.txt", "w");
    if (file1 == NULL) {
        printf("Error in opening the file!\n");
        exit(1);
    }

    //for (i = 0; i < N; i++) xx[i] = 0.;  --> Cold Start
    for (i = 0; i < N; i++) xx[i] = 5*(xx[i] - 0.5);

    Nsweeps = 100;
    for (i = 0; i < Nsweeps; i++) {
        Metropolis();
        fprintf(file1, "%d\t%lf\n", i, Action());
    }

    fclose(file1);
    printf("Data saved in ThermalizationTime.txt\n\n");


    printf("Correlation for Harmonic Oscillator:\n");
    
    int j, k;
    Nsweeps = 1000;

    double **Matrix = (double **)malloc(N * sizeof(double *));
    if (Matrix == NULL) {
        perror("Error allocating memory for Matrix");
        return EXIT_FAILURE;
    }

    for (i = 0; i < N; i++) {
        Matrix[i] = (double *)malloc(Nsweeps * sizeof(double));
        if (Matrix[i] == NULL) {
            perror("Error allocating memory for Matrix[i]");
            return EXIT_FAILURE;
        }
    }

    for (i = 0; i < Nsweeps; i++) {
        Metropolis();

        for ( j = 0; j < N; j++) {
            Matrix[j][i] = ConfCorrelator(j);
        }
    }

    double Corr[N] = {0}, varCorr[N] = {0};

    for (j = 0; j < N; j++) {
        for ( i = 0; i < Nsweeps; i++) {
            Corr[j] += Matrix[j][i];
        }
        Corr[j] /= Nsweeps;
    }

    for(j = 0; j < N; j++) {
        for (i = 0; i < Nsweeps; i++) {
            varCorr[j] += (Matrix[j][i] - Corr[j])*(Matrix[j][i] - Corr[j]);
        }
        varCorr[j] /= Nsweeps - 1;
    }

    double *Gamma = (double *)calloc(20, sizeof(double));
    if (Gamma == NULL) {
        perror("Error allocating memory for Gamma");
        return EXIT_FAILURE;
    }

    for (k = 0; k < 20; k++) {
        Gamma[k] = 0;
        for (i = 0; i < Nsweeps - k; i++) {
            Gamma[k] += (Matrix[1][i] - Corr[1])*(Matrix[1][(i+k)] - Corr[1]);
        }
        Gamma[k] /= (Nsweeps - k)*varCorr[1];
    }

    FILE *file2 = fopen("Correlation.txt", "w");
    if (file2 == NULL) {
        printf("Error in opening the file!\n");
        exit(1);
    } 

    for (j = 0; j < N; j++) {
        fprintf(file2, "%d\t%lf\n", j, Corr[j]);
    }
    
    fclose(file2);
    printf("Data saved in Correlation.txt\n");

    FILE * file3 = fopen("Gamma.txt", "w");
    if (file3 == NULL) {
        printf("Error in opening the file!\n");
        exit(1);
    }

    for (k = 0; k < 20; k++) {
        fprintf(file3, "%d\t%lf\n", k, Gamma[k]);
    }

    fclose(file3);
    printf("Data saved in Gamma.txt\n");


    for (i = 0; i < N; i++) {
    free(Matrix[i]);
    exit(0);
    }
    free(Matrix);
    free(Gamma);

   exit(0);
}

