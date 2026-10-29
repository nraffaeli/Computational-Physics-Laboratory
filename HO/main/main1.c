/*
*******************************************************************************
*
* File main1.c
*
* Harmonic Oscillator.
* This program calculates the correlation function and the autocorrelation in
* Markov time.
*
******************************************************************************/

#define MAIN_PROGRAM
#define MAIN1_C

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
   int Nsweeps;
   double Gamma0;
   double MeanCorrelator[N];
   double *Correlator;
   double value;
   FILE *fCorr;
   FILE *fRead;
   FILE *fAuto;

   rlxd_init(1, SEED);
   ranlxd(xx, N);

   for (j = 0; j < N; j++) xx[j] = 0.0;
   for (i = 0; i < NTHERM; i++) Metropolis();

   Nsweeps = 100000;
   maxLag = Nsweeps / 10;
   if (maxLag > 100) {
      maxLag = 100;
   }

   fCorr = fopen("Correlator.txt", "w");
   if (fCorr == NULL) {
      perror("Error opening Correlator.txt");
      free(Correlator);
      return EXIT_FAILURE;
   }
   for (i = 0; i < Nsweeps; i++) {
      Metropolis();
      for (j = 0; j < N; j++) {
         value = ConfCorrelator(j);
         fprintf(fCorr, "%.15e%s", value, j == N - 1 ? "\n" : " ");
      }
   }
   fclose(fCorr);

   fRead = fopen("Correlator.txt", "r");
   if (fRead == NULL) {
      perror("Error opening Correlator.txt for reading");
      free(Correlator);
      return EXIT_FAILURE;
   }

   fAuto = fopen("Autocorrelator.txt", "w");
   if (fAuto == NULL) {
      perror("Error opening Autocorrelator.txt");
      fclose(fRead);
      free(Correlator);
      return EXIT_FAILURE;
   }

   Correlator = (double *)malloc(Nsweeps * sizeof(double));
   if (Correlator == NULL) {
      perror("Error allocating memory for correlator data");
      return EXIT_FAILURE;
   }

   for (j = 0; j < N; j++) {
      rewind(fRead);
      MeanCorrelator[j] = 0.0;
      for (i = 0; i < Nsweeps; i++) {
         for (k = 0; k < N; k++) {
            if (fscanf(fRead, "%lf", &value) != 1) {
               perror("Error reading Correlator.txt");
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
      MeanCorrelator[j] /= Nsweeps;

      Gamma0 = Gamma(Correlator, MeanCorrelator, Nsweeps, 0, j);
      for (i = 0; i < maxLag; i++) {
         fprintf(fAuto, "%.15e%s", Gamma(Correlator, MeanCorrelator, Nsweeps, i, j) / Gamma0,
              i == maxLag - 1 ? "\n" : " ");
      }
   }

   fclose(fAuto);
   fclose(fRead);
   free(Correlator);

   return EXIT_SUCCESS;
}

