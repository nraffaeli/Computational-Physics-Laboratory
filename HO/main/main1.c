/*
*******************************************************************************
*
* File main1.c
*
* Harmonic Oscillator.
* This program calculates the correlation function and the autocorrelation in 
* Markovian time. 
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

   rlxd_init(1, 12345);
   ranlxd(xx, N);

   int i, j, k; // Markovian time, physical time, autocorrelation lag

   for (j = 0; j < N; j++) xx[j] = 0.0;
   for (i = 0; i < TAUth; i++) Metropolis();


   int Nsweeps = 100000;

   double *MarkovCorr = (double *)malloc(Nsweeps * sizeof(double));
   if (MarkovCorr == NULL) {
      perror("Error allocating memory for MarkovCorr");
      return EXIT_FAILURE;
   }

   double CorrBar[N] = {0.0};

   for (i = 0; i < Nsweeps; i++) {
      Metropolis();
      MarkovCorr[i] = ConfCorrelator(1);
      for (j = 0; j < N; j++) {
         CorrBar[j] += ConfCorrelator(j);
      }
   }
   for (j = 0; j < N; j++) {
      CorrBar[j] /= Nsweeps;
   }

   FILE *fCorr = fopen("Correlator.txt", "w");
   if (fCorr == NULL) {
      perror("Error opening Correlator.txt");
      free(MarkovCorr);
      return EXIT_FAILURE;
   }
   for (j = 0; j < N; j++) {
      fprintf(fCorr, "%d %f\n", j, CorrBar[j]);
   }
   fclose(fCorr);

   double MarkovMean = 0.0;
   for (i = 0; i < Nsweeps; i++) {
      MarkovMean += MarkovCorr[i];
   }
   MarkovMean /= Nsweeps;

   double Gamma0 = 0.0;
   for (i = 0; i < Nsweeps; i++) {
      double fluctuation = MarkovCorr[i] - MarkovMean;
      Gamma0 += fluctuation * fluctuation;
   }
   Gamma0 /= Nsweeps;

   FILE *fAuto = fopen("Autocorrelator.txt", "w");
   if (fAuto == NULL) {
      perror("Error opening Autocorrelator.txt");
      free(MarkovCorr);
      return EXIT_FAILURE;
   }
   if (Gamma0 == 0.0) {
      fprintf(stderr, "Cannot calculate autocorrelation: zero variance.\n");
   } else {
      int maxLag = 1000;
      for (k = 0; k < maxLag; k++) {
         double gamma = 0.0;
         for (i = 0; i < Nsweeps - k; i++) {
            gamma += (MarkovCorr[i] - MarkovMean) *
                     (MarkovCorr[i + k] - MarkovMean);
         }
         gamma /= Nsweeps - k;
         fprintf(fAuto, "%d %f\n", k, gamma / Gamma0);
      }
   }
   fclose(fAuto);

   free(MarkovCorr);

   exit(0);
}

