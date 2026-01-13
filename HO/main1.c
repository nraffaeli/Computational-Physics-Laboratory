
/*******************************************************************************
*
* File main1.c
*
* Harmonic Oscillator.
* This program simulates a 1D harmonic oscillator using the Metropolis algorithm.
* It calculates the correlation function and the energy of the system.
* The program uses a Markov chain to sample the configuration space and to
* decorrelate the data in Markovian time is implemented a jackknife 
* clusterization.
*
* Author: Nicolò Raffaeli
*
*******************************************************************************/

#define MAIN_PROGRAM

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

   int i, j, k; // tempo Markoviano, tempo fisico, indice di cluster

   for (j = 0; j < N; j++) xx[j] = 5*(xx[j] - 0.5);
   for (i = 0; i < TAUth; i++) Metropolis();


   int Nsweeps = 100000, Dbin = 10*TAUm;
   int Nbin = (int)floor((double)Nsweeps / Dbin);

   double **Corr = (double **)malloc(N * sizeof(double *));
   if (Corr == NULL) {
      perror("Error allocating memory for Corr");
      return EXIT_FAILURE;
   }
   for (j = 0; j < N; j++) {
      Corr[j] = (double *)malloc(Nbin * sizeof(double));
      if (Corr[j] == NULL) {
         perror("Error allocating memory for Corr[j]");
         return EXIT_FAILURE;
      }
   }
   for (j = 0; j < N; j++) {
      for (k = 0; k < Nbin; k++) {
         Corr[j][k] = 0.0;
      }
   }
   
   double MeanBinCorr[N], kBinCorr[N];

   for (k = 0; k < Nbin; k++) {
      
      for (j = 0; j < N; j++) {
         MeanBinCorr[j] = 0.0;
         kBinCorr[j] = 0.0;
      }
      for (i = 0; i < Dbin; i++) {
         Metropolis();
         for (j = 0; j < N; j++) {
            MeanBinCorr[j] += ConfCorrelator(j); 
            if (i == 0) kBinCorr[j] = ConfCorrelator(j);
         }
      }
      for (j = 0; j < N; j++) {
         MeanBinCorr[j] /= Dbin;
         Corr[j][k] = MeanBinCorr[j] - (kBinCorr[j] - MeanBinCorr[j])/(Dbin - 1);
      }
   }

   double CorrBar[N];
   
   for (j = 0; j < N; j++) {
      CorrBar[j] = 0.0;
      for (k = 0; k < Nbin; k++) {
         CorrBar[j] += Corr[j][k];
      }
      CorrBar[j] /= Nbin;
   }

   double **Energy = (double **)malloc(N * sizeof(double *));
   if (Energy == NULL) {
      perror("Error allocating memory for Corr");
      return EXIT_FAILURE;
   }
   for (j = 0; j < N; j++) {
      Energy[j] = (double *)malloc(Nbin * sizeof(double));
      if (Energy[j] == NULL) {
         perror("Error allocating memory for Corr[i]");
         return EXIT_FAILURE;
      }
   }
   for (j = 0; j < N; j++) {
      for (k = 0; k < Nbin; k++) {
         Energy[j][k] = 0.0;
      }
   }

   for (k = 0; k < Nbin; k++) {
      for (j = 0; j < N; j++) { 
         Energy[j][k] = acosh(0.5*(Corr[(j + 1)%N][k] + Corr[(j - 1 + N)%N][k])/Corr[j][k]);
      }
   }

   double EnergyBar[N], varEnergy[N];

   for (j = 0; j < N; j++) {
      varEnergy[j] = 0.0;
      EnergyBar[j] = acosh(0.5*(CorrBar[(j + 1)%N] + CorrBar[(j - 1 + N)%N])/CorrBar[j]);

      for (k = 0; k < Nbin; k++) {
         varEnergy[j] += (Energy[j][k] - EnergyBar[j])*(Energy[j][k] - EnergyBar[j]);
      }
      varEnergy[j] *= (Nbin - 1)/Nbin; 
   }


   for (j = 0; j < N; j++) free(Corr[j]); 
   free(Corr);
  
   for (j = 0; j < N; j++) free(Energy[j]);    
   free(Energy);

   
   FILE *f = fopen("DataHO.txt", "w");
   if(f == NULL) {
      printf("Error opening file!\n");
      exit(1);
   }

   for (j = 0; j < N; j++) {
      fprintf(f, "%d %f %f\n", j, EnergyBar[j], varEnergy[j]);
   }
   
   fclose(f);
   printf("Data saved in DataHO.txt\n");

   exit(0);
}


