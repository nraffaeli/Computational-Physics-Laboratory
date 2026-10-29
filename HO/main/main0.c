/*
*******************************************************************************
*
* File main0.c
*
* Compute thermalization time for two differnt starting configuration,
* a cold start and a hot start.
*
******************************************************************************/

#define MAIN_PROGRAM
#define MAIN0_C

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "global.h"
#include "random.h"
#include "ho.h"

int main(void)
{
    int i, j; // Markov time, physical time
    FILE *fc;
    FILE *fh;


    rlxd_init(1, SEED);
    ranlxd(xx, N);

    /* Cold start */
    fc = fopen("ColdThermalization.txt", "w");
    if (fc == NULL) {
        printf("Error in opening the file!\n");
        exit(1);
    }
    for (j = 0; j < N; j++) xx[j] = 0.;
    fprintf(fc, "%d\t%lf\n", 0, Action());
    for (i = 0; i < NTHERM; i++) {
        Metropolis();
        fprintf(fc, "%d\t%lf\n", i + 1, Action());
    }

    /* Hot start */
    fh = fopen("HotThermalization.txt", "w");
    if (fh == NULL) {
        printf("Error in opening the file!\n");
        exit(1);
    }
    for (j = 0; j < N; j++) xx[j] = 10.0 * (xx[j] - 0.5);
    fprintf(fh, "%d\t%lf\n", 0, Action());
    for (i = 0; i < NTHERM; i++) {
        Metropolis();
        fprintf(fh, "%d\t%lf\n", i + 1, Action());
    }

    fclose(fc);
    fclose(fh);

    return EXIT_SUCCESS;
}