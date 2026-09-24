/*
*******************************************************************************
*
* File main1.c
*
* In this code we check the termalization time for two differnt starting 
* configurations, a cold start and a hot start. 
*
******************************************************************************/

#define MAIN_PROGRAM
#define MAIN0_C

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

    FILE *fc = fopen("ColdThermalization.txt", "w");
    if (fc == NULL) {
        printf("Error in opening the file!\n");
        exit(1);
    }

    FILE *fh = fopen("HotThermalization.txt", "w");
    if (fh == NULL) {
        printf("Error in opening the file!\n");
        exit(1);
    }

    int i, j; /* Markov time, physical time */

    /* Cold start */
    for (j = 0; j < N; j++) xx[j] = 0.;
    fprintf(fc, "%d\t%lf\n", 0, Action());
    for (i = 0; i < TAUth; i++) {
        Metropolis();
        fprintf(fc, "%d\t%lf\n", i + 1, Action());
    }

    /* Hot start */
    for (j = 0; j < N; j++) xx[j] = 10.0 * (xx[j] - 0.5);
    fprintf(fh, "%d\t%lf\n", 0, Action());
    for (i = 0; i < TAUth; i++) {
        Metropolis();
        fprintf(fh, "%d\t%lf\n", i + 1, Action());
    }

    fclose(fc);
    fclose(fh);

    return 0;
}