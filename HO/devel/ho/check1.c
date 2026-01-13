
/*******************************************************************************
*
* File check1.c
*
* Check for the action of AO (A Try)
*
* Author: Nicolò Raffaeli
*
*******************************************************************************/
#define MAIN_PROGRAM
#define CHECK1_C

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "global.h"
#include "random.h"
#include "ho.h"

int main(void)
{

   int i, n = 10;
   double vtry[] = {1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.4, 1.3, 1.2, 1.1, 1};

   for (i = 0; i < N; i++) xx[i] = 0.;
   
   printf("Print of the vector of %d element:\n", n);

   for (i = 0; i < n+1; i++) {
      xx[i] = vtry[i];
      printf("%lf \n", xx[i]);
   }
   printf("\n");
   printf("Action value to check:\n-->\t%lf\n\n", Action() - M);


   rlxd_init(1, 12345);
   ranlxd(xx, N);

   double Rif = DeltaAction(xx[n]+0.5, n), S1 = Action(), S2;
   
   xx[n] = xx[n] + 0.5;
   S2 = Action();

   printf("Check for DeltaAction:\nIf the value is close to zero then the routine is correct.\n-->\t%lf\n\n", S2 - S1 - Rif);

   exit(0);
}


