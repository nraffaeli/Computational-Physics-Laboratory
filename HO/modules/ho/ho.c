/*******************************************************************************
*
* File ho.c
* 
* Harmonic Oscillator modules implemented
*
*******************************************************************************/
#define HO_C

#include <math.h>
#include "ho.h"
#include "global.h"
#include "random.h"

double Action(void)
{  
    int i;
    double s = .0;

    for (i = 0; i < N; i++)
    {
        s += 0.5*M*(xx[(i+1)%N]-xx[i])*(xx[(i+1)%N]-xx[i]) + 0.5*M*W*W*xx[i]*xx[i];
    }
            
    return s;
}

double DeltaAction(double varx, int n)
{        
    return - M*(xx[n+1] + xx[n-1])*(varx - xx[n]) + M*(1 + W*W/2.)*(varx*varx - xx[n]*xx[n]);
}

void Metropolis(void)
{
    int i;
    double dx = .0, dS = .0, r[2*N]; 
    ranlxd(r, 2*N);
    
    for (i = 0; i < N; i++)
    {
        dx = 2*DELTA*(r[i] - 0.5);
        dS = DeltaAction(xx[i] + dx, i);
        
        if (r[i+N] < exp(-dS))
        {
            xx[i] += dx;
        }
    }
}

double ConfCorrelator(int delta)
{
    int i;
    double res = .0;
    
    for (i = 0; i < N; i++)
    {
        res += xx[i]*xx[(i+delta)%N];
    }
    
    return res/N;
}