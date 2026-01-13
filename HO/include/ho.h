#ifndef HO_H
#define HO_H

extern double Action(void);

extern double DeltaAction(double varx, int n);

extern void Metropolis(void);

extern double ConfCorrelator(int delta);

#endif
