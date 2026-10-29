/*
*******************************************************************************
*
* File check1.c
*
* Measure the Metropolis acceptance rate and autocorrelation time for
* several proposal widths.
*
*******************************************************************************/
#define MAIN_PROGRAM
#define CHECK1_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"
#include "random.h"
#include "ho.h"

static int MetropolisCheck(double proposal_width)
{
	int i;
	int accepted;
	double dx;
	double dS;
	double r[2 * N];

	accepted = 0;
	ranlxd(r, 2 * N);

	for (i = 0; i < N; i++) {
		dx = 2.0 * proposal_width * (r[i] - 0.5);
		dS = DeltaAction(xx[i] + dx, i);

		if (dS <= 0.0 || r[i + N] < exp(-dS)) {
			xx[i] += dx;
			accepted++;
		}
	}

	return accepted;
}

int main(void)
{
	int i, j, k; // Markov time, physical time, proposal index
    int maxLag;
	int Nproposal, Nsweeps;
	int accepted;
	double acceptance_rate;
	double total_accepted;
	double min_width, max_width, step_log;
	double *proposal_widths;
	double *Correlator;
	double MeanCorrelator[1];
	double Gamma0;
	FILE *fAcceptance, *fTau;

	Nsweeps = 10000;
	Nproposal = 100;
	maxLag = Nsweeps / 10;
	if (maxLag > 100) {
		maxLag = 100;
	}
	min_width = 0.1;
	max_width = 100.0;
	step_log = (log(max_width) - log(min_width)) / (Nproposal - 1);

    proposal_widths = (double *)malloc(Nproposal * sizeof(double));
    if (proposal_widths == NULL) {
        perror("Error allocating memory for proposal widths");
        return EXIT_FAILURE;
    }

	for (k = 0; k < Nproposal; k++) {
		proposal_widths[k] = exp(log(min_width) + k * step_log);
	}

	fAcceptance = fopen("MetropolisAcceptance.txt", "w");
	if (fAcceptance == NULL) {
		perror("Error opening MetropolisAcceptance.txt");
		free(proposal_widths);
		return EXIT_FAILURE;
	}

	fTau = fopen("MetropolisTau.txt", "w");
	if (fTau == NULL) {
		perror("Error opening MetropolisTau.txt");
		free(proposal_widths);
		fclose(fAcceptance);
		return EXIT_FAILURE;
	}

    Correlator = (double *)malloc(Nsweeps * sizeof(double));
    if (Correlator == NULL) {
        perror("Error allocating memory for correlator data");
        free(proposal_widths);
        fclose(fAcceptance);
        fclose(fTau);
        return EXIT_FAILURE;
    }

	for (k = 0; k < Nproposal; k++) {
		rlxd_init(1, SEED);

		for (j = 0; j < N; j++) {
			xx[j] = 0.0;
		}

		for (i = 0; i < NTHERM; i++) {
			MetropolisCheck(proposal_widths[k]);
		}

		total_accepted = 0.0;
		for (i = 0; i < Nsweeps; i++) {
			accepted = MetropolisCheck(proposal_widths[k]);
			total_accepted += accepted;
            Correlator[i] = ConfCorrelator(0);
		}

		MeanCorrelator[0] = 0.0;
        for (i = 0; i < Nsweeps; i++) {
            MeanCorrelator[0] += Correlator[i];
        }
        MeanCorrelator[0] /= Nsweeps;

        Gamma0 = Gamma(Correlator, MeanCorrelator, Nsweeps, 0, 0);
        for (i = 0; i < maxLag; i++) {
            fprintf(fTau, "%.15e%s", Gamma(Correlator, MeanCorrelator, Nsweeps, i, 0) / Gamma0,
                i == maxLag - 1 ? "\n" : " ");
        }

		acceptance_rate = total_accepted / ((double)Nsweeps * N);
		fprintf(fAcceptance, "%.15e %.15e\n",
				proposal_widths[k], acceptance_rate);
	}

	free(proposal_widths);
	fclose(fAcceptance);
	fclose(fTau);

	return EXIT_SUCCESS;
}
