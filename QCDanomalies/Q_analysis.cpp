/*
 * Q_analysis.cpp
 * This program analyzes the QCD anomalies data from four different WFtime 
 * points and generates histograms for each time point. It calculates the mean and 
 * variance of the Chi values at t100, and saves the results to a file.
 * Data in define as to be specify fram log file.
 * 
 * To compile this code, use the following command:
 * g++ Q_analysis.cpp $(root-config --cflags --libs) -o Q_analysis
 *  
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include "TH1F.h"
#include "TCanvas.h"

#define N0 20
#define N1 20
#define N2 20
#define N3 20

#define t_0 6.2

using namespace std;

int main(void) {

    TH1F *h000 = new TH1F("h000", "Q(t000)", 100, -3, 3);
    TH1F *h010 = new TH1F("h010", "Q(t010)", 100, -3, 3);
    TH1F *h025 = new TH1F("h025", "Q(t025)", 100, -3, 3);
    TH1F *h100 = new TH1F("h100", "Q(t100)", 100, -3, 3);
    TCanvas *c = new TCanvas("c", "QCD Anomalies Analysis");

    ifstream infile_000("Qt000.txt");
    ifstream infile_010("Qt010.txt");
    ifstream infile_025("Qt025.txt");
    ifstream infile_100("Qt100.txt");
    if (!infile_000.is_open() || !infile_010.is_open() || !infile_025.is_open() || !infile_100.is_open()) {
        cerr << "Error opening input files" << endl;
        return 1;
    }

    double value000, value010, value025, value100;
    double sum = 0, sumsq;
    int count = 0;
    while (infile_000 >> value000 && infile_010 >> value010 && infile_025 >> value025 && infile_100 >> value100) {
        h000->Fill(value000);
        h010->Fill(value010);
        h025->Fill(value025);
        h100->Fill(value100);

        sum += value100*value100;
        sumsq += value100*value100*value100*value100;
        count++;
        
    }

    double mean = sum/count;
    double variance = (sumsq/count - (mean*mean));
    double meanerr = sqrt(variance/count);

    infile_000.close();
    infile_010.close();
    infile_025.close();
    infile_100.close();

    h000->SetStats(0);
    h010->SetStats(0);
    h025->SetStats(0);
    h100->SetStats(0);

    h000->SetFillColor(kBlue);
    h010->SetFillColor(kBlue);
    h025->SetFillColor(kBlue);
    h100->SetFillColor(kBlue);

    h000->SetLineColor(kBlue);
    h010->SetLineColor(kBlue);
    h025->SetLineColor(kBlue);
    h100->SetLineColor(kBlue);

    c->Divide(2,2);

    c->cd(1);
    h000->Draw();

    c->cd(2);
    h010->Draw();

    c->cd(3);
    h025->Draw();

    c->cd(4);
    h100->Draw();

    c->SaveAs("Q20_WF.png");

    double chi = t_0*t_0*mean/(N0*N1*N2*N3);
    double chierr = t_0*t_0*meanerr/(N0*N1*N2*N3); 

    FILE *outfile = fopen("chi_extrapolation.txt", "a");
    if (!outfile) {
        fprintf(stderr, "Error opening output file chi_extrapolation.txt\n");
        return 1;
    }
    fprintf(outfile, "%.17lf\t%.17lf%.17lf\t%.17lf\n", chi, chierr, t_0, .0);
    fclose(outfile);

    delete h000;
    delete h010;
    delete h025;
    delete h100;
    delete c;

    return 0;   
}