/*
 * chi_analysis.cpp
 * This program serve to contiuum extrapolation of the chi value.
 * 
 * To compile this code, use the following command:
 * g++ chi_analysis.cpp $(root-config --cflags --libs) -o chi_analysis
 * 
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include "TGraphErrors.h"
#include "TAxis.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"

using namespace std;

int main(void){

    ifstream infile("chi_extrapolation.txt");
    vector<double> x, y, ex, ey;
    double xv, yv, exv, eyv;
    while (infile >> yv >> eyv >> xv >> exv) {
        x.push_back(1/xv);
        y.push_back(yv);
        ex.push_back(exv);
        ey.push_back(eyv);
    }
    infile.close();

    int n = x.size();
    TGraphErrors *g = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
    TCanvas *c = new TCanvas("c", "QCD Anomalies Analysis");

    g->SetTitle("");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue);
    g->Draw("AP");
    g->GetXaxis()->SetTitle("a^{2}/t_{0}");
    g->GetYaxis()->SetTitle("#chi");
    g->GetXaxis()->SetRangeUser(-0.1, 0.5);

    TF1 *f_fit = new TF1("fit", "pol1");
    TFitResultPtr fit_result = g->Fit(f_fit, "S");
    if(fit_result->IsValid() == true){
        cout << "Il fit è avvenuto con successo \n" << endl;
    }

    f_fit->SetLineColor(kRed);
    f_fit->Draw("same");

    c->SaveAs("chi_extrapolation.png");

    delete f_fit;
    delete g;
    delete c;

    return 0;
}

