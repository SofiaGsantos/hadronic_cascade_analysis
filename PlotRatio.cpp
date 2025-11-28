#include <iostream>
#include <fstream>
#include <vector>
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace std;

int main() {

    vector<double> multiplicity;
    vector<double> ratios;

    // ----------------- read multiplicity -----------------
    {
        ifstream fin("multiplicity.txt");
        if (!fin.is_open()) {
            cerr << "Erro ao abrir multiplicity.txt\n";
            return 1;
        }
        double x;
        while (fin >> x) {
            multiplicity.push_back(x);
        }
    }

    // ----------------- read ratios -----------------
    {
        ifstream fin("ratios.txt");
        if (!fin.is_open()) {
            cerr << "Erro ao abrir ratios.txt\n";
            return 1;
        }
        double y;
        while (fin >> y) {
            ratios.push_back(y);
        }
    }

    // ----------------- check sizes -----------------
    if (multiplicity.size() != ratios.size()) {
        cerr << "Listas com tamanhos diferentes!\n";
        return 1;
    }

    int N = multiplicity.size();

    // ----------------- plot -----------------
    TCanvas *c1 = new TCanvas("c1", "Ratio vs Multiplicity", 900, 700);

    TGraph *gr = new TGraph(N, multiplicity.data(), ratios.data());
    gr->SetTitle("K*/K^{-} ;(dN_{ch}/d#eta)^{1/3};K*/K^{-}");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);

    gr->Draw("AP"); 

    c1->SetGrid();
    c1->SaveAs("ratio_vs_multiplicity.png");

    return 0;
}
