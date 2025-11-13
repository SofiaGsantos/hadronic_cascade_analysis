#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

// ROOT headers
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TColor.h"

double eta(double px, double py, double pz) {
    double p = std::sqrt(px*px + py*py + pz*pz);
    if (std::abs(p - pz) != 0 ) {
        return 0.5 * std::log((p + pz) / (p - pz));
    } else {
        return -1000; 
    }        
}

void plot(const std::map<double, double>& Kstar_by_time,
          const std::map<double, double>& K_by_time) {
    
    // ---  K*0 ---
    std::vector<double> time_Kstar, counts_Kstar;
    for (const auto& par : Kstar_by_time) {
        time_Kstar.push_back(par.first);
        counts_Kstar.push_back(par.second);
    }

    int N_kstar = time_Kstar.size();
    TGraph* g_kstar = new TGraph(N_kstar, &time_Kstar[0], &counts_Kstar[0]);
    g_kstar->SetMarkerStyle(20);
    g_kstar->SetMarkerColor(kBlue);
    g_kstar->SetLineColor(kBlue);
    g_kstar->SetTitle("K*^{0} Temporal Evolution");
    g_kstar->GetXaxis()->SetTitle("Time (fm/c)");
    g_kstar->GetYaxis()->SetTitle("Counts");

    // ---  K- ---
    std::vector<double> tempos_k, contagens_k;
    for (const auto& par : K_by_time) {
        tempos_k.push_back(par.first);
        contagens_k.push_back(par.second);
    }

    int N_k = tempos_k.size();
    TGraph* g_k = new TGraph(N_k, &tempos_k[0], &contagens_k[0]);
    g_k->SetMarkerStyle(21);
    g_k->SetMarkerColor(kRed);
    g_k->SetLineColor(kRed);
    g_k->SetTitle("K^{-} Temporal Evolution");
    g_k->GetXaxis()->SetTitle("Time (fm/c)");
    g_k->GetYaxis()->SetTitle("Counst");

    TCanvas* c = new TCanvas("c", "K^{-} and K*^{0} Temporal Evolution", 1200, 600);
    gStyle->SetOptStat(0);
    c->Divide(2, 1); 

    c->cd(1);
    g_kstar->Draw("APL");


    c->cd(2);
    g_k->Draw("APL");

    c->SaveAs("kstar_kminus_.pdf");
}




int main() {
    std::string lista_arquivos = "lista.txt";
    std::ifstream lista(lista_arquivos.c_str());
    if (!lista.is_open()) {
        std::cerr << "Error to open list files.\n";
        return 1;
    }

    std::map<double, int> Kstar_by_time;
    std::map<double, int> K_by_time;
    int evento = 0;

    const int PDG_KSTAR0 = 313;
    const int PDG_K_MINUS = -321;

    std::string nome_arquivo;
    while (std::getline(lista, nome_arquivo)) {
        if (nome_arquivo.empty()) continue;

        std::ifstream arquivo(nome_arquivo.c_str());
        if (!arquivo.is_open()) {
            std::cerr << "file don't found: " << nome_arquivo << "\n";
            continue;
        }

        std::string linha;
        while (std::getline(arquivo, linha)) {
            if (linha.empty()) continue;

            std::istringstream iss(linha);
            std::vector<std::string> colunas;
            std::string item;
            while (iss >> item)
                colunas.push_back(item);
            if (colunas.size() ==4){
                evento++;
            }

            if (colunas.size() == 12) {
                try {
                    std::string tempo_str = colunas[11];
                    std::replace(tempo_str.begin(), tempo_str.end(), 'D', 'E');
                    double tempo = std::stod(tempo_str);
                    int pdg = std::stoi(colunas[1]);
                    double px = std::stod(colunas[3]);
                    double py = std::stod(colunas[4]);
                    double pz = std::stod(colunas[5]);

                    if (pdg == PDG_KSTAR0 && eta(px, py, pz) < 0.5 ){ 
                        Kstar_by_time[tempo]++;
                    }
                    if (pdg == PDG_K_MINUS  && eta(px, py, pz) < 0.5){ 
                        K_by_time[tempo]++;
                    }
                } catch (...) {
                    continue;
                }
            }
        }

        arquivo.close();

    }

    std::cout << "Número total de samples: " << evento << "\n";
    std::map<double, double> Kstar_by_time_normalizado;
    std::map<double, double> K_by_time_normalizado;

    for (const auto& par : Kstar_by_time) {
        Kstar_by_time_normalizado[par.first] = static_cast<double>(par.second) / evento;
    }

    for (const auto& par : K_by_time) {
        K_by_time_normalizado[par.first] = static_cast<double>(par.second) / evento;
    }


    plot(Kstar_by_time_normalizado, K_by_time_normalizado);

    return 0;
}

