#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

// ROOT headers
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"

// --------------------------------------------------------------------
// plot the (regeneration-decay)  graph for a set of resonances
// --------------------------------------------------------------------
void plot_all(const std::map<std::string, std::map<double, double>>& data_map) {
    TCanvas* c = new TCanvas("c", "Gain/Loss Ratio", 1200, 600);
    gStyle->SetOptStat(0);

    int color_index = 2; 
    TLegend* legend = new TLegend(0.70, 0.10, 0.90, 0.30);

    bool first = true;
    for (const auto& kv : data_map) {
        const std::string& name = kv.first;
        const std::map<double, double>& series = kv.second;

        std::vector<double> time, ratio;
        time.reserve(series.size());
        ratio.reserve(series.size());
        for (const auto& par : series) {
            time.push_back(par.first);
            ratio.push_back(par.second);
        }

        if (time.empty()) {
            std::cerr << "Warning: Empty ratio series for " << name << "!\n";
            continue;
        }

        TGraph* graph = new TGraph((int)time.size(), time.data(), ratio.data());
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(color_index);
        graph->SetLineColor(color_index);
        graph->SetTitle("Net Gain Rate");
        graph->GetXaxis()->SetTitle("Simulation Time (fm/c)");
        graph->GetYaxis()->SetTitle("#frac{dN_{reg}}{dt} - #frac{dN_{dec}}{dt}");
        graph->GetXaxis()->SetLimits(0, 100);

        if (first) {
            graph->Draw("AP*");  
            first = false;
        } else {
            graph->Draw("P SAME");
        }

        legend->AddEntry(graph, name.c_str(), "lp");
        color_index++;
    }

    legend->Draw();
    c->SaveAs("net_gain_resonances_smash.pdf");
}

// --------------------------------------------------------------------
// normalization
// input: map<double,int> or map<double,double>.
// output: map<bin_center, sum/bin_size>
// --------------------------------------------------------------------

template <typename T>
std::map<double, double> bin_and_sum(const std::map<double, T>& input, double bin_size) {
    std::map<int, double> bins; 

    for (auto it = input.begin(); it != input.end(); ++it) {
        double time = it->first;
        double value = static_cast<double>(it->second);
        int bin_index = static_cast<int>(std::floor(time / bin_size));
        bins[bin_index] += value; 
    }

    std::map<double, double> output;
    for (auto it = bins.begin(); it != bins.end(); ++it) {
        int bin_index = it->first;
        double sum = it->second;
        double bin_center = (static_cast<double>(bin_index) + 0.5) * bin_size;
        double norm_sum = sum / bin_size; 
        output[bin_center] = norm_sum;
    }

    return output;
}

// --------------------------------------------------------------------
// calculate de net gain rate
// --------------------------------------------------------------------
std::map<double, double> netgain(const std::map<double, double>& gain_binned,
                                      const std::map<double, double>& decay_binned) {
    std::map<double, double> netgain;
    for (const auto& kv : decay_binned) {
        double t = kv.first;
        double d = kv.second;
        auto itg = gain_binned.find(t);
        double g = (itg != gain_binned.end()) ? itg->second : 0.0;
        netgain[t] = g - d;
        }
    
    return netgain;
}

int main() {
    //code input: files list
    std::string files_list = "teste.txt";
    std::ifstream list(files_list.c_str());
    if (!list.is_open()) {
        std::cerr << "Error to open files list.\n";
        return 1;
    }

    std::map<double, int> gainKstar;      std::map<double, int> decayKstar;
    std::map<double, int> gainRho;        std::map<double, int> decayRho;
    std::map<double, int> gainPhi;        std::map<double, int> decayPhi;
    std::map<double, int> gainLambdastar; std::map<double, int> decayLambdastar;
    std::map<double, int> gainSigmastar;  std::map<double, int> decaySigmastar;

    const double bin_size = 1.0;

    std::string file;
    while (std::getline(list, file)) {
        if (file.empty()) continue;

        std::ifstream text_file(file.c_str());
        if (!text_file.is_open()) {
            std::cerr << "error to open: " << file << "\n";
            continue;
        }

        std::string line;
        while (std::getline(text_file, line)) {
            if (line.empty()) continue;

            std::istringstream iss(line);
            std::vector<std::string> columns;
            std::string item;
            while (iss >> item) columns.push_back(item);

            //a line containing particle information has 22 columns
            if (columns.size() == 22) {
                try {
                    double current_time   = std::stod(columns[0]);
                    double last_coll_time = std::stod(columns[17]);
                    int    pdg            = std::stoi(columns[9]);
                    double formtime       = std::stod(columns[13]);

                    // K* (313)
                    if (pdg == 313) {
                        if (last_coll_time == current_time) decayKstar[current_time]++;
                        if (formtime       == current_time) gainKstar[current_time]++;
                    }
                    // Rho (113)
                    if (pdg == 113) {
                        if (last_coll_time == current_time) decayRho[current_time]++;
                        if (formtime       == current_time) gainRho[current_time]++;
                    }
                    // Phi (333)
                    if (pdg == 333) {
                        if (last_coll_time == current_time) decayPhi[current_time]++;
                        if (formtime       == current_time) gainPhi[current_time]++;
                    }
                    // Lambda* (3124)
                    if (pdg == 3124) {
                        if (last_coll_time == current_time) decayLambdastar[current_time]++;
                        if (formtime       == current_time) gainLambdastar[current_time]++;
                    }
                    // Sigma* (3214)
                    if (pdg == 3214) {
                        if (last_coll_time == current_time) decaySigmastar[current_time]++;
                        if (formtime       == current_time) gainSigmastar[current_time]++;
                    }
                } catch (const std::exception& e) {
                    std::cerr << "Error parsing line: " << e.what() << "\n";
                }
            }
        }
        text_file.close();
    }

    // ----------------------------------------------------------------
    // 1) binning 
    // ----------------------------------------------------------------
    auto binned_gainKstar   = bin_and_sum(gainKstar,      bin_size);
    auto binned_decayKstar  = bin_and_sum(decayKstar,     bin_size);

    auto binned_gainRho     = bin_and_sum(gainRho,        bin_size);
    auto binned_decayRho    = bin_and_sum(decayRho,       bin_size);

    auto binned_gainPhi     = bin_and_sum(gainPhi,        bin_size);
    auto binned_decayPhi    = bin_and_sum(decayPhi,       bin_size);

    auto binned_gainLambda  = bin_and_sum(gainLambdastar, bin_size);
    auto binned_decayLambda = bin_and_sum(decayLambdastar,bin_size);

    auto binned_gainSigma   = bin_and_sum(gainSigmastar,  bin_size);
    auto binned_decaySigma  = bin_and_sum(decaySigmastar, bin_size);

    // ----------------------------------------------------------------
    // 2) calculation of the net gain rate
    // ----------------------------------------------------------------
    std::map<double, double> netgainKstar   = netgain(binned_gainKstar,  binned_decayKstar);
    std::map<double, double> netgainRho     = netgain(binned_gainRho,    binned_decayRho);
    std::map<double, double> netgainPhi     = netgain(binned_gainPhi,    binned_decayPhi);
    std::map<double, double> netgainLambda  = netgain(binned_gainLambda, binned_decayLambda);
    std::map<double, double> netgainSigma   = netgain(binned_gainSigma,  binned_decaySigma);

    // prepare for plotting
    std::map<std::string, std::map<double,double>> data_map;
    data_map["K*"]       = netgainKstar;
    data_map["#Rho"]     = netgainRho;
    data_map["#Phi"]     = netgainPhi;
    data_map["#Lambda*"] = netgainLambda;
    data_map["#Sigma*"]  = netgainSigma;

    plot_all(data_map);

    return 0;
}


