// This analysis is inspired on investigations made in https://arxiv.org/abs/2101.07302 
//----------------------------------------------------------------------------------------
// The objective is inderstanding the machanisms behind K*/K- supression
// specifically, we want to know de role played by rescatering in the observed K*/K- supression

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>

// ROOT
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"

using namespace std;

//identification of columns
// --------------------------------------------------------------------------------
const int ID   = 10;
const int PDG  = 9;
// --------------------------------------------------------------------------------

//support structure
struct ParticleInfo {
    int id{};
    int pdg{};
};


// functions

// 1) helps to identify a decay
// in the file, an interaction is identify by "# interaction" in the beginning of the line
static inline bool starts_with(const string &s, const string &prefix) {         //const string &s > Read-only reference
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

//2) split_tokens take a string (line) and split it into words separated by spaces
//ex:   "  K*  313   resonance   0.892 "  
//     ["K*", "313", "resonance", "0.892"]
static vector<string> split_tokens(const string &line) {
    istringstream iss(line);
    vector<string> tk;
    string t;
    while (iss >> t) tk.push_back(t);
    return tk;
}

//3) take a line that contains something like:
//in 3 out 2 type 10
//and extract the numeric values.
static bool parse_header_fields(const string& hdr, int& inN, int& outN, int& type) {
    inN = outN = -1; type = -1;
    auto tk = split_tokens(hdr);
    for (size_t i = 0; i + 1 < tk.size(); ++i) {
        if (tk[i] == "in")   { inN  = stoi(tk[i+1]); }
        if (tk[i] == "out")  { outN = stoi(tk[i+1]); }
        if (tk[i] == "type") { type = stoi(tk[i+1]); }
    }
    return (inN >= 0 && outN >= 0 && type >= 0);
}

static int parse_event_line(const string& hdr) {
    int event_number = -1;
    auto tk = split_tokens(hdr);
    for (size_t i = 0; i + 1 < tk.size(); ++i) {
        if (tk[i] == "event" && tk[i+2] == "in")  { event_number  = stoi(tk[i+1]); }
    }
    return (event_number);
}

static int parse_Kminus_final(const string& hdr, int num) {
    auto tk = split_tokens(hdr);
    int finalparticles = 0;
    for (size_t i = 0; i + 1 < tk.size(); ++i) {
        if (tk[i] ==  "out" ) { finalparticles = stoi(tk[i+1]); }
    }
    return (finalparticles);
}

//4)add information to the "particle" structure according to the particle's line
static bool parse_particle_row(const string& line, ParticleInfo& p) {
    auto cols = split_tokens(line);
    if ((int)cols.size() <= ID) return false;
    try {
        p.id        = stoi(cols[ID]);
        p.pdg       = stoi(cols[PDG]);
    } catch (...) { return false; }
    return true;
}

//5) draw the graph of the ratio K*/K-
// void plot(const map<double, double>& probability) {
//     if (probability.empty()) {
//         cerr << "Warning: Empty probability map!\n";
//         return;
//     }
//     vector<double> time, prob;
//     time.reserve(probability.size());
//     prob.reserve(probability.size());
//     for (const auto& kv : probability) {
//         time.push_back(kv.first);
//         prob.push_back(kv.second);
//     }

//     TCanvas* c = new TCanvas("c", " K^{*0}/K^{-}", 1200, 600);
//     gStyle->SetOptStat(0);

//     TGraph* g = new TGraph(time.size(), &time[0], &prob[0]);
//     g->SetMarkerStyle(20);
//     g->SetMarkerColor(kRed);
//     g->GetXaxis()->SetTitle("Time [fm/c]");
//     g->GetYaxis()->SetTitle("#frac{K^{*0}}{K^{-}}");
//     g->GetXaxis()->SetLimits(0, 120);
//     g->Draw("AP");

//     c->SaveAs("Ratio.pdf");
// }

int main() {
    ifstream list("test.txt");
    if (!list.is_open()) {
        cerr << "Error to open files list.\n";
    }

    string filename;
    map <int,double> ratios;

    while (getline(list, filename)) {
        if (filename.empty()) continue;

        ifstream fin(filename);
        if (!fin.is_open()) {
            cerr << "Error to open: " << filename << "\n";
            continue;
        }
        string line;
        int event = 0;
        int j = 0;
        vector<vector<int>> active_daughter_ids;            // IDs of daughters still "alive" for detection

        while (getline(fin, line)) {

            if ( event == j ){ 


                // ---------------------------------  large block for detecting decays  --------------------------------------
                if (starts_with(line, "# interaction")) { 

                    // read header
                    int inN=-1, outN=-1, type=-1;
                    if (!parse_header_fields(line, inN, outN, type)) continue;
                    //cout << "encontrou linha de interação \n";

                    // read inN + outN lines 
                    vector<int> partspdg;
                    vector<int> partsid;
                    partspdg.reserve(inN + outN);
                    partsid.reserve(inN + outN);

                    for (int i = 0; i < inN + outN; i++) {
                        string pline;
                        if (!getline(fin, pline)) break;
                        ParticleInfo p;
                        if (!parse_particle_row(pline, p)) continue;
                        if (type == 5 ) {
                            int pdg = p.pdg;
                            int id = p.id;
                            partspdg.push_back(pdg);
                            partsid.push_back(id);
                        }
                    
                    }

                    if (partspdg[0] == 313 || partspdg[0] == -313 || partspdg[0] == 323 || partspdg[0] == -323) {  //only K*(892): more abundant that decays into Kπ
                        // daughter
                        int d1 = partsid[1];
                        int d2 = partsid[2];

                        // tracking daughters: map and vector of active IDs
                        active_daughter_ids.push_back({d1, d2});
                    } 
                }   


                // ---------------------------------  large block to count the particles that reached the detectors  --------------------------------------
            
                if (starts_with(line, "# event " + to_string(event) + " out")){
                    vector <int> final_Kminus;             //IDs of K- that reach the detectors

                    for (int i=0; i < parse_Kminus_final(line, event); ){
                        string pline;
                        if (!getline(fin, pline)) break;
                        ParticleInfo p;
                        if (!parse_particle_row(pline, p)) continue;         
                        if (p.pdg == 321 || p.pdg == -321){                  
                            final_Kminus.push_back(p.id);
                            int id = p.id;

                            for (size_t k = 0; k < active_daughter_ids.size(); ) {
                                if (active_daughter_ids[k][0] == id || active_daughter_ids[k][1] == id) {
                                    active_daughter_ids.erase(active_daughter_ids.begin() + k);
                                } else {
                                    k++;   
                                }
                            }
                        }
                    }

                    double ratio = (double(active_daughter_ids.size()) / double(final_Kminus.size()));
                    ratios[j] = ratio;
                    active_daughter_ids.clear();
                    j+=1;
                }
            }
        }
    }
    
    //------------------------ media on the oversamples ----------------------------
    double sum = 0.0;
    for(const auto& par : ratios) {
        sum += par.second;     
    }
    double average = sum / ratios.size();
    cout << "a razão é " << average;
    return(average);
} 
    