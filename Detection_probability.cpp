#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <unordered_map>
#include <numeric>
#include <algorithm>

// ROOT
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TStyle.h"

// --------------------------------------------------------------------------------
const int TIME_COL     = 0;
const int ID_COL       = 10;
const int PDG_COL      = 9;
const int NCOLL_COL    = 12;
const int FORMTIME_COL = 13;
// --------------------------------------------------------------------------------

struct ParticleInfo {
    double time{};
    int id{};
    int pdg{};
    int ncoll{};
    double form_time{};
};

struct DecayRecord {
    double time{};
    int d1{};
    int d2{};
    bool counted{false}; 
};

// functions


// 1) helps to identify a decay
// in the file, an interaction is identify by "# interaction" in the beginning of the line
static inline bool starts_with(const std::string &s, const std::string &prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

//ltrim remove o espaço em branco do início de uma string
static inline std::string ltrim(const std::string &s) {
    size_t i = s.find_first_not_of(" \t\r\n");
    return (i == std::string::npos) ? std::string() : s.substr(i);
}

//2) split_tokens take a string (line) and split it into words separated by spaces
//ex:   "  K*  313   resonance   0.892 "  
//     ["K*", "313", "resonance", "0.892"]
static std::vector<std::string> split_tokens(const std::string &line) {
    std::istringstream iss(line);
    std::vector<std::string> tk;
    std::string t;
    while (iss >> t) tk.push_back(t);
    return tk;
}

//3) take a line that contains something like:
//in 3 out 2 type 10
//and extract the numeric values.
static bool parse_header_fields(const std::string& hdr, int& inN, int& outN, int& type) {
    inN = outN = -1; type = -1;
    auto tk = split_tokens(hdr);
    for (size_t i = 0; i + 1 < tk.size(); ++i) {
        if (tk[i] == "in")   { inN  = std::stoi(tk[i+1]); }
        if (tk[i] == "out")  { outN = std::stoi(tk[i+1]); }
        if (tk[i] == "type") { type = std::stoi(tk[i+1]); }
    }
    return (inN >= 0 && outN >= 0 && type >= 0);
}

static bool parse_particle_row(const std::string& line, ParticleInfo& p) {
    auto cols = split_tokens(line);
    int need = std::max({TIME_COL, ID_COL, PDG_COL, NCOLL_COL, FORMTIME_COL});
    if ((int)cols.size() <= need) return false;
    try {
        p.time      = std::stod(cols[TIME_COL]);
        p.id        = std::stoi(cols[ID_COL]);
        p.pdg       = std::stoi(cols[PDG_COL]);
        p.ncoll     = std::stoi(cols[NCOLL_COL]);
        p.form_time = std::stod(cols[FORMTIME_COL]);
    } catch (...) { return false; }
    return true;
}

// plot
void plot(const std::map<double, double>& probability) {
    if (probability.empty()) {
        std::cerr << "Warning: Empty probability map!\n";
        return;
    }
    std::vector<double> time, prob;
    time.reserve(probability.size());
    prob.reserve(probability.size());
    for (const auto& kv : probability) {
        time.push_back(kv.first);
        prob.push_back(kv.second);
    }

    TCanvas* c = new TCanvas("c", " K^{*0}/K^{-}", 1200, 600);
    gStyle->SetOptStat(0);

    TGraph* g = new TGraph(time.size(), &time[0], &prob[0]);
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kRed);
    g->GetXaxis()->SetTitle("Time [fm/c]");
    g->GetYaxis()->SetTitle("#frac{K^{*0}}{K^{-}")
    g->GetXaxis()->SetLimits(0, 120);
    g->Draw("AP");

    c->SaveAs("Ratio.pdf");
}

int main() {
    std::ifstream list("teste.txt");
    if (!list.is_open()) {
        std::cerr << "Error to open files list (teste.txt).\n";
    }

    const double bin_size = 1.0;
    std::map<double, int> total;    // total decays
    std::map<double, int> detected; // detected decays

    std::string filename;
    while (std::getline(list, filename)) {
        if (filename.empty()) continue;

        std::ifstream fin(filename);
        if (!fin.is_open()) {
            std::cerr << "Error to open: " << filename << "\n";
            continue;
        }

        // Estructures by file
        std::vector<DecayRecord> decays;                 // decay list
        std::unordered_map<int,int> id_to_decay_index;   // id (daughter) -> index in decays
        std::vector<int> active_daughter_ids;            // IDs of daughters still "alive" for detection

        std::string line;
        while (std::getline(fin, line)) {
            if (!starts_with(line, "# interaction")) continue;

            // read header
            int inN=-1, outN=-1, type=-1;
            if (!parse_header_fields(line, inN, outN, type)) continue;

            // read inN + outN lines 
            std::vector<ParticleInfo> parts;
            parts.reserve(inN + outN);

            for (int i = 0; i < inN + outN; ) {
                std::streampos pos = fin.tellg();
                std::string pline;
                if (!std::getline(fin, pline)) break;
                auto trimmed = ltrim(pline);
                if (trimmed.empty()) continue;
                if (starts_with(trimmed, "# interaction")) {
                    // for the case that there ir a header in a wrong place
                    fin.seekg(pos);
                    break;
                }

                ParticleInfo p;
                if (!parse_particle_row(pline, p)) continue;

                // *** check new interaction ***
                auto it = id_to_decay_index.find(p.id);
                if (it != id_to_decay_index.end()) {
                    int idx = it->second;
                    if (idx >= 0 && idx < (int)decays.size() && decays[idx].counted) {
                        // if interacts, decreases the count by 1 
                        detected[decays[idx].time]--;
                        decays[idx].counted = false;

                        // remove both IDs 
                        int d1 = decays[idx].d1, d2 = decays[idx].d2;
                        id_to_decay_index.erase(d1);
                        id_to_decay_index.erase(d2);

                        ;
                        active_daughter_ids.erase(
                            std::remove(active_daughter_ids.begin(), active_daughter_ids.end(), d1),
                            active_daughter_ids.end()
                        );
                        active_daughter_ids.erase(
                            std::remove(active_daughter_ids.begin(), active_daughter_ids.end(), d2),
                            active_daughter_ids.end()
                        );
                    }
                }

                parts.push_back(p);
                ++i;
            }

            // if DECAY (type==5) and K*
            if (type == 5 && (int)parts.size() >= 3) {
                // mother
                const ParticleInfo& mother = parts[0];
                int mpdg = mother.pdg;

                if (mpdg == 313 || mpdg == -313 || mpdg == 323 || mpdg == -323) {
                    // daughter
                    const ParticleInfo& d1 = parts[1];
                    const ParticleInfo& d2 = parts[2];

                    double decay_time = mother.time;
                    total[decay_time]++;

                    // assumes detectable
                    detected[decay_time]++;

                    // registers decay
                    DecayRecord rec;
                    rec.time = decay_time;
                    rec.d1 = d1.id;
                    rec.d2 = d2.id;
                    rec.counted = true;
                    int idx = (int)decays.size();
                    decays.push_back(rec);

                    // tracking daughters: map and vector of active IDs
                    id_to_decay_index[d1.id] = idx;
                    id_to_decay_index[d2.id] = idx;
                    active_daughter_ids.push_back(d1.id);
                    active_daughter_ids.push_back(d2.id);
                }
            }
        } 

        fin.close();
    } 

// Raw K*/K- over time
std::map<double, double> raw_probability;
for (std::map<double,int>::const_iterator it = total.begin(); it != total.end(); ++it) {
    double t = it->first;
    int dec = it->second;
    int det = (detected.count(t) ? detected[t] : 0);
    raw_probability[t] = (dec > 0) ? double(det) / double(dec) : 0.0;
}

// Binning 
std::map<int, std::vector<double>> bins;
for (std::map<double,double>::const_iterator it = raw_probability.begin(); it != raw_probability.end(); ++it) {
    double t = it->first;
    double p = it->second;
    int b = (t >= 0.0) ? int(t / bin_size) : int((t - bin_size + 1.0) / bin_size);
    bins[b].push_back(p);
}

std::map<double, double> binned_probability;
for (std::map<int,std::vector<double>>::const_iterator it = bins.begin(); it != bins.end(); ++it) {
    int b = it->first;
    const std::vector<double>& vec = it->second;
    double avg = vec.empty() ? 0.0 : std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
    double center = (b + 0.5) * bin_size;
    binned_probability[center] = avg;
}

    plot(binned_probability);
    return 0;
}






