//auxiliary functions for analysis

#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;   

//identification of columns
// --------------------------------------------------------------------------------
const int TIME = 0;
const int PX   = 6;
const int PY   = 7;
const int PZ   = 8;
const int Q    = 11;
const int ID   = 10;
const int PDG  = 9;
// --------------------------------------------------------------------------------

//support structure
struct ParticleInfo {
    double time{};
    double px{};
    double py{};
    double pz{};
    int charge{};
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


//4) extracts the final number of particles from each event
static int parse_Kminus_final(const string& hdr) {
    auto tk = split_tokens(hdr);
    int finalparticles = 0;
    for (size_t i = 0; i + 1 < tk.size(); ++i) {
        if (tk[i] ==  "out" ) { finalparticles = stoi(tk[i+1]); }
    }
    return (finalparticles);
}

//5) add information to the "particle" structure according to the particle's line
static bool parse_particle_row(const string& line, ParticleInfo& p) {
    auto cols = split_tokens(line);
    if ((int)cols.size() <= Q) return false;
    try {
        p.time     = stod(cols[TIME]);
        p.px       = stod(cols[PX]);
        p.py       = stod(cols[PY]);
        p.pz       = stod(cols[PZ]);
        p.charge   = stoi(cols[Q]);
        p.id       = stoi(cols[ID]);
        p.pdg      = stoi(cols[PDG]);
    } catch (...) { return false; }
    return true;
}

//6) calculation of eta
double eta(double px, double py, double pz)
{
    double p = sqrt(px*px + py*py + pz*pz);
    return 0.5 * log( (p + pz) / (p - pz) );
}

//7) calculation of pT
double compute_pT(double px, double py) {
    return std::sqrt(px*px + py*py);
}

#endif