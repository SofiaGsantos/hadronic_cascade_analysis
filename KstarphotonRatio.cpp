// The objective is inderstanding the machanisms behind ressonance supression
// specifically, we want to know identify the supression of a different channel (K* -> K⁰ + γ )

#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include "utils.h"

using namespace std;

const double pTmin = 0.2;   
const double pTmax = 5.0;

int main() {
    ifstream list("test.txt");
    if (!list.is_open()) {
        cerr << "Error to open files list.\n";
    }

    string filename;
    vector <double> ratios;

    while (getline(list, filename)) {
        if (filename.empty()) continue;

        vector <double> y_values;
        ifstream fin(filename);
        if (!fin.is_open()) {
            cerr << "Error to open: " << filename << "\n";
            continue;
        }
        string line;
        int event_number = 0;
        vector<vector<int>> daughter_ids;            

        while (getline(fin, line)) {

            // ---------------------------------  large block for detecting decays  --------------------------------------
            if (starts_with(line, "# interaction")) { 

                // read header
                int inN=-1, outN=-1, type=-1;
                if (!parse_header_fields(line, inN, outN, type)) continue;

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
                
                if (partspdg[0] == 313 || partspdg[0] == -313 || partspdg[0] == 323 || partspdg[0] == -323 || partspdg[0] == 10313 || partspdg[0] == -10313 || partspdg[0] == 20313 || partspdg[0] == - 20313 ) {  //only K*(892): more abundant that decays into Kπ
                    // daughter
                    int d1 = partsid[1];
                    int d2 = partsid[2];

                    // tracking daughters: map and vector of active IDs
                    daughter_ids.push_back({d1, d2});
                } 
            }   


                // ---------------------------------  large block to count the particles that reached the detectors  --------------------------------------
            
             if (starts_with(line, "# event " + to_string(event_number) + " out")){
                vector <int> final_photon;             //IDs of K- that reach the detectors
                vector <int> final_daughters; 
                cout << "found the end of the event \n";
                for (int i=0; i < parse_Kminus_final(line); ){
                    string pline;
                    if (!getline(fin, pline)) break;
                    ParticleInfo p;
                    if (!parse_particle_row(pline, p)) continue;  
                    double et = eta(p.px, p.py, p.pz);
                    double pT = compute_pT (p.px, p.py);

                    if (p.pdg == 311 || p.pdg == -311 || p.pdg == 22  && ( et > -0.5 && et < 0.5) && (pT >= pTmin && pT <= pTmax )){ 
                        final_daughters.push_back(p.id);              
                    }
                    
                    if (p.pdg == 22 || p.pdg == 22 && ( et > -0.5 && et < 0.5) && (pT >= pTmin && pT <= pTmax )){
                        final_photon.push_back(p.pdg);
                    }
                }

                size_t k = 0;
                while (k < daughter_ids.size()) {

                    int id1 = daughter_ids[k][0];
                    int id2 = daughter_ids[k][1];

                    bool ok1 = std::find(final_daughters.begin(), final_daughters.end(), id1) != final_daughters.end();
                    bool ok2 = std::find(final_daughters.begin(), final_daughters.end(), id2) != final_daughters.end();

                    if (ok1 && ok2) {
                        // both daughter survived
                        k++;
                    } else {
                        cout << "remove invalid daughter: (" << id1 << ", " << id2 << ")\n";
                        daughter_ids.erase(daughter_ids.begin() + k);
                    }
                }

                cout << "final_photons size : " << final_photon.size() << "\n";
                double ratio = (double(daughter_ids.size()) / double(final_photon.size()));
                y_values.push_back(ratio);
                daughter_ids.clear();
                event_number++;
            }

        }
        //------------------------ media on the oversamples ----------------------------
        double sum = 0;
        for(const auto& par : y_values) {
            sum += par;     
        }
        double average = sum / y_values.size();
        cout << "the ratio is " << average;
        ratios.push_back(average);
         
    }
ofstream fout("ratios_channel2.txt");
for (auto r : ratios) {
    fout << r << "\n";
}
fout.close();
return(-1);
} 
    