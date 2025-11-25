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
#include "utils.h"

using namespace std;


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
        int event_number = 0;
        vector<vector<int>> active_daughter_ids;            // IDs of daughters still "alive" for detection

        while (getline(fin, line)) {

            if (starts_with(line, "# event")){
                cout << line << "\n";
                cout << "looking for " << event_number << "\n";
            }

            if (starts_with(line, "# event " + to_string(event_number) + " in")){ 
                cout << "event number : " << event_number << "\n";
            }

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

                if (partspdg[0] == 313 || partspdg[0] == -313 || partspdg[0] == 323 || partspdg[0] == -323) {  //only K*(892): more abundant that decays into Kπ
                    // daughter
                    int d1 = partsid[1];
                    int d2 = partsid[2];

                    // tracking daughters: map and vector of active IDs
                    active_daughter_ids.push_back({d1, d2});
                } 
            }   


                // ---------------------------------  large block to count the particles that reached the detectors  --------------------------------------
            
             if (starts_with(line, "# event " + to_string(event_number) + " out")){
                vector <int> final_Kminus;             //IDs of K- that reach the detectors
                cout << "found the end of the event \n";
                for (int i=0; i < parse_Kminus_final(line); ){
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
                ratios[event_number] = ratio;
                active_daughter_ids.clear();
                 event_number++;
             }

        }
    }

    //------------------------ media on the oversamples ----------------------------
    double sum = 0;
    for(const auto& par : ratios) {
        sum += par.second;     
    }
    double average = sum / ratios.size();
    cout << "a razão é " << average;
    return(average);
} 
    