//The purpose of this code is obtain the multiplicity of the events 
//Here we consider just charged particles with -0.5< η < 0.5

#include <iostream>
#include <vector>
#include <cmath>
#include "utils.h"

using namespace std;

//calculation of eta
double eta(double px, double py, double pz)
{
    double p = sqrt(px*px + py*py + pz*pz);
    return 0.5 * log( (p + pz) / (p - pz) );
}


int main() {
    ifstream list("test.txt");
    if (!list.is_open()) {
        cerr << "Error to open files list.\n";
    }

    string filename;
    vector<int> particles_per_event;
    while (getline(list, filename)) {
        if (filename.empty()) continue;

        ifstream fin(filename);
        if (!fin.is_open()) {
            cerr << "Error to open: " << filename << "\n";
            continue;
        }
        string line;
        int charged_particles = 0;
        int event_number = 0 ;
        while (getline(fin, line)) {

            if (starts_with(line, "# event " + to_string(event_number) + " out")){
                event_number++;
                cout << "found the end of the event \n";
                cout << "there are "<< parse_Kminus_final(line) << " particles \n";
                for (int i=0; i < parse_Kminus_final(line); i++){
                    string pline;
                    if (!getline(fin, pline)) break;
                    ParticleInfo p;
                    if (!parse_particle_row(pline, p)) continue;
                    double et = eta(p.px, p.py, p.pz);
                    if (p.charge != 0){ // && et > -0.5 && et < 0.5
                        charged_particles++;
                    }
                }
                particles_per_event.push_back(charged_particles);
                charged_particles = 0;
            }
        }
    }
    //------------------------ media on the oversamples ----------------------------
    double sum = 0;
    for(const auto& par : particles_per_event) {
        sum += par;     
    }
    double average = sum / particles_per_event.size();
    cout << "the multiplicity is " << average;
    return(average);
}