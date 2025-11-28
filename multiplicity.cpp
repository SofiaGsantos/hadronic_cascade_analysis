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

    vector<double> multiplicity;
    string filename;

    while (getline(list, filename)) {
        if (filename.empty()) continue;

        ifstream fin(filename);
        if (!fin.is_open()) {
            cerr << "Error to open: " << filename << "\n";
            continue;
        }
        cout << "open file " << filename << "\n";

        string line;
        int event_number = 0;
        vector<double> x_values;

        while (getline(fin, line)) {

            if (starts_with(line, "# event " + to_string(event_number) + " out")){
                
                cout << "found the end of an event \n";
                int nch = 0;

                for (int i=0; i < parse_Kminus_final(line); i++){
                    string pline;
                    if (!getline(fin, pline)) break;
                    ParticleInfo p;
                    if (!parse_particle_row(pline, p)) continue;
                    if (p.charge == 0) continue;
                    double et = eta(p.px, p.py, p.pz);
                    if ( et > -0.5 && et < 0.5){
                        nch++;
                    }
                }
                double x = pow(nch, 1.0/3.0);
                x_values.push_back(x);

                event_number++;
            }
        }
        //------------------------ media on the oversamples ----------------------------
        double sum = 0;
        for(const auto& par : x_values) {
            sum += par;     
        }
        double average = sum / x_values.size();
        cout << "the multiplicity is " << average;
        multiplicity.push_back(average);
    }
ofstream fout("multiplicity.txt");
for (auto m : multiplicity) {
    fout << m << "\n";
}
fout.close();
return -1;
}