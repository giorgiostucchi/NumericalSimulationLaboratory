#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

int main()
{

    //*********Exercise 1.2.2************
    Random rnd;
    rnd.SetSeed("Primes", "seed.in"); // loads the seed data from file

    int M = 10000;                                                             // number of realisations
    vector<int> reps{1, 2, 10, 100};                                           // different numbers of repetitions
    vector<string> distributions{"standard", "exponential", "Cauchy Lorentz"}; // this vector will be used to cycle over the distributions without copying the code

    for (unsigned int k = 0; k < distributions.size(); k++)
    { // cycle over distribution type
        cout << "Writing data for " + (distributions.at(k)) + " distribution ..." << endl;
        for (unsigned int j = 0; j < reps.size(); j++)
        { // cycle over requested number of repetitions for each distribution
            cout << "Computing averages for N = " + to_string(reps.at(j)) << endl;
            string filename = to_string(M) + (distributions.at(k)) + "_distr_N=" + to_string(reps.at(j)) + ".txt";
            fstream f;
            f.open(filename, ios::out);
            for (int index = 0; index < M; index++)
            { // cycle for M number of realisations of each repetition
                double prog_sum = 0;
                for (int i = 0; i < reps.at(j); i++)
                { // cycle over repetition, different for each distribution
                    if (distributions.at(k) == "standard")
                    {
                        prog_sum += rnd.Gauss(0., 1.);
                    }
                    else if (distributions.at(k) == "exponential")
                    {
                        prog_sum += rnd.Expo(1.);
                    }
                    else if (distributions.at(k) == "Cauchy Lorentz")
                    {
                        prog_sum += rnd.Cauchy_Lorentz(1., 0.);
                    }
                    else
                    {
                        cerr << "PROBLEM!" << endl;
                    }
                }
                f << prog_sum / reps.at(j) << endl;
            }
            f.close();
        }
    }

    cout << "Done!" << endl;

    rnd.SaveSeed();
    return 0;
}