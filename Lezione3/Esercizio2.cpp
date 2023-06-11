#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

int main()
{

    //*********EXERCISE 2*********** Sampling of discretized GBM
    Random rnd;
    rnd.SetSeed("Primes", "seed.in"); // loads the seed data from file

    int M = 1000000;    // number of throws
    int N = 100;        // number of blocks
    int L = int(M / N); // throws per block

    vector<string> exercises{"call", "put"}; // this vector will be used to cycle over the exercises without copying the code

    for (unsigned int k = 0; k < exercises.size(); k++)
    {
        cout << "Writing " + (exercises.at(k)) + " option prices ..." << flush;

        vector<double> av;        // stores the averages computed on each block
        vector<double> av2;       // stores the squared average of each block
        double sum_averages = 0;  // needed to calculate the estimate of the average as the number of blocks increases
        double sum_averages2 = 0; // needed to calculate the estimate of the uncertainty on the average as the number of blocks increases
        double S_0 = 100.;        // asset price at start
        double T = 1.;            // delivery time
        int nstep = 100;
        double dt = T / nstep;
        double strike_price = 100.;
        double r = 0.1;
        double volatility = 0.25;

        string filename = (exercises.at(k)) + "_option_prices.txt";
        fstream f;
        f.open(filename, ios::out);

        for (int i = 1; i <= N; i++) // cycle over blocks
        {
            double sum_block = 0;
            for (int j = 0; j < L; j++) // cycle within the block
            {
                double S_T = S_0;
                for (int index = 0; index < nstep; index++) // compute the S_T in steps of dt
                {
                    S_T *= exp((r - pow(volatility, 2.) / 2.) * dt + volatility * rnd.Gauss(0., 1.) * sqrt(dt));
                }
                if ((exercises.at(k) == "put"))
                {
                    sum_block += max(0., strike_price - S_T) * exp(-r * T); // progressive sum within each block
                }
                else if ((exercises.at(k) == "call"))
                {
                    sum_block += max(0., S_T - strike_price) * exp(-r * T); // progressive sum within each block
                }
                else
                {
                    cerr << "PROBLEM!" << endl;
                }
            }
            av.push_back(sum_block / double(L));            // store the average on the block
            av2.push_back(sum_block * sum_block / (L * L)); // store the average on the block squared
            sum_averages = sum_averages * double(i - 1);    // remove the division of the last cycle (not needed for the first cycle of course)
            sum_averages += av.at(i - 1);                   // add average on last block
            sum_averages = sum_averages / double(i);        // find improved estimate of average
            sum_averages2 = sum_averages2 * double(i - 1);  // do the same for the squared average
            sum_averages2 += av2.at(i - 1);
            sum_averages2 = sum_averages2 / double(i);
            double sigma_mean = 0;
            if ((i - 1) != 0) // calculate the uncertainty on the mean when N is not 0
            {
                sigma_mean = sqrt((sum_averages2 - pow(sum_averages, 2)) / (i - 1));
            }
            f << sum_averages << " " << sigma_mean << endl; // save data on output file
        }

        f.close();
        cout << "  Done!" << endl;
    }

    rnd.SaveSeed(); // if we want to run the executable for more than one time, we'll have to substitute the seed.in file with the contents of the seed.out file
    return 0;
}