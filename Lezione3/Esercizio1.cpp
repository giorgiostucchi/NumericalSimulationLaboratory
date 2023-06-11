#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cfloat>
#include "random.h"

using namespace std;

int main()
{

    //*********EXERCISE 1.3************   Direct sampling vanilla option calculation
    Random rnd;
    rnd.SetSeed("Primes", "seed.in"); // loads the seed data from file

    int M = 100000;     // number of throws
    int N = 100;        // number of blocks
    int L = int(M / N); // throws per block
    double S_0 = 100.;  // asset price at start
    double T = 1.;      // delivery time
    double strike_price = 100.;
    double r = 0.1;
    double volatility = 0.25;

    vector<double> av;         // stores the averages computed on each block
    vector<double> av2;        // stores the squared average of each block
    vector<double> av_put;     // stores the averages computed on each block
    vector<double> av_put2;    // stores the squared average of each block
    double call_averages = 0;  // needed to calculate the estimate of the average as the number of blocks increases
    double call_averages2 = 0; // needed to calculate the estimate of the uncertainty on the average as the number of blocks increases
    double put_averages = 0;   // needed to calculate the estimate of the average as the number of blocks increases
    double put_averages2 = 0;  // needed to calculate the estimate of the uncertainty on the average as the number of blocks increases

    string filename = "call_prices.txt";
    string filename2 = "put_prices.txt";
    fstream f, g;
    g.open(filename2, ios::out);
    f.open(filename, ios::out);

    for (int i = 1; i <= N; i++) // cycle over blocks
    {
        double call_sum = 0.;
        double put_sum = 0.;
        for (int j = 0; j < L; j++) // cycle within block
        {
            double S_T = S_0 * exp((r - pow(volatility, 2.) / 2.) * T + volatility * rnd.Gauss(0., 1.) * sqrt(T));
            double call = max(0., S_T - strike_price) * exp(-r * T);
            call_sum += call;
            double put = max(0., strike_price - S_T) * exp(-r * T);
            put_sum += put;
        }
        av.push_back(call_sum / L);                    // average on the block
        av2.push_back(call_sum * call_sum / (L * L));  // average on the block squared
        call_averages = call_averages * double(i - 1); // remove the division of the last cycle (not needed for the first cycle of course)
        call_averages += av.at(i - 1);                 // add average on last block
        call_averages = call_averages / double(i);     // find improved estimate of average
        call_averages2 = call_averages2 * double(i - 1);
        call_averages2 += av2.at(i - 1);
        call_averages2 = call_averages2 / double(i);

        av_put.push_back(put_sum / L);                  // average on the block
        av_put2.push_back(put_sum * put_sum / (L * L)); // average on the block squared
        put_averages = put_averages * double(i - 1);    // remove the division of the last cycle (not needed for the first cycle of course)
        put_averages += av_put.at(i - 1);               // add average on last block
        put_averages = put_averages / double(i);        // find improved estimate of average
        put_averages2 = put_averages2 * double(i - 1);
        put_averages2 += av_put2.at(i - 1);
        put_averages2 = put_averages2 / double(i);
        double sigma_mean = 0;
        double sigma_mean_put = 0;
        if ((i - 1) != 0) // calculate the uncertainty on the mean when N is not 0
        {
            sigma_mean = sqrt((call_averages2 - pow(call_averages, 2)) / (i - 1));
            sigma_mean_put = sqrt((put_averages2 - pow(put_averages, 2)) / (i - 1));
        }
        f << call_averages << " " << sigma_mean << endl;    // save data on output file
        g << put_averages << " " << sigma_mean_put << endl; // save data on output file
    }

    f.close();
    g.close();

    rnd.SaveSeed();
    return 0;
}