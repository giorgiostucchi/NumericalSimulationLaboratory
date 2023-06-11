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

    //*********EXERCISE 1.3************ 
    //This code estimates the value of pi using the Buffon's needle experiment and calculates the average
    // and uncertainty as the number of blocks increases
    Random rnd;
    rnd.SetSeed("Primes", "seed.in"); // loads the seed data from file

    int M = 100000;     // number of throws
    int N = 100;        // number of blocks
    int L = int(M / N); // throws per block
    double l = 0.5;     // length of needle
    double d = 1.0;     // distance between lines

    vector<double> av;        // stores the averages computed on each block
    vector<double> av2;       // stores the squared average of each block
    double sum_averages = 0;  // needed to calculate the estimate of the average as the number of blocks increases
    double sum_averages2 = 0; // needed to calculate the estimate of the uncertainty on the average as the number of blocks increases

    string filename = "pi_averages&uncertainties.txt";
    fstream f;
    f.open(filename, ios::out);

    for (int i = 1; i <= N; i++)
    {
        int hits = 0; // success counter
        for (int j = 0; j < L; j++)
        {
            double end_height = rnd.Rannyu() * d;  // extract coordinate of one endpoint of segment
            double angle = rnd.Rannyu() * DBL_MAX; // extract angle randomly
            double proj = l * sin(angle);          // compute the projection of the segment on the vertical component
            if (floor(end_height / d) != floor((end_height + proj) / d))
            { // check if the segment endpoints lie on different sides of a line
                hits++;
            }
        }
        double pi_estimate = (2 * L * l) / (double(hits) * d);
        cout << pi_estimate << endl;
        av.push_back(pi_estimate);                   // average on the block
        av2.push_back(pi_estimate * pi_estimate);    // average on the block squared
        sum_averages = sum_averages * double(i - 1); // remove the division of the last cycle (not needed for the first cycle of course)
        sum_averages += av.at(i - 1);                // add average on last block
        sum_averages = sum_averages / double(i);     // find improved estimate of average
        sum_averages2 = sum_averages2 * double(i - 1);
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

    rnd.SaveSeed();
    return 0;
}