#include <cmath>
#include <iostream>
#include <iomanip>
#include "Integrals_complete.h"
#include "FunzioneBase.h"
#include <vector>
#include "random.h"
#include <fstream>

using namespace std;

//This program performs the integration with uniform distribution sampling

int main()
{

  Integranda f; // Instanciate Class for Function to be integrated
  MeanMC Integratore(0., 1., "Primes", "seed.in");

  int M = 10000;      // number of throws
  int N = 100;        // number of blocks
  int L = int(M / N); // throws per block

  vector<double> av;        // stores the averages computed on each block
  vector<double> av2;       // stores the squared average of each block
  double sum_averages = 0;  // needed to calculate the estimate of the average as the number of blocks increases
  double sum_averages2 = 0; // needed to calculate the estimate of the uncertainty on the average as the number of blocks increases

  string filename = "integrals&uncertainties.txt";
  ofstream g;
  g.open(filename, ios::out);

  for (int i = 1; i <= N; i++)
  {
    double sum_block = 0;
    for (int j = 0; j < L; j++)
    {
      sum_block += Integratore.Integra(1000, f); // progressive sum within each block
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
    g << sum_averages << " " << sigma_mean << endl; // save data on output file
  }

  return 0;
}