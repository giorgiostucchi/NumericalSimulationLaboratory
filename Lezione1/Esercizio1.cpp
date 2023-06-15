#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

double GetChiSquared(vector<double> bins, double n, double m)
{
   double chi_squared = 0;
   for (int i = 0; i < m; i++)
   {
      chi_squared += pow((double(bins.at(i)) - n / m), 2) / (n / m);
   }
   return chi_squared;
};

int main()
{

   //*********EXERCISES 1.1.1 & 1.1.2************
   Random rnd;
   rnd.SetSeed("Primes", "seed.in"); // loads the seed data from file

   int M = 100000;     // number of throws
   int N = 100;        // number of blocks
   int L = int(M / N); // throws per block

   vector<int> exercises{1, 2}; // this vector will be used to cycle over the exercises without copying the code

   for (unsigned int k = 0; k < exercises.size(); k++)
   {
      cout << "Writing data for exercise 1.1." << to_string(exercises.at(k)) << " ...";

      vector<double> av;        // stores the averages computed on each block
      vector<double> av2;       // stores the squared average of each block
      double sum_averages = 0;  // needed to calculate the estimate of the average as the number of blocks increases
      double sum_averages2 = 0; // needed to calculate the estimate of the uncertainty on the average as the number of blocks increases

      string filename = "data_exercise_1.1." + to_string(exercises.at(k)) + ".txt";
      fstream f;
      f.open(filename, ios::out);

      for (int i = 1; i <= N; i++)
      {
         double sum_block = 0;
         for (int j = 0; j < L; j++)
         {
            if ((exercises.at(k) == 1))
            {
               sum_block += rnd.Rannyu(); // progressive sum within each block for exercise 1.1
            }
            else if ((exercises.at(k) == 2))
            {
               sum_block += pow((rnd.Rannyu() - 0.5), 2); // progressive sum within each block for exercise 1.2
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

   //**********EXERCISE 1.3*************
   double n = 1E4; // # throws each time
   double m = 1E2; // # bins in the interval [0,1]

   cout << "Calculating chi_squared for exercise 1.1.3 ...";

   fstream f;
   f.open("data_100_chi_squared.txt", ios::out);

   for (int i = 0; i < 100; i++) // calculate chi squared 100 times
   {
      vector<double> bins(100, 0); // inizialize a vector of size 100 with zeros
      for (int j = 0; j < n; j++)
      {
         double extraction = rnd.Rannyu() * 100.; // fill the bins by extracting within the [0,1] range and multiplying by 100
         bins.at(floor(extraction))++;
      }
      f << GetChiSquared(bins, n, m) << endl; // compute and save chi squared on file
   }

   f.close();
   cout << "  Done!" << endl;

   rnd.SaveSeed(); // if we want to run the executable for more than one time, we'll have to substitute the seed.in file with the contents of the seed.out file
   return 0;
}