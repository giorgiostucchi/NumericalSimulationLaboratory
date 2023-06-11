#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "vectoralgebra.h"

using namespace std;

//This Program performs simulation of Random Walks (discrete and continuous)

int main()
{

   //*********EXERCISES 2.2.1 & 2.2.2************
   Random rnd;
   rnd.SetSeed("Primes", "seed.in"); // loads the seed data from file

   int I = 100;
   int M = 10000;      // number of throws
   int N = 100;        // number of blocks
   int L = int(M / N); // throws per block

   vector<int> exercises{1, 2}; // this vector will be used to cycle over the exercises without copying the code

   for (unsigned int k = 0; k < exercises.size(); k++)
   {
      cout << "Writing data for exercise 2.2." << to_string(exercises.at(k)) << " ...";

      vector<vector<double>> av;             // stores the averages computed on each block
      vector<vector<double>> av2;            // stores the squared average of each block
      vector<double> sum_averages(100, 0.);  // needed to calculate the estimate of the average as the number of blocks increases
      vector<double> sum_averages2(100, 0.); // needed to calculate the estimate of the uncertainty on the average as the number of blocks increases
      vector<double> sigma_mean(100, 0.);

      string filename = "data_exercise_2.2." + to_string(exercises.at(k)) + ".txt";
      fstream f;
      f.open(filename, ios::out);

      for (int i = 1; i <= N; i++)
      { // cycle over number of blocks
         vector<double> sum_block(100, 0.);
         for (int j = 0; j < L; j++)
         { // cycle within the block of 100 RWs
            vector<double> position{0., 0., 0.};
            for (int steps = 0; steps < I; steps++)
            { // compute the single random walk (I=100 steps)
               if ((exercises.at(k) == 1))
               {
                  double extract = double(floor(3. * rnd.Rannyu()));             // extract x,y,z
                  double direction = double(floor(2. * rnd.Rannyu())) * 2. - 1.; // extract positive or negative movement
                  position.at(extract) += direction;
                  sum_block.at(steps) += pow((mod(position)), 2.); // compute rsd at step i for single random walk and store it
               }
               else if ((exercises.at(k) == 2))
               {
                  double phi = rnd.Rannyu() * 2. * M_PI;
                  double theta = rnd.Rannyu() * M_PI;
                  position.at(0) += sin(theta) * cos(phi);
                  position.at(1) += sin(theta) * sin(phi);
                  position.at(2) += cos(theta);
                  sum_block.at(steps) += pow((mod(position)), 2.);
               }
               else
               {
                  cerr << "PROBLEM!" << endl;
               }
            }
         }
         av.push_back((sum_block / double(L)) ^ 0.5); // store the averages of the 100 RW in the block for each step
         cout << "sum_block[0] " << (sum_block.at(0)) << endl;
         av2.push_back((sum_block / double(L)));      // averages squared
         sum_averages = sum_averages * double(i - 1); // remove the division of the last cycle (not needed for the first cycle of course)
         sum_averages += av.at(i - 1);                // add average on last block
         sum_averages = sum_averages / double(i);     // find improved estimate of average
         cout << "sum _ averages " << sum_averages.at(0) << endl;
         sum_averages2 = sum_averages2 * double(i - 1); // do the same for the squared average
         sum_averages2 += av2.at(i - 1);
         sum_averages2 = sum_averages2 / double(i);
         if (i == 100)
         {
            sigma_mean = ((sum_averages2 - (sum_averages ^ 2.)) / double((i - 1))) ^ 0.5;
            for (unsigned int it = 0; it < sum_averages.size(); it++)
               f << sum_averages.at(it) << " " << sigma_mean.at(it) << endl;
         }
      }

      f.close();
      cout << "  Done!" << endl;
   }

   rnd.SaveSeed(); // if we want to run the executable for more than one time, we'll have to substitute the seed.in file with the contents of the seed.out file
   return 0;
}