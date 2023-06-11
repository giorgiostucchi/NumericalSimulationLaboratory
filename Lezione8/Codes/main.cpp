#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "main.h"
#include "WaveFunctionBase.h"
#include <algorithm>

using namespace std;

int main(int argc, char* argv[])
{
    // Check if the correct number of arguments was passed
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <sigma> <mu>" << std::endl;
    return 1;
  }

  // Convert the command-line arguments to double
  try {
    sigma = std::stod(argv[1]);
    mu = std::stod(argv[2]);
  } catch (...) {
    std::cerr << "Invalid argument format" << std::endl;
    return 1;
  }

  Input();
  for (int iblk = 1; iblk <= nblk; iblk++) // Simulation
  {
    Reset(iblk); // Reset block averages
    for (int istep = 1; istep <= nstep; istep++)
    {
      MetropolisIntegration();
      Accumulate(); // Update block averages
    }
    Averages(iblk); // Print results for current block
  }
  return 0;
}

void Input(void)
{
  ifstream ReadInput, Primes, Seed;

  // Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2;
  Primes.close();

  Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  Seed.close();

  delta = 3.;
  npoints = 1000; //numero punti per ogni stima dell'integrale
  nblk = 100; //numero di blocchi
  nstep = 1000; // numero di stime dell'integrale per blocchi

  cout << "The program performs Metropolis Integration" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl
       << endl;

  f.Setmu(mu);
  f.Setsigma(sigma);

  // Prepare arrays for measurements
  iint = 0; //integral slot (just one in this case)

  n_props = 1; // Number of observables

  return;
}

void MetropolisIntegration() // Properties measurement
{
  double x_old = 0;

  // Metropolis algorithm
  for (int i = 0; i < npoints; i++)
  {
    x = x_old + rnd.Rannyu(-delta, delta);
    double Prob = f.norm(x) / f.norm(x_old);
    if (rnd.Rannyu() < Prob)
    {
      x_old = x;
      accepted++;
    }
    walker[iint] += f.hamiltonian(x_old);
    attempted++;
  }
  walker[iint] /= npoints;
  return;
}

void Reset(int iblk) // Reset block averages
{

  if (iblk == 1)
  {
    for (int i = 0; i < n_props; ++i)
    {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for (int i = 0; i < n_props; ++i)
  {
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Accumulate(void) // Update block averages
{

  for (int i = 0; i < n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) // Print results for current block
{

  ofstream Int_file;
  const int wd = 12;

  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << double(accepted) / double(attempted) << endl
       << endl;
  Int_file.open("output_int.dat", ios::app);

  stima_int = blk_av[iint] / blk_norm;
  glob_av[iint] += stima_int;
  glob_av2[iint] += stima_int * stima_int;
  err_int = Error(glob_av[iint], glob_av2[iint], iblk);

  // Potential energy per particle
  Int_file << setw(wd) << iblk << setw(wd) << stima_int << setw(wd) << glob_av[iint] / (double)iblk << setw(wd) << err_int << endl;

  cout << "----------------------------" << endl
       << endl;

  Int_file.close();
}

double Error(double sum, double sum2, int iblk)
{
  if (iblk == 1)
  {
    return 0;
  }
  else
  {
    return sqrt(fabs(sum2 / (double)iblk - pow(sum / (double)iblk, 2)) / (double)(iblk - 1));
  }
}