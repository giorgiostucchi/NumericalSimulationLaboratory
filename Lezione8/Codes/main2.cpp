#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "main2.h"
#include "WaveFunctionBase.h"
#include "vectoralgebra.h"
#include <cstdlib>
#include <algorithm>

using namespace std;

int main()
{
    CleanPreviousData(); // Clean data from previous executions of the program
    Input(); // Set up the parameters of the simulations (npoints, nblocks, etc..)
    double betamin = 10.;
    delta_annealing = 0.05;
    for(int cycle = 0; cycle < 30; cycle++)
    {
        beta = betamin + cycle * 100.; //update temperature
        cout << cycle + 1 << endl;
        MetropolisParametersUpdate(); // Perform metropolis moves on the parameters sigma and mu (50 proposed steps)
    }
    OptimizedEnergyRun(); //run with optimized parameters (program main.exe)
    FillHisto(); // Fill the histogram of the distribution
    PrintResults(); // Print All of the remaining Results
    return 0;
}

void OptimizedEnergyRun(void){
    std::string command = "./main.exe " + std::to_string(sigmas.back()) + " " + std::to_string(mus.back());
  int result = system(command.c_str());

  // Check if the program ran successfully
  if (result != 0)
  {
    std::cerr << "Error running main.exe" << std::endl;
  }
}

void CleanPreviousData(void){
    int ret = system("make clean_dat");
    if (ret != 0) {
        std::cerr << "Error: Makefile command failed with error code " << ret << std::endl;
    }
}

void FillHisto(void)
{   
    ofstream out; 
    out.open("histo.dat");
    f.Setmu(mus.back());
    f.Setsigma(sigmas.back());
    int max = 50000;
    double y = 1.;
    double y_old = 0.;

    // Metropolis algorithm
    for (int i = 0; i < max; i++)
    {
        y = y_old + rnd.Rannyu(-delta, delta);
        double Prob = f.norm(y) / f.norm(y_old);
        if (rnd.Rannyu() < Prob)
        {
            y_old = y;
        }
        out << y_old << endl;
    }
    out.close();
}

void MetropolisParametersUpdate(void)
{
    attempted_met = 0;
    accepted_met = 0;

    double sigma_old = sigma;
    double mu_old = mu;
    double H_old = IntegralCalculation();
    cout << "Initial mu = " << mu << " initial sigma = " << sigma << endl;

    int annealing_steps = 50;
    for (int i = 0; i < annealing_steps; i++)
    {
        sigma = sigma_old + rnd.Rannyu(-delta_annealing, delta_annealing);
        mu = mu_old + rnd.Rannyu(-delta_annealing, delta_annealing);
        f.Setmu(mu);
        f.Setsigma(sigma);
        double H_new = IntegralCalculation();

        double Prob = exp(-beta * (H_new - H_old));

        if (rnd.Rannyu() < Prob)
        {
            mu_old = mu;
            sigma_old = sigma;
            H_old = H_new;
            accepted_met++;
        }
        else
        {
            sigma = sigma_old;
            mu = mu_old;
            f.Setmu(mu_old);
            f.Setsigma(sigma_old);
        }
        attempted_met++;
    }

    energies.push_back(H_old);
    error_energies.push_back(err_int);
    mus.push_back(mu_old);
    sigmas.push_back(sigma_old);

    cout << "Final Energy = " << H_old << " with error " << err_int << endl;
    cout << "Acceptance rate: " << (double)accepted_met / (double)attempted_met << endl;
    cout << "Scostamento massimo: " << (double)accepted_met * delta_annealing << endl;
    cout << "Values for beta = " << beta << ": final mu =" << mu_old << " and final sigma = " << sigma_old << endl
         << endl;

    return;
}

void PrintResults(void){
    ofstream Ene, Param;

    Ene.open("output_energies_annealing.dat", ios::app);
    Param.open("output_mus_sigmas.dat", ios::app);
    for (unsigned int i = 0; i < energies.size(); i++)
    {
        Ene << energies.at(i) << " " << error_energies.at(i) << endl;
        Param << sigmas.at(i) << " " << mus.at(i) << endl;
    }
    Ene.close();
    Param.close();
}

double IntegralCalculation(void)
{
    for (int iblk = 1; iblk <= nblk; iblk++) // Simulation
    {
        Reset(iblk); // Reset block averages
        for (int istep = 1; istep <= nstep; istep++)
        {
            MetropolisIntegration();
            Accumulate(); // Update block averages
        }
        Averages(iblk);
    }
    return glob_av[iint] / (double)nblk;
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
    npoints = 1000; // numero punti per ogni stima dell'integrale
    mu = 0.8;
    sigma = 0.8;
    nblk = 20;   // numero di blocchi
    nstep = 100; // numero di stime dell'integrale per blocchi

    cout << "The program performs Metropolis Integration" << endl;
    cout << "Moves parameter = " << delta << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl
         << endl;

    f.Setmu(mu);
    f.Setsigma(sigma);

    // Prepare arrays for measurements
    iint = 0; // integral value

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

    stima_int = blk_av[iint] / blk_norm;
    glob_av[iint] += stima_int;
    glob_av2[iint] += stima_int * stima_int;
    err_int = Error(glob_av[iint], glob_av2[iint], iblk);

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