/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  int counter = 0;
  Input(); //Initialize the system: set the starting temp, number of spins, etc reading from input.dat
  double T0 = 0.5; //initial temp
  double T1 = 2.0; //final temp
  int n_d = 15; //number of intermediate temperatures
  double d_T = (T1 - T0) / double(n_d);
  for (int i = 0; i <= n_d; i++)
  {
    double T = T0 + i * d_T;
    Update(T); //set new temperature and load final configuration of previous step
    Equilibration(1500, counter); //equilibrate the system so that is reaches a configuration corresponding to the new temperature
    //the equilibration steps were included for safety, but they are not stricly necessary in this case, see notebook for explanation
    Simulation(); //perform metropolis algorithm and makes measurements on the system
  }                                        
  return 0;
}

void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; i++)
  {
    s[i] = -1;
}
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;

ConfFinal();
}

void Simulation()
{
for (int iblk = 1; iblk <= nblk; ++iblk) // Simulation
{
    Reset(iblk); // Reset block averages
    for (int istep = 1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); // Update block averages
    }
    Averages(iblk); // Print results for current block
}
ConfFinal(); // Write final configuration
}

void Equilibration(int eq_step, int &counter)
{ 
  string option;
  if(metro==1) option = "metropolis";
  else option = "gibbs";
  ofstream inst_ene;
  ofstream inst_mag;
  if(counter == 0){
    inst_ene.open("inst_ene." + option + ".dat", ios::out); // Open files outside the loop
    inst_mag.open("inst_mag." + option + ".dat", ios::out);
  }else{
  inst_ene.open("inst_ene." + option + ".dat", ios::app); // Open files outside the loop
  inst_mag.open("inst_mag." + option + ".dat", ios::app);
  }

  for (int i = 1; i <= eq_step; i++) // Simulation
  {
    Move(metro);
    Measure();
    if(counter == 0){
    inst_ene << i << " " << walker[iu] / double(nspin) << endl;
    inst_mag << i << " " << walker[im] / double(nspin) << endl;
    }
  }
  counter++;
  inst_ene.close(); // Close files after the loop
  inst_mag.close();
  ConfFinal(); // Write final configuration
}

void Update(double T)
{
  temp = T;
  beta = 1.0 / temp;
  cout << "Temperature = " << temp << endl;

    ifstream ReadConfig("config.final");
    if (!ReadConfig.good())
    {
      cerr << "error in opening config.final" << endl;
      exit(-1);
    };
    for (int i = 0; i < nspin; ++i)
    {
      ReadConfig >> s[i];
    }
}

void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for (int i = 0; i < nspin; ++i)
  {
    // Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu() * nspin);
    sm = s[o];
    energy_old = Boltzmann(sm, o);
    energy_new = Boltzmann(-sm, o);
    double delta_energy = -energy_old + energy_new;
    attempted++;

    if (metro == 1) // Metropolis
    {
      double a = min(1., exp(-beta * delta_energy));
      double r = rnd.Rannyu();
      if (r < a)
      {
        s[o] = -sm;
        accepted++;
      }
    }
    else // Gibbs sampling
    {
      accepted++;                                      // Gibbs sampling has a=1
      double p = 1. / (1. + exp(beta * delta_energy)); // flipping probability
      double r = rnd.Rannyu();
      if (r < p)
        s[o] = -sm;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * (s[Pbc(ip - 1)] + s[Pbc(ip + 1)]) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

  // cycle over spins
  for (int i = 0; i < nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i + 1)] - 0.5 * h * (s[i] + s[Pbc(i + 1)]);
    m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
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

  ofstream Ene, Heat, Mag, Chi;
  const int wd = 12;
  string option;
  if(metro==1) option = "metropolis";
  else option = "gibbs";

  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted / attempted << endl << endl;

  Ene.open("output.ene.0." + option, ios::app);
  Heat.open("output.heat.0." + option, ios::app);
  Mag.open("output.mag.0." + option, ios::app);
  Chi.open("output.chi.0." + option, ios::app);
  stima_u = blk_av[iu] / blk_norm / (double)nspin;                           //energy
  stima_c = beta*beta * (blk_av[ic]/blk_norm - stima_u*stima_u * (double)nspin * (double)nspin) / (double)nspin;            //heat capacity
  stima_m = blk_av[im]/blk_norm;                                           //magnetization
  stima_x = beta * blk_av[ix]/blk_norm;                                   //susceptibility
  glob_av[iu] += stima_u;
  glob_av2[iu] += stima_u * stima_u;
  glob_av[ic] += stima_c;
  glob_av2[ic] += stima_c * stima_c;
  glob_av[im] += stima_m;
  glob_av2[im] += stima_m * stima_m;
  glob_av[ix] += stima_x;
  glob_av2[ix] += stima_x * stima_x;
  err_u = Error(glob_av[iu], glob_av2[iu], iblk);
  err_c = Error(glob_av[ic], glob_av2[ic], iblk);
  err_m = Error(glob_av[im], glob_av2[im], iblk);
  err_x = Error(glob_av[ix], glob_av2[ix], iblk);
  Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double)iblk << setw(wd) << err_u << endl;
  Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic] / (double)iblk << setw(wd) << err_c << endl;
  if(h!=0) Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im] / (double)iblk << setw(wd) << err_m << endl;
  Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix] / (double)iblk << setw(wd) << err_x << endl;
  Ene.close();
  Heat.close();
  Mag.close();
  Chi.close();


  // INCLUDE YOUR CODE HERE

  cout << "----------------------------" << endl
       << endl;
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl
       << endl;
  WriteConf.open("config.final");
  for (int i = 0; i < nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i) // Algorithm for periodic boundary conditions
{
  if (i >= nspin)
    i = i - nspin;
  else if (i < 0)
    i = i + nspin;
  return i;
}

double Error(double sum, double sum2, int iblk)
{
  if (iblk == 1)
    return 0.0;
  else
    return sqrt((sum2 / (double)iblk - pow(sum / (double)iblk, 2)) / (double)(iblk - 1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
