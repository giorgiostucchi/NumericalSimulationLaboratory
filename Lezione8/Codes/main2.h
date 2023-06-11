//Random numbers
#include "random.h"
#include "WaveFunctionBase.h"
#include <vector>
#include <algorithm>
using namespace std;

int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iint;
double walker[m_props];
double mu,sigma;
double x;
double beta;

// averages
double blk_av[m_props], blk_norm;
double glob_av[m_props], glob_av2[m_props];
double stima_int;
double err_int; 

// simulation
int nstep, nblk, npoints;
double delta, delta_annealing;
int accepted, attempted;
int accepted_met = 0, attempted_met = 0;
double rms = 1.;

//results
vector<double> energies;
vector<double> error_energies;
vector<double> mus;
vector<double> sigmas;



//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void OptimizedEnergyRun(void);
void CleanPreviousData(void);
double IntegralCalculation(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void MetropolisIntegration(void);
void MetropolisParametersUpdate(void);
double Error(double,double,int);
void PrintResults(void);
void FillHisto(void);
WaveFunction f;