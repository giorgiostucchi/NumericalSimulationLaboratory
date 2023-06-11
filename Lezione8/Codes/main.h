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

// averages
double blk_av[m_props], blk_norm;
double glob_av[m_props], glob_av2[m_props];
double stima_int;
double err_int; 

// simulation
int nstep, nblk, npoints;
double delta;
int accepted, attempted;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void MetropolisIntegration(void);
double Error(double,double,int);
WaveFunction f;