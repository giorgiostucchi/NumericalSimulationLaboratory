#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "classes.h"
#include "vectoralgebra.h"
#include <cstdlib>
#include <algorithm>

using namespace std;

//The program performs a minimization by means of a genetic algorithm. in particular we minimize the distance
//between cities generated randomly on a map.

int main(int argc, char **argv)
{

    if (argc != 3)
    {
        cout << "first argument: <0> circle map, <1> map on square; second argument: evolution steps " << endl;
        return -1;
    }

    int opt = atoi(argv[1]);
    int steps = atoi(argv[2]);

    CleanPreviousData(opt); //Clean data from previous runs of the program (for the selected map)

    int n_ind = 1000; //number of individuals
    int n_cities = 34; //number of cities / genes

    Map myMap(n_cities, opt); //Istanciate the random map
    Population myPop(n_ind, myMap); //create the population

    myPop.CheckBounds(); //check that the population abides by the rules we want (starts with 1, contains integers in range 1-ncities)
    myPop.SortByFitness(); //order the population by cost function performance

    myPop.FullEvolution(steps); //evolve the population

    myMap.PrintMap();//print map
    myPop.PrintBestPath();//print final optimized path

    return 0;
}