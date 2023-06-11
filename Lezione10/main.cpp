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
//cities generated randomly on a map. we do 10 (this can be changed) independent simulations. we specify that they are independent 
//with the "int program = 1", which tells the classes that we are using this specific program and not
// the parallelized one. The classe are in common.

int main(int argc, char **argv)
{

    if (argc != 4)
    {
        cout << "first argument: <0> circle map, <1> map on square, <2> AmericanCapitals; second argument: evolution steps ; third argument: number of independent process" << endl;
        return -1;
    }

    int program = 1;

    int opt = atoi(argv[1]);
    int steps = atoi(argv[2]);
    int processes = atoi(argv[3]);

    CleanPreviousData(opt, program); //Clean data from previous runs of the program (for the selected map)
    
    int n_ind = 1000; //number of individuals
    int n_cities = 34; //number of cities / genes

    Map myMap(n_cities, opt, program);  //Istanciate the random map
    myMap.PrintMap();

    for(int i = 0; i < processes; i++){
        Population myPop(n_ind, myMap, i); //create the population

        myPop.CheckBounds(); //check that the population abides by the rules we want (starts with 1, contains integers in range 1-ncities)
        myPop.SortByFitness(); //order the population by cost function performance

        myPop.FullEvolution(steps, i); //evolve the population

        myPop.PrintBestPath(i);//print final optimized path
    }

    return 0;
}