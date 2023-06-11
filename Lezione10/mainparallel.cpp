#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "classes.h"
#include "vectoralgebra.h"
#include <cstdlib>
#include <stdio.h>
#include <mpi.h>
#include <algorithm>

using namespace std;

//The program performs a minimization by means of a genetic algorithm. in particular we minimize the distance
//cities generated randomly on a map or uploaded from file. we do simulations that exchange their best genes every totalsteps/10 steps. we specify that they are communicating processes 
//with the "int program = 0", which tells the classes that we are using this specific program. The classe are in common. TO COMPILE: DEACTIVATE CONDA.
//TO RUN: mpiexec -np 10 ./mainparallel.exe 0 100   (for example)

int main(int argc, char **argv)
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int program = 0; // tells the functions in the classes that we're using the current program (parallel)
    int opt;         // which map we want to use
    int steps;       //how many evolution steps
    vector<int> receiver_sequence(size, 0); //needed to organize the communication

    if (rank == 0)
    {
        if (argc != 3)
        {
            cout << "first argument: <0> circle map, <1> map on square, <2> American Capitals; second argument: evolution steps " << endl;
            return -1;
        }

        opt = atoi(argv[1]);
        steps = atoi(argv[2]);

        CleanPreviousData(opt, program);
    }
    MPI_Bcast(&opt, 1, MPI_INT, 0, MPI_COMM_WORLD); //communicate map option and number of steps to all subprocesses
    MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); 

    int n_ind = 1000;
    int n_cities = 34;
    int n_migr = steps / 10; //how many exchanges during the full evolution

    if (rank == 0)
    {
        cout << "The program uses Genetic Algorithms to solve the travelling salesman problem" << endl;
        cout << size << " populations are created and evolve independently, except for a random exchange of their best individuals every " << n_migr << " generations." << endl;
        cout << "Each population has " << steps << " generations." << endl
             << endl;
    }

    Map myMap(n_cities, opt, 0);
    MPI_Barrier(MPI_COMM_WORLD);
    Population myPop(n_ind, myMap, rank + 1); //by including the rank in the generation process, we create different populations for each process

    Random rnd; //the only one using this is actually rank 0 later. however, if I included this line in a if(rank==0) statement, when I regenerate the receiver sequence it would always end with the same number, because the random number generator doesn't use the "primes" in the first extraction
    rnd.SetSeed("Primes", "seed.in");

    myPop.CheckBounds();
    myPop.SortByFitness();

    for (int k = 0; k < (steps / n_migr); k++)
    {
        if (rank == 0)
        {
            iota(receiver_sequence.begin(), receiver_sequence.end(), 0);
            rnd_shuffle_complete(receiver_sequence.begin(), receiver_sequence.end(), rnd); //generate a random receiver sequence
        }
        MPI_Bcast(receiver_sequence.data(), receiver_sequence.size(), MPI_INT, 0, MPI_COMM_WORLD); //send the sequence to everyone
        MPI_Barrier(MPI_COMM_WORLD);//wait for all

        myPop.FullEvolution(n_migr, rank); //evolve the process

        MPI_Barrier(MPI_COMM_WORLD); //wait for all before exchanging individuals with best genes

        vector<int> best_genes = myPop.GetIndividual(0).GetGenes(); //save best genes                                                                             
        vector<int> buf(best_genes.size()); //create receiver memory allocation
        int tag1 = 1;

        for(int i = 0; i <  size; i++){
            if(rank == i) MPI_Send(best_genes.data(), best_genes.size(), MPI_INT, receiver_sequence.at(i), tag1, MPI_COMM_WORLD); //send genes
            if(rank == receiver_sequence.at(i)) MPI_Recv(buf.data(), best_genes.size(), MPI_INT, i, tag1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //receive genes
        }                                                                                   
        myPop.GetIndividual(0).SetGenes(buf); //Set new genes for best individual (zero-th) of the population of the process 
    }

    if (rank == 0)
    {
        myPop.PrintMap(); //print map, only once
    }

    myPop.PrintBestPath(rank); //each process print its optimized path

    MPI_Finalize();

    return 0;
}