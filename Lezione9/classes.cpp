#include "classes.h"
#include <iostream>
#include <fstream>
#include <algorithm> // std::random_shuffle
#include <vector>    // std::vector
#include <numeric>
#include <cmath>
#include <random>
#include "vectoralgebra.h"
#include <functional>

using namespace std;

void CleanPreviousData(int option){
    string option_string;
    if(option == 0)  option_string = "circumference";
    if(option == 1)  option_string = "square";
    string command = "rm *_" + option_string +".dat";
    int ret = system(command.c_str());
    if (ret != 0) {
        std::cerr << "Error: Makefile command failed with error code " << ret << endl<< endl;
    }
}

void rnd_shuffle(vector<int>::iterator first, vector<int>::iterator last, Random &rnd)
{
    int n = (last - first);
    for (int i = n - 1; i > 0; --i)
    {
        swap(first[i], first[int(rnd.Rannyu() * (i + 1))]);
    }
}

//Map Methods

Map::Map(int ncities, int opt)
{
    m_ncities = ncities;
    m_opt = opt;
    string infile_primes = "Primes";
    string infile_seed = "seed.in";
    m_gen.SetSeed(infile_primes, infile_seed);

    // circle option
    if (opt == 0)
    {
        for (int i = 0; i < ncities; i++)
        {
            double theta = 2 * M_PI * m_gen.Rannyu();
            double x = cos(theta);
            double y = sin(theta);
            m_cities_coord.push_back(make_pair(x, y));
        }
        // cout << m_cities_coord.at(0).first << "  " << m_cities_coord.at(0).second << endl;
        // cout << GetNthCoordinates(0).first << " " << GetNthCoordinates(0).second << endl;
        cout << "cities generated on an unitary circle" << endl;
    }

    // random in square option
    else if (opt == 1)
    {
        for (int i = 0; i < ncities; i++)
        {
            double x = m_gen.Rannyu(-1, 1);
            double y = m_gen.Rannyu(-1, 1);
            m_cities_coord.push_back(make_pair(x, y));
        }
        cout << "cities generated in a square" << endl;
    }

    else
    {
        cout << "map option error" << endl;
        exit(-1);
    };
};

void Map::PrintMap(void){
    string option_string;
    int option = m_opt;
    if(option == 0)  option_string = "circumference";
    if(option == 1)  option_string = "square";
    string filenamebest = "map_" + option_string + ".dat";

    ofstream best;

    best.open(filenamebest, ios::app);
    for(int i = 0; i < m_ncities; i++){
        best << i << " " << GetNthCoordinates(i).first << " " <<  GetNthCoordinates(i).second << endl;
    }
    best.close();
    cout << "Map on " + option_string +" printed " << endl;
}


//Individual Methods

void Individual::Fitness(void)
{
    double sum = 0.;
    for (unsigned int i = 0; i < m_genes.size(); i++)
    {
        sum += NormEuclidean(m_map->GetNthCoordinates(m_genes.at(i) - 1) - m_map->GetNthCoordinates(m_genes.at((i + 1) % m_genes.size()) - 1));
    }
    m_fitness = sum;
};


//Population Methods

Population::Population(int individuals, Map &myMap) : m_map(myMap)
{
    string infile_primes = "Primes";
    string infile_seed = "seed.in";
    m_gen.SetSeed(infile_primes, infile_seed);

    m_ind = individuals;
    m_ngenes = myMap.GetNumberCities();

    for (int i = 0; i < m_ind; i++)
    {
        vector<int> vect(m_ngenes);
        iota(vect.begin(), vect.end(), 1);
        myVectUtils::Print(vect);
        rnd_shuffle(vect.begin() + 1, vect.end(), m_gen);
        myVectUtils::Print(vect);
        cout << endl;
        Individual indiv(vect, &m_map);
        // indiv.PrintGenes();

        m_population.push_back(indiv);
    }
};

void Population::CheckBounds(void)
{
    for (int i = 0; i < m_ind; i++)
    {
        vector<int> temp_genome = m_population.at(i).GetGenes();

        // Check if all elements are within the range of 1 to m_ngenes
        for (auto it = temp_genome.begin(); it != temp_genome.end(); ++it)
        {
            if (*it < 1 || *it > m_ngenes)
            {
                cerr << "Error: Genome element " << *it << " is out of range." << endl;
                exit(-1);
            }
        }

        // Check if the first element is 1
        if (temp_genome.front() != 1)
        {
            cerr << "Error: First element of genome is not 1." << endl;
            exit(-1);
        }
    }
    cout << "population generated correctly within bounds" << endl;
}

void Population::SortByFitness(void)
{
    // Calculate fitness for each individual
    for (auto &indiv : m_population)
    {
        indiv.Fitness();
    }

    // Sort population by fitness in ascending order
    std::sort(m_population.begin(), m_population.end(),
              [](const Individual &indiv1, const Individual &indiv2)
              {
                  return indiv1.GetFitness() < indiv2.GetFitness();
              });
    cout << "population sorted by fitness" << endl;
}

void Population::EvolutionaryStep(void)
{
    vector<Individual> new_population;

    for(int i = 0; i < m_ind; i++){
        pair<Individual, Individual> parents = Selection();
        pair<Individual, Individual> crossed_parents = Crossover(parents);
        Individual son1 = Mutation(crossed_parents.first);
        Individual son2 = Mutation(crossed_parents.second);

        new_population.push_back(son1);
        new_population.push_back(son2);
    }

    m_population = new_population;
}

void Population ::PrintFitness()
{
    string option_string;
    int option = m_map.GetOption();
    if(option == 0)  option_string = "circumference";
    if(option == 1)  option_string = "square";
    string filenamebest = "best_fitness_" + option_string + ".dat";
    string filenamebesthalf = "besthalf_fitness_" + option_string + ".dat";

    ofstream best, besthalf;

    best.open(filenamebest, ios::app);
    best << GetIndividual(0).GetFitness() << endl;
    best.close();

    besthalf.open(filenamebesthalf, ios::app);
    besthalf << GetFitnessHalf() << endl;
    besthalf.close();

}

void Population ::PrintBestPath()
{
    string option_string;
    int option = m_map.GetOption();
    if(option == 0)  option_string = "circumference";
    if(option == 1)  option_string = "square";
    string filenamebest = "best_path_" + option_string + ".dat";

    ofstream best;

    vector<int> sequence = GetIndividual(0).GetGenes();

    best.open(filenamebest, ios::app);
    for(int i = 0; i < m_ngenes; i++){
        best << sequence.at(i) - 1 << " " <<  sequence.at((i+1)%m_ngenes ) - 1 << endl;
    }
    best.close();
    cout << "Final Optimized Path printed " << endl;
}

double Population ::GetFitnessHalf()
{
    double sum = 0;
    double sum2 = 0;
    for (int i = 0; i < (m_ind / 2); i++){
        sum += GetIndividual(i).GetFitness();
        sum2 += pow(GetIndividual(i).GetFitness(), 2);
    }
    double ave = sum / double(m_ind / 2);
    double ave2 = sum2 / double(m_ind / 2);
    double variance = ave2 - ave*ave;
    cout << variance << endl;
    return sum / double(m_ind / 2);
}

void Population ::FullEvolution(int n_steps)
{
    for (int i = 0; i < n_steps; i++)
    {
        cout << "Step " << i + 1 << endl;
        EvolutionaryStep();
        CheckBounds();
        SortByFitness();
        PrintFitness();
    }
    m_gen.SaveSeed();
}

pair<Individual, Individual> Population::Selection(void)
{
    double r = m_gen.Rannyu();
    double s = m_gen.Rannyu();
    double p_exp = 2;
    int j = int(m_ind * pow(r, p_exp));
    int k = int(m_ind * pow(s, p_exp));
    pair<Individual, Individual> selected = make_pair(m_population[j], m_population[k]);

    return selected;
}

Individual Population::Mutation(Individual indiv)
{
    double r = m_gen.Rannyu();
    vector<int> temp_genes = indiv.GetGenes();

    if (r < 0.1) // mutation 1, swapping a pair of genes
    {
        swap(temp_genes.at(int(m_gen.Rannyu(1, m_ngenes))), temp_genes.at(int(m_gen.Rannyu(1, m_ngenes))));
    }
    if (r < 0.3 && r > 0.2) // mutation 2, shifting a block of m+1 genes
    {
        int start = int(m_gen.Rannyu(1, m_ngenes - 1));
        int m = int(m_gen.Rannyu(0, m_ngenes - start - 2));
        int n = int(m_gen.Rannyu(1, m_ngenes - m - start - 1));
        rotate(temp_genes.begin() + start, temp_genes.begin() + start + m + 1, temp_genes.begin() + start + m + n + 1);
    }
    if (r < 0.5 && r > 0.4) // mutation 3, exchanging position of two blocks of size m
    {
        int m = int(m_gen.Rannyu(0, temp_genes.size() / 2 - 1));  // index tra 1 e N/2 (escluso)
        int idx1 = int(m_gen.Rannyu(1, temp_genes.size() / 2 - m));        // index tra 1 e N/2 - (m+1)   dove m+1 sono le celle del blocco
        int idx2 = int(m_gen.Rannyu(idx1 + m + 1, temp_genes.size() - m)); // index tra (m+1) + idx1 e N - (m+1) incluso
        swap_ranges(temp_genes.begin() + idx1, temp_genes.begin() + idx1 + m + 1, temp_genes.begin() + idx2);
    }
    if (r < 0.7 && r > 0.6) // mutation 4, reversing the order of a subvector within the whole genome
    {
        int start = int(m_gen.Rannyu(1, m_ngenes - 1));
        int end = int(m_gen.Rannyu(start, m_ngenes));
        reverse(temp_genes.begin() + start, temp_genes.begin() + end + 1);
    }

    Individual new_indiv(temp_genes, &m_map);
    return new_indiv;
}

pair<Individual, Individual> Population::Crossover(pair<Individual, Individual> parents)
{
    double r = m_gen.Rannyu();
    if (r > 0.5)
    {

        int m = int(m_gen.Rannyu(1, m_ngenes));
        vector<int> par1 = parents.first.GetGenes();
        vector<int> par2 = parents.second.GetGenes();
        vector<int> son1;
        vector<int> son2;

        son1.insert(son1.begin(), par1.begin(), par1.begin() + m);
        son2.insert(son2.begin(), par2.begin(), par2.begin() + m);

        // add tail to son 1
        vector<int>::iterator it1;

        for (std::vector<int>::iterator it = par2.begin(); it != par2.end(); ++it)
        {
            int new_value = *it;
            it1 = find(son1.begin(), son1.end(), new_value);
            if (it1 == son1.end())
                son1.push_back(new_value);
        }

        // add tail to son 2
        vector<int>::iterator it2;

        for (std::vector<int>::iterator it = par1.begin(); it != par1.end(); ++it)
        {
            int new_value = *it;
            it1 = find(son2.begin(), son2.end(), new_value);
            if (it1 == son2.end())
                son2.push_back(new_value);
        }

        Individual firstson(son1, &m_map);
        Individual secondson(son2, &m_map);

        parents = make_pair(firstson, secondson);
    }
    return parents;
}
