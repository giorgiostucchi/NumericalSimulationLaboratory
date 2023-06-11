#ifndef __classes__
#define __classes__

#include <algorithm> // std::random_shuffle
#include <vector>    // std::vector
#include <numeric>
#include "random.h"
#include <string>
#include <utility>
#include "vectoralgebra.h"
#include <iostream>

using namespace std;

void CleanPreviousData(int); //Clean data from previous runs of the program (for the selected map)
void rnd_shuffle(vector<int>::iterator first, vector<int>::iterator last, Random &rnd); //shuffle a vector

class Map
{

public:
  Map(int ncities, int opt);

  int GetNumberCities()
  {
    return m_ncities;
  }
  pair<double, double> GetNthCoordinates(int n)
  {
    return m_cities_coord.at(n);
  }
  int GetOption(){return m_opt;}
  void PrintMap(void);

private:
  vector<pair<double, double>> m_cities_coord;
  int m_opt;
  int m_ncities;
  Random m_gen;
};

class Individual
{
public:
  ~Individual() {}
  Individual(vector<int> orders, Map *myMap) : m_genes(orders), m_map(myMap), m_fitness(0.0) {}
  vector<int> GetGenes() const { return m_genes; }
  void PrintGenes() const { myVectUtils::Print(m_genes); }
  void Fitness();
  double GetFitness() const { return m_fitness; }
  void SetGenes(vector<int> new_genes) { m_genes = new_genes; }

private:
  vector<int> m_genes;
  Map *m_map;
  double m_fitness;
};

class Population
{

public:
  ~Population() {}
  Population(int, Map &);
  void CheckBounds(void);
  Individual GetIndividual(int n) { return m_population.at(n); };
  void SortByFitness(void);
  pair<Individual, Individual> Selection();
  pair<Individual, Individual> Crossover(pair<Individual, Individual>);
  Individual Mutation(Individual);
  void EvolutionaryStep();
  void FullEvolution(int nsteps);
  void PrintFitness(void);
  double GetFitnessHalf(void);
  void PrintBestPath();

/*   these were just here to check whether the mutations worked correctly
  void TrialMutation2(int L)
  {
    Individual ind = GetIndividual(L);
    myVectUtils::Print(ind.GetGenes());
    vector<int> temp_genes = ind.GetGenes();
    int start = int(m_gen.Rannyu(1, m_ngenes - 1)); //ok
    int m = int(m_gen.Rannyu(0, m_ngenes - start - 2));
    int n = int(m_gen.Rannyu(1, m_ngenes - m - start - 1));
    cout <<  *(temp_genes.begin() + start)<< " " << *(temp_genes.begin() + start + m) << " "<< *(temp_genes.begin() + start + m + n -1) << endl;
    cout << *(temp_genes.begin() + 31) << endl;
    rotate(temp_genes.begin() + start, temp_genes.begin() + start + m + 1, temp_genes.begin() + start +m + n + 1);


    myVectUtils::Print(temp_genes);
    cout << "start ,  celle spostate ,  spostamento  = " << start << "  " << m+1 << "  " << n << endl;
    m_gen.SaveSeed();
  };

  void TrialMutation3(int L)
  {
    Individual ind = GetIndividual(L);
    myVectUtils::Print(ind.GetGenes());
    vector<int> temp_genes = ind.GetGenes();
    int m = 14; //int(m_gen.Rannyu(0, temp_genes.size() / 2 - 1));  // index tra 1 e N/2 (escluso)
    int idx1 = int(m_gen.Rannyu(1, temp_genes.size() / 2 - m)); // index tra 1 e N/2 - (m+1)   dove m+1 sono le celle del blocco
    int idx2 = int(m_gen.Rannyu(idx1 + m + 1, temp_genes.size() - m)); // index tra (m+1) + idx1 e N - (m+1) incluso
    cout << "indexes : " << idx1 << "  " << idx2 << endl;
    swap_ranges(temp_genes.begin() + idx1, temp_genes.begin() + idx1 + m + 1, temp_genes.begin() + idx2); 
    myVectUtils::Print(temp_genes);
    m_gen.SaveSeed();
  };

  void TrialMutation4(int L)
  {
    Individual ind = GetIndividual(L);
    myVectUtils::Print(ind.GetGenes());
    vector<int> temp_genes = ind.GetGenes();
    int start = int(m_gen.Rannyu(1, m_ngenes - 1));  //massimo 30 
    int end = int(m_gen.Rannyu(start, m_ngenes)); // se start e' 30 allora e' 31 elemento
    cout << " begin and end : " << start << "  " << end << endl;
    reverse(temp_genes.begin()+start, temp_genes.begin()+end+1);
    myVectUtils::Print(temp_genes);
    m_gen.SaveSeed();
  };
*/


private:
  int m_ind;
  int m_ngenes;
  Map m_map;
  vector<Individual> m_population;
  Random m_gen;
};


#endif