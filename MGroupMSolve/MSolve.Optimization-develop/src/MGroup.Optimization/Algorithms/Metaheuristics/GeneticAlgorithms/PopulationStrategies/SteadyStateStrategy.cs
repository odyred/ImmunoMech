using System;
using MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations;
using MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations;
using MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;

namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.PopulationStrategies
{
    class SteadyStateStrategy<T> : IPopulationStrategy<T>
    {
        public Individual<T>[] CreateNextGeneration(Individual<T>[] originalPopulation, ISelectionStrategy<T> selection, 
                                                    IRecombinationStrategy<T> recombination, IMutationStrategy<T> mutation)
        {
            throw new NotImplementedException();
        }
    }
}
