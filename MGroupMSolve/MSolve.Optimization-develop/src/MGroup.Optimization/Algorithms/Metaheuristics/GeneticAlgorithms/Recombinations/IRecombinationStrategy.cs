using MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections;

namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Recombinations
{
    public interface IRecombinationStrategy<T>
    {
        Individual<T>[] Apply(ISelectionStrategy<T> selection, Individual<T>[] population, int offspringsCount);
    }
}
