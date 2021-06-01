namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections
{
    public interface ISelectionStrategy<T>
    {
        Individual<T>[][] Apply(Individual<T>[] population, int parentGroupsCount, int parentsPerGroup);
    }
}
