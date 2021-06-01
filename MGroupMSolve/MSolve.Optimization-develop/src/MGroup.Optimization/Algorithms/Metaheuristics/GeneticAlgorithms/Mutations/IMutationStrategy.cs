namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations
{
    public interface IMutationStrategy<T>
    {
        void Apply(Individual<T>[] population);
    }
}
