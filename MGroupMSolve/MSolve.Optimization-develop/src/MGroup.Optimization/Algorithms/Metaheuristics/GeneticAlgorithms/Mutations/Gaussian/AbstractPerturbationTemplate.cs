using MGroup.Optimization.Commons;

namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Mutations.Gaussian
{
    public abstract class AbstractPerturbationTemplate : IMutationStrategy<double>
    {
        public void Apply(Individual<double>[] population)
        {
            foreach (var individual in population)
            {
                double[] chromosome = individual.Chromosome;
                individual.Chromosome = VectorOperations.Add(individual.Chromosome, ComputePerturbations());
            }
        }

        protected abstract double[] ComputePerturbations();
    }
}
