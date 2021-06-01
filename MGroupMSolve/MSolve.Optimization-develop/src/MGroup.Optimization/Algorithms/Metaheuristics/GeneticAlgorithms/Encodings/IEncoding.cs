namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings
{
    public interface IEncoding<T>
    {
        T[] ComputeGenotype(double[] phenotype);
        double[] ComputePhenotype(T[] genotype);
        //int[] IntegerPhenotype(T[] genotype);
    }
}
