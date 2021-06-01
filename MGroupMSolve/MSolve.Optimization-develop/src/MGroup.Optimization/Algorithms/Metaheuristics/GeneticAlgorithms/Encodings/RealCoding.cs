namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings
{
    class RealCoding : IEncoding<double>
    {
        public double[] ComputeGenotype(double[] phenotype)
        {
            return phenotype;
        }

        public double[] ComputePhenotype(double[] genotype)
        {
            return genotype;
        }
    }
}
