using MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections.FitnessScaling;
using MGroup.Optimization.Commons;
using Troschuetz.Random;

namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Selections
{
    public class RemainderStochasticSamplingWithReplacement<T> : AbstractRemainderSelection<T>
    {
        public RemainderStochasticSamplingWithReplacement(IFitnessScalingStrategy<T> fitnessScaling) :
            base(fitnessScaling, RandomNumberGenerationUtilities.troschuetzRandom)
        {
        }

        public RemainderStochasticSamplingWithReplacement(IFitnessScalingStrategy<T> fitnessScaling, 
                                                          IGenerator randomNumberGenerator):
            base(fitnessScaling, randomNumberGenerator)
        {
        }

        protected override sealed Individual<T>[] RemainderPhase(Individual<T>[] population, double[] fractionalParts, 
            int extrasParentsCount)
        {
            Roulette roulette = Roulette.CreateFromPositive(fractionalParts, RNG);
            var extraParents = new Individual<T>[extrasParentsCount];
            for (int i = 0; i < extrasParentsCount; ++i)
            {
                extraParents[i] = population[roulette.SpinWheelWithBall()];
            }
            return extraParents;
        }
    }
}
