using System;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Optimization.Algorithms.Metaheuristics.DifferentialEvolution;
using MGroup.Optimization.Benchmarks.Unconstrained;
using MGroup.Optimization.Convergence;
using MGroup.Optimization.Problems;
using Xunit;

namespace MGroup.Optimization.Tests
{
    public static class TestDE
    {
        [Fact]
        public static void Run()
        {
            int seed = 1;
            var rng = new Random(seed);

            OptimizationProblem optimizationProblem = new Rosenbrock();

            var builder = new DifferentialEvolutionAlgorithm.Builder(optimizationProblem);
            builder.PopulationSize = 100;
            builder.MutationFactor = 0.6;
            builder.CrossoverProbability = 0.9;
            builder.ConvergenceCriterion = new MaxFunctionEvaluations(100000);
            builder.RandomNumberGenerator = rng;
            IOptimizationAlgorithm de = builder.Build();

            IOptimizationAnalyzer analyzer = new OptimizationAnalyzer(de);
            analyzer.Optimize();

            double expectedFitness = 0.0;
            var expectedDesign = Vector.CreateWithValue(optimizationProblem.Dimension, 1.0);
            Assert.Equal(expectedFitness, de.BestFitness, 10);
            Assert.True(Vector.CreateFromArray(de.BestPosition).Equals(expectedDesign, 1E-6));
        }
    }
}
