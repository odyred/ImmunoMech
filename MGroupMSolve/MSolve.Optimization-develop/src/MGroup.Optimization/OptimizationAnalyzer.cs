using MGroup.Optimization.Problems;

namespace MGroup.Optimization
{
    public class OptimizationAnalyzer : IOptimizationAnalyzer
    {
        public OptimizationProblem optimizationProblem;
        public IOptimizationAlgorithm optimizationAlgorithm;

        public OptimizationAnalyzer(IOptimizationAlgorithm optimizationAlgorithm)
        {
            this.optimizationAlgorithm = optimizationAlgorithm;
        }

        void IOptimizationAnalyzer.Optimize()
        {
            optimizationAlgorithm.Solve();
        }
    }
}
