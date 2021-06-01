using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Optimization.Algorithms.GradientBased.OC.Convergence
{
    public class MaxIterationsConvergence : IOptimalityCriteriaConvergence
    {
        private readonly int maxIterations;

        /// <summary>
        /// </summary>
        /// <param name="maxIterations">
        /// This convergence criterion is not satisfied if 0 &lt;= currentIteration &lt; <paramref name="maxIterations"/>.
        /// </param>
        public MaxIterationsConvergence(int maxIterations)
        {
            this.maxIterations = maxIterations;
        }

        public bool HasConverged(int currentIteration, double currentObjectiveFunction, IVectorView nextDesignVariables)
            => currentIteration >= maxIterations; // Other criteria may cause the > case.
    }
}
