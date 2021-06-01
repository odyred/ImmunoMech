namespace MGroup.Optimization.Convergence
{
    // Convergence criteria must be checked at the END of each iteration of the optimization process
    public interface IConvergenceCriterion
    {
        bool HasConverged(IOptimizationAlgorithm algorithm);
    }
}
