using MGroup.Optimization.Algorithms.GradientBased.OC.Bisection;
using MGroup.Optimization.Algorithms.GradientBased.OC.Convergence;
using MGroup.Optimization.Problems;

namespace MGroup.Optimization.Algorithms.GradientBased.OC
{
    public class OptimalityCriteriaBuilder
    {
        public double InitialBisectionLimitLower { get; set; } = 0.0;
        public double InitialBisectionLimitUpper { get; set; } = 1E5;
        public IBisectionConvergence BisectionConvergence { get; set; } = new BisectionConvergenceAbsoluteChange(1E-4);
        public double DampingCoeff { get; set; } = 0.5;
        public double MoveLimit { get; set; } = 0.2;
        public IOptimalityCriteriaConvergence OptimalityCriteriaConvergence { get; set; }
            = new DesignVariableChangeConvergence(1E-2);

        public OptimalityCriteria BuildOptimizer(DifferentiableObjectiveFunction objective, EqualityConstraint constraint,
            double boundLower, double boundUpper)
            => new OptimalityCriteria(objective, constraint, boundLower, boundUpper, OptimalityCriteriaConvergence,
                InitialBisectionLimitLower, InitialBisectionLimitUpper, BisectionConvergence, DampingCoeff, MoveLimit);
    }
}
