

//TODO: if the design variables correspond to the next iteration, then the initial design variables of the whole optimization
//      cannot be stored. Effectively this means that at least two iterations are executed. Is this a problem?

using MGroup.LinearAlgebra.Reduction;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Optimization.Algorithms.GradientBased.OC.Convergence
{
    public class DesignVariableChangeConvergence : IOptimalityCriteriaConvergence
    {
        private readonly double changeTolerance;
        private IVectorView currentDesignVariables;

        /// <summary>
        /// </summary>
        /// <param name="changeTolerance">
        /// If |x(t+1) - x(t)| &gt; <paramref name="changeTolerance"/> for any design variable x, then this convergence 
        /// criterion is not satisfied.
        /// </param>
        public DesignVariableChangeConvergence(double changeTolerance)
        {
            this.changeTolerance = changeTolerance;
        }

        public bool HasConverged(int currentIteration, double currentObjectiveFunction, IVectorView nextDesignVariables)
        {
            bool result;
            if (currentIteration == 0) result = false;
            else result = nextDesignVariables.Subtract(currentDesignVariables).MaxAbsolute() <= changeTolerance;
            currentDesignVariables = nextDesignVariables.Copy(); //TODO: Perhaps this copy can be avoided
            return result;
        }
    }
}
