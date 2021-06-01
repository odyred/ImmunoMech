using MGroup.Optimization.Problems;

namespace MGroup.Optimization.Constraints.Penalties
{
    public interface IPenaltyStatic
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="fitness"></param>
        /// <param name="design">The current design</param>
        /// <returns>The penalized value for the specified design</returns>
        double Evaluate(double fitness, IDesign design);
    }
}
