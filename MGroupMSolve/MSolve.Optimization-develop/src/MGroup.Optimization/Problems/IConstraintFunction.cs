using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Optimization.Problems
{
    //TODO: use a design class for this.
    public delegate double EqualityConstraint(Vector x);

    public interface IConstraintFunction
    {
        double Evaluate(double[] x);
    }
}
