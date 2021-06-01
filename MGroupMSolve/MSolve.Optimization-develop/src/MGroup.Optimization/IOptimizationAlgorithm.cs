namespace MGroup.Optimization
{
    public interface IOptimizationAlgorithm
    {
        double BestFitness
        {
            get;
        }

        double[] BestPosition
        {
            get;
        }

        int CurrentIteration
        { 
            get; 
        }

        int CurrentFunctionEvaluations
        {
            get;
        }

        void Solve();
    }
}