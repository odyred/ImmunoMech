using MGroup.Stochastic.Interfaces;
using System;

namespace MGroup.Stochastic
{
    public class MonteCarlo
    {
        public int NoOfIterations { get; }
        public ISystemRealizer SystemRealizer { get; }
        public ISystemResponseEvaluator SystemResponseEvaluator { get; }
        public double[] MonteCarloMeanValue;
        public double[] MonteCarloStandardDeviation;

        /// <summary>Initializes a new instance of the <see cref="MonteCarlo"/> class.</summary>
        /// <param name="noOfIterations">The no of iterations.</param>
        /// <param name="systemRealizer">The system realizer.</param>
        /// <param name="systemResponseEvaluator">The system response evaluator.</param>
        public MonteCarlo(int noOfIterations, ISystemRealizer systemRealizer, ISystemResponseEvaluator systemResponseEvaluator)
        {
            NoOfIterations = noOfIterations;
            SystemRealizer = systemRealizer;
            SystemResponseEvaluator = systemResponseEvaluator;
        }

        /// <summary>Evaluates 1st and 2nd statistical moments of the predesignated response.</summary>
        public void Evaluate()
        {
            SystemRealizer.Realize(0);
            int systemResponseDimension = SystemResponseEvaluator.Evaluate(0).Length;
            var systemResponse = new double[systemResponseDimension];
            MonteCarloMeanValue = new double[systemResponseDimension];
            MonteCarloStandardDeviation = new double[systemResponseDimension];
            systemResponse = SystemResponseEvaluator.Evaluate(0);
            for (int i = 0; i < systemResponse.Length; i++)
            {
                MonteCarloMeanValue[i] += systemResponse[i];
            }

            for (int i = 1; i < NoOfIterations; i++)
            {
                SystemRealizer.Realize(i);
                systemResponse = SystemResponseEvaluator.Evaluate(i);
                for (int j = 0; j < systemResponse.Length; j++)
                {
                    MonteCarloMeanValue[j] += systemResponse[j];
                }
            }

            for (int i = 0; i < systemResponse.Length; i++)
            {
                MonteCarloMeanValue[i] = MonteCarloMeanValue[i] / NoOfIterations;
                MonteCarloStandardDeviation[i] = (systemResponse[i] - MonteCarloMeanValue[i]) * (systemResponse[i] - MonteCarloMeanValue[i]);
                MonteCarloStandardDeviation[i] = Math.Sqrt(MonteCarloStandardDeviation[i] / NoOfIterations);
            }
        }
    }
}
