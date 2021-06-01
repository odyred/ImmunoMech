using MGroup.Analyzers;
using MGroup.FEM.Entities;
using MGroup.Problems;
using MGroup.Solvers.Direct;
using MGroup.Stochastic.Interfaces;
using MGroup.Stochastic.Structural.StochasticRealizers;

namespace MGroup.Stochastic.Structural
{
    public class StructuralSystemEvaluator : ISystemRealizer, ISystemResponseEvaluator
    {
        public double YoungModulus { get; }
        public IStochasticDomainMapper DomainMapper;
        public KarhunenLoeveCoefficientsProvider StochasticRealization { get; }
        public GiannisModelBuilder ModelBuilder { get; }
        private Model currentModel;
        int karLoeveTerms = 4;
        double[] domainBounds = new double[2] { 0, 1 };
        double sigmaSquare = .01;
        double meanValue = 1;
        int partition = 21;
        double correlationLength = 1.0;
        bool isGaussian = true;
        int PCorder = 1;
        bool midpointMethod = true;

        /// <summary>Initializes a new instance of the <see cref="StructuralSystemEvaluator"/> class.</summary>
        /// <param name="youngModulus">The young modulus.</param>
        /// <param name="domainMapper">The domain mapper.</param>
        public StructuralSystemEvaluator(double youngModulus, IStochasticDomainMapper domainMapper)
        {
            YoungModulus = youngModulus;
            DomainMapper = domainMapper;
            ModelBuilder = new GiannisModelBuilder();
            StochasticRealization = new KarhunenLoeveCoefficientsProvider(partition, youngModulus, midpointMethod,
                isGaussian, karLoeveTerms, domainBounds, sigmaSquare, correlationLength);
        }

        /// <summary>Realizes the specified iteration.</summary>
        /// <param name="iteration">The iteration.</param>
        public void Realize(int iteration)
        {
            currentModel = ModelBuilder.GetModel(StochasticRealization, DomainMapper, iteration);
        }



        /// <summary>Evaluates the specified iteration.</summary>
        /// <param name="iteration">The iteration.</param>
        /// <returns></returns>
        public double[] Evaluate(int iteration)
        {
            var solver = new SkylineSolver.Builder().BuildSolver(currentModel);
            var provider = new ProblemStructural(currentModel, solver);
            var childAnalyzer = new LinearAnalyzer(currentModel, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(currentModel, solver, provider, childAnalyzer);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            return new[] { solver.LinearSystems[0].Solution[56], solver.LinearSystems[0].Solution[58] };
        }

    }
}
