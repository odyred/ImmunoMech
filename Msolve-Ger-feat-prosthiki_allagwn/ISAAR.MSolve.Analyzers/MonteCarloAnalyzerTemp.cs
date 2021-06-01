using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Runtime;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.ModelUpdaters;
using ISAAR.MSolve.Solvers.Interfaces;
using Troschuetz.Random.Distributions.Continuous;
using System.IO;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.Analyzers
{ 
    public class MonteCarloAnalyzerTemp : IAnalyzer
    {
        private int currentSimulation = -1;
        private readonly int blockSize = 5;
        private readonly int simulations;
        private readonly int simulationStartFrom = 0;
        private readonly int randomFileSimulations;
        private readonly int DOFToMonitor;
        private readonly IDictionary<int, ILinearSystem> subdomains;
        //private readonly IDictionary<int, IMatrix2D<double>> matrices;
        private readonly IDictionary<int, IMatrix2D>[] matrices;
        private readonly IDictionary<int, SkylineMatrix2D> factorizedMatrices = new Dictionary<int, SkylineMatrix2D>();
        private readonly Model model;
        private readonly Dictionary<int, IAnalyzerLog[]> logs = new Dictionary<int, IAnalyzerLog[]>();
        private readonly IAnalyzerProvider provider;
        private double[][] randomNumbers;
        //private readonly StiffnessMatrixProductionMode stiffnessMatrixProductionMode = StiffnessMatrixProductionMode.Normal;
        private readonly string stiffnessMatrixPath = String.Empty;
        private readonly string randomsReadFileName = String.Empty;
        //private readonly double[] stochasticDomain;
        private IAnalyzer childAnalyzer;
        private IAnalyzer parentAnalyzer = null;
        private readonly string fileNameForLogging = "monteCarlo";
        private readonly IStochasticCoefficientsProvider coefficientsProvider;
        private readonly IStochasticModelBuilder stochasticModelBuilder;
        private readonly List<int> matrixOrder = new List<int>();
        private readonly List<double> matrixMagnitudes = new List<double>();

        public IDictionary<int, SkylineMatrix2D> FactorizedMatrices { get { return factorizedMatrices; } }


        public MonteCarloAnalyzerTemp(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains,
            IStochasticModelBuilder stochasticModelBuilder, int DOFToMonitor, int simulations)
        {
            this.childAnalyzer = embeddedAnalyzer;
            this.provider = provider;
            this.model = model;
            this.subdomains = subdomains;
            this.DOFToMonitor = DOFToMonitor;
            this.simulations = simulations;
            this.childAnalyzer.ParentAnalyzer = this;
            //this.matrices = new Dictionary<int, IMatrix2D<double>>(subdomains.Count);
            //this.matrices = new Dictionary<int, IMatrix2D>[expansionOrder + 1];
            this.stochasticModelBuilder = stochasticModelBuilder;
            //this.typeOfUncertainty = typeOfUncertainty;
            //this.stochasticDomain = stochasticDomain;
        }

        public MonteCarloAnalyzerTemp(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains,
            IStochasticModelBuilder stochasticModelBuilder, int DOFToMonitor, int simulations, string fileNameForLogging)
            : this(model, provider, embeddedAnalyzer, subdomains, stochasticModelBuilder, DOFToMonitor, simulations)
        {
            this.fileNameForLogging = fileNameForLogging;
        }

        public MonteCarloAnalyzerTemp(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains,
            IStochasticModelBuilder stochasticModelBuilder, int DOFToMonitor, int simulations, string fileNameForLogging, string stiffnessMatrixPath)
            : this(model, provider, embeddedAnalyzer, subdomains, stochasticModelBuilder, DOFToMonitor, simulations, fileNameForLogging)
        {
            this.stiffnessMatrixPath = stiffnessMatrixPath;
        }

        public MonteCarloAnalyzerTemp(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains,
            IStochasticModelBuilder stochasticModelBuilder, int DOFToMonitor, int simulations, int blockSize, string fileNameForLogging, string stiffnessMatrixPath)
            : this(model, provider, embeddedAnalyzer, subdomains, stochasticModelBuilder, DOFToMonitor, simulations, fileNameForLogging, stiffnessMatrixPath)
        {
            this.blockSize = blockSize;
        }

        public MonteCarloAnalyzerTemp(Model model, IAnalyzerProvider provider, IAnalyzer embeddedAnalyzer, IDictionary<int, ILinearSystem> subdomains,
            IStochasticModelBuilder stochasticModelBuilder, int DOFToMonitor, int simulations, int blockSize, string fileNameForLogging, string stiffnessMatrixPath, string randomsReadFileName,
            int simulationStartFrom)
            : this(model, provider, embeddedAnalyzer, subdomains, stochasticModelBuilder, DOFToMonitor, simulations, blockSize, fileNameForLogging, stiffnessMatrixPath)
        {
            this.randomsReadFileName = randomsReadFileName;
            this.simulationStartFrom = simulationStartFrom;
        }

        #region IAnalyzer Members

        public Dictionary<int, IAnalyzerLog[]> Logs { get { return logs; } }
        public IAnalyzer ParentAnalyzer
        {
            get { return parentAnalyzer; }
            set { parentAnalyzer = value; }
        }

        public IAnalyzer ChildAnalyzer
        {
            get { return childAnalyzer; }
            set { childAnalyzer = value; }
        }

        public void BuildMatrices()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
            provider.Reset();
            childAnalyzer.BuildMatrices();
        }


        public void Initialize()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");

        }

        public void Solve()
        {
            if (childAnalyzer == null) throw new InvalidOperationException("Monte Carlo analyzer must contain an embedded analyzer.");
            SolveNormal();
        }

        private void SolveNormal()
        {
            int dofNo = DOFToMonitor;
            string[] values = new string[simulations];
            double[] numberValues = new double[simulations];
            var times = new Dictionary<string, TimeSpan>();
            times.Add("all", TimeSpan.Zero);
            times.Add("element", TimeSpan.Zero);
            times.Add("factorize", TimeSpan.Zero);
            times.Add("solution", TimeSpan.Zero);
            var start = DateTime.Now;
            var iterationCount = new int[simulations - simulationStartFrom];
            int totalIterations = 0;
            for (int i = simulationStartFrom; i < simulationStartFrom + simulations; i++)
            {
                currentSimulation = i;
                stochasticModelBuilder.BuildStochasticModel();
                var e = DateTime.Now;
                BuildMatrices();
                times["element"] += DateTime.Now - e;

                e = DateTime.Now;
                childAnalyzer.Initialize();
                times["factorize"] += DateTime.Now - e;
                GCSettings.LatencyMode = GCLatencyMode.LowLatency;
                e = DateTime.Now;
                childAnalyzer.Solve();

                times["solution"] += DateTime.Now - e;
                GCSettings.LatencyMode = GCLatencyMode.Batch;
                values[i - simulationStartFrom] = subdomains[1].Solution[dofNo].ToString();
                numberValues[i - simulationStartFrom] = subdomains[1].Solution[dofNo];
            }
            MonteCarloMeanValue = numberValues.Average();
            double sumOfSquaresOfDifferences = numberValues.Select(val => (val - MonteCarloMeanValue) * (val - MonteCarloMeanValue)).Sum();
            MonteCarloStandardDeviation = Math.Sqrt(sumOfSquaresOfDifferences / numberValues.Length);
            times["all"] = DateTime.Now - start;

        }

        public double MonteCarloMeanValue { get; set; }
        public double MonteCarloStandardDeviation { get; set; }
        private void MakePreconditioner(int simulation)
        {
            int matrixNo = matrixOrder[simulation + blockSize / 2];
            string name = String.IsNullOrWhiteSpace(stiffnessMatrixPath) ? "K" : stiffnessMatrixPath;
            string path = Path.GetDirectoryName(name);
            string nameOnly = Path.GetFileNameWithoutExtension(name);
            string ext = Path.GetExtension(name);

            foreach (var sub in subdomains)
            {
                var m = new SkylineMatrix2D(new int[0]);
                m.ReadFromFile(String.Format("{0}\\{1}Sub{3}Sim{4}{2}", path, nameOnly, ext, sub.Key, matrixNo));
                m.Factorize(1e-8, new List<IVector>(), new List<int>());
                if (factorizedMatrices.ContainsKey(sub.Key))
                    factorizedMatrices[sub.Key] = m;
                else
                    factorizedMatrices.Add(sub.Key, m);
            }
        }

        private void SolveWithOrder()
        {
            int dofNo = model.Subdomains[0].GlobalNodalDOFsDictionary[6051][DOFType.Y];
            string[] values = new string[simulations];
            //var fileName = String.Format(@"{0}-{1}-{2}.txt", fileNameForLogging, expansionOrder, simulationStartFrom);
            //StreamWriter sw = File.CreateText(fileName);
            //sw.Dispose();
            var times = new Dictionary<string, TimeSpan>();
            times.Add("all", TimeSpan.Zero);
            times.Add("element", TimeSpan.Zero);
            times.Add("factorize", TimeSpan.Zero);
            times.Add("solution", TimeSpan.Zero);
            var start = DateTime.Now;
            var e = start;
            for (int i = simulationStartFrom; i < simulationStartFrom + simulations; i++)
            {
                if (i % blockSize == 0)
                {
                    e = DateTime.Now;
                    MakePreconditioner(i);
                    times["factorize"] += DateTime.Now - e;
                }
                currentSimulation = matrixOrder[i];
                //coefficientsProvider.RandomVariables = randomNumbers[currentSimulation];
                e = DateTime.Now;
                BuildMatrices();
                times["element"] += DateTime.Now - e;
                //if (stiffnessMatrixProductionMode == StiffnessMatrixProductionMode.StoreToDisk) continue;

                e = DateTime.Now;
                childAnalyzer.Initialize();
                childAnalyzer.Solve();
                times["solution"] += DateTime.Now - e;
                values[i] = subdomains[1].Solution[dofNo].ToString();

                //using (sw = File.AppendText(fileName))
                //{
                //    sw.WriteLine(values[i]);
                //}

            }
            times["all"] = DateTime.Now - start;
            var s = new List<string>();
            s.Add(times["all"].ToString());
            s.Add(times["element"].ToString());
            s.Add(times["factorize"].ToString());
            s.Add(times["solution"].ToString());
            File.WriteAllLines(String.Format(@"{0}-Times-{1}-{2}.txt", fileNameForLogging, blockSize, simulationStartFrom), s);
        }


        #endregion
    }
}
