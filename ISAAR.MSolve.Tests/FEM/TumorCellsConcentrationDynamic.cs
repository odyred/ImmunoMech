using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Readers;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public class TumorCellsConcentrationDynamic
    {
        private const int subdomainID = 0;

        [Fact]
        private static void RunTest()
        {
            Model model = CreateModel();
            IVectorView solution = SolveModel(model);
            Assert.True(CompareResults(solution));
        }

        private static bool CompareResults(IVectorView solution)
        {
            var comparer = new ValueComparer(1E-5);

            //                                                   dofs:   1,   2,   4,   5,   7,   8
            var expectedSolution = Vector.CreateFromArray(new double[] { 150, 200, 150, 200, 150, 200 });
            int numFreeDofs = 6;
            if (solution.Length != 6) return false;
            for (int i = 0; i < numFreeDofs; ++i)
            {
                if (!comparer.AreEqual(expectedSolution[i], solution[i])) return false;
            }
            return true;
        }

        private static Model CreateModel()
        {
            Model model = new Model();
            string filename = "..\\..\\..\\InputFiles\\Cantilever2D.txt";
            ComsolMeshReader modelReader = new ComsolMeshReader(model, filename);
            modelReader.CreateModelFromFile();
            return model;
        }

        private static IVectorView SolveModel(Model model)
        {
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var provider = new ProblemThermal(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new ThermalDynamicAnalyzer(model, solver, provider, childAnalyzer, 0.5, 0.5, 1000);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[subdomainID].Solution;
        }
    }
}

