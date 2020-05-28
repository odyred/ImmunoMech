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
using System.IO;
using ISAAR.MSolve.FEM.Loading.SurfaceLoads;
using static ISAAR.MSolve.FEM.Loading.SurfaceLoads.WeakDirichlet;
using ISAAR.MSolve.FEM.Loading;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System;

namespace ISAAR.MSolve.Tests.FEM
{
    public class TumorCellsConcentrationDynamic2
    {
        private const int subdomainID = 0;
        [Fact]
        private static void RunTest()
        {
            Model model = CreateModel().Item1;
            ComsolMeshReader2 modelReader = CreateModel().Item2;
            IVectorView solution = SolveModel(model,modelReader);
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

        private static Tuple<Model,ComsolMeshReader2> CreateModel()
        {
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "9hexa.mphtxt");
            ComsolMeshReader2 modelReader = new ComsolMeshReader2(filename);
            Model model = modelReader.CreateModelFromFile();
            //Boundary Conditions
            var flux = new FluxLoad(0);
            var dir1 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, 1);
            });
            var dir2 = new DirichletDistribution(list => {
                return Vector.CreateWithValue(list.Count, 0);
            });
            var weakDirichlet1 = new WeakDirichlet(dir1,modelReader.diffusionCoeff);
            var weakDirichlet2 = new WeakDirichlet(dir2, modelReader.diffusionCoeff);

            var dirichletFactory1 = new SurfaceLoadElementFactory(weakDirichlet1);
            var dirichletFactory2 = new SurfaceLoadElementFactory(weakDirichlet2);
            var fluxFactory = new SurfaceLoadElementFactory(flux);

            //model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 1 });
            //model.NodesDictionary[1].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 1 });
            //model.NodesDictionary[2].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 1 });
            //model.NodesDictionary[3].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 1 });
            //model.NodesDictionary[40].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0 });
            //model.NodesDictionary[41].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0 });
            //model.NodesDictionary[42].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0 });
            //model.NodesDictionary[43].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0 });

            int[] boundaryIDs = new int[] { 0, };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Element element in modelReader.elementBoundaries[boundaryID])
                {
                    IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    var dirichletElement1 = dirichletFactory1.CreateElement(CellType.Quad4, nodes);
                    model.SurfaceLoads.Add(dirichletElement1);
                    //var surfaceElement = new SurfaceLoadElement();
                    //element.ID = TriID;
                    //surfaceElement.ElementType = DirichletElement1;
                    //model.SubdomainsDictionary[0].Elements.Add(dirichletElement1);
                    //model.ElementsDictionary.Add(TriID, surfaceElement);

                    //model.NodesDictionary[surfaceElement.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100 });
                }
            }
            boundaryIDs = new int[] { 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Element element in modelReader.elementBoundaries[boundaryID])
                {
                    IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
                    var dirichletElement2 = dirichletFactory2.CreateElement(CellType.Quad4, nodes);
                    model.SurfaceLoads.Add(dirichletElement2);
                }
            }

            //boundaryIDs = new int[] { 1, 2, 3, 4 };
            //foreach (int boundaryID in boundaryIDs)
            //{
            //    foreach (Element element in modelReader.elementBoundaries[boundaryID])
            //    {
            //        IReadOnlyList<Node> nodes = (IReadOnlyList<Node>)element.Nodes;
            //        var fluxElement = fluxFactory.CreateElement(CellType.Quad4, nodes);
            //        model.SurfaceLoads.Add(fluxElement);
            //    }
            //}

            return new Tuple<Model,ComsolMeshReader2>(model,modelReader);
        }

        private static IVectorView SolveModel(Model model, ComsolMeshReader2 modelReader)
        {
            double[] temp0 = new double[model.Nodes.Count];
            int[] boundaryIDs = new int[] { 0 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Element element in modelReader.elementBoundaries[boundaryID])
                {
                    foreach (Node node in element.Nodes)
                    {
                        temp0[node.ID] = 1;
                    }
                }
            }
            boundaryIDs = new int[] { 5 };
            foreach (int boundaryID in boundaryIDs)
            {
                foreach (Element element in modelReader.elementBoundaries[boundaryID])
                {
                    foreach (Node node in element.Nodes)
                    {
                        temp0[node.ID] = 0;
                    }
                }
            }
            Vector initialTemp = Vector.CreateFromArray(temp0);
            var builder = new DenseMatrixSolver.Builder();
            //builder.IsMatrixPositiveDefinite = false;
            var solver = builder.BuildSolver(model);
            var provider = new ProblemConvectionDiffusion(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new ConvectionDiffusionDynamicAnalyzer(model, solver, provider, childAnalyzer, 2.25e-5, 5, initialTemp);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return parentAnalyzer.temperature[subdomainID];
 //           return solver.LinearSystems[subdomainID].Solution;
        }
    }
}

