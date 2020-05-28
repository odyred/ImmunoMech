﻿using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public class ConvectionDiffusion1DBenchmark
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
            var comparer = new ValueComparer(2E-2);

            //                                                   dofs:   1,   2,   4,   5,   7,   8
            var expectedSolution = Vector.CreateFromArray(new double[] { 0.999035307401193,
                0.997371094382290,
                0.993518945712248,
                0.985437312426501,
                0.970000503922635,
                0.943065362762755,
                0.900046845798943,
                0.837421108911808,
                0.753470801649483,
                0.488434721438016});
            int numFreeDofs = 10;
            if (solution.Length != 10) return false;
            for (int i = 0; i < numFreeDofs; ++i)
            {
                if (!comparer.AreEqual(expectedSolution[i], solution[i])) return false;
            }
            return true;
        }

        private static Model CreateModel()
        {
            var model = new Model();

            // Subdomains
            model.SubdomainsDictionary.Add(0, new Subdomain(subdomainID));

            // Material
            double k = 1.0;
            double U = 0;
            double L = 0;
            double h = 1;
            double crossSectionArea = 1;
            var material = new ConvectionDiffusionMaterial(k, U, L);

            // Nodes
            int numNodes = 10;
            for (int i = 0; i < numNodes; i++)
            {
                model.NodesDictionary.Add(i, new Node(id: i, x: i * h));
            }

            //Elements
            var numElements = numNodes - 1;
            for (int i = 0; i < numElements; i++)
            {
                Node[] startEndNodes = { model.NodesDictionary[i], model.NodesDictionary[i+1] };
                var elementType = new ConvectionDiffusionRod(startEndNodes, crossSectionArea, material);
                var elementWrapper = new Element() { ID = i, ElementType = elementType };
                foreach (var node in startEndNodes) elementWrapper.AddNode(node);
                model.ElementsDictionary[elementWrapper.ID] = elementWrapper;
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }

            // Dirichlet BC
            //double kappa = 1;
            //double startNodeTemperature = 1;
            //double endNodeTemperature = 0;
            //model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = kappa * k / h *startNodeTemperature });
            //model.NodesDictionary[9].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = kappa * k / h * endNodeTemperature });

            // Dirichlet BC
            double kappa = 1;
            double startNodeTemperature = 1;
            double endNodeTemperature = 0;
            model.Loads.Add(new Load() { Amount = kappa * k / h * startNodeTemperature, Node = model.NodesDictionary[0], DOF = ThermalDof.Temperature });
            model.Loads.Add(new Load() { Amount = kappa * k / h * endNodeTemperature, Node = model.NodesDictionary[9], DOF = ThermalDof.Temperature });
            //model.Loads.Add(new Load() { Amount = q / 2.0, Node = model.NodesDictionary[8], DOF = ThermalDof.Temperature });


            ////////model.TimeDependentNodalLoads.Add(new SteadyNodalLoad(q / 2.0) { Node = model.NodesDictionary[2], DOF = DOFType.Temperature });

            return model;
        }

        private static IVectorView SolveModel(Model model)
        {
            double[] temp0 = new double[] { 1,0,0,0,0,0,0,0,0,0 };
            Vector initialTemp = Vector.CreateFromArray(temp0);
            var solver = (new DenseMatrixSolver.Builder()).BuildSolver(model);
            //Gmres solver = (new DenseMatrixSolver.Builder()).BuildSolver(model);
            var provider = new ProblemConvectionDiffusion(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new ConvectionDiffusionDynamicAnalyzer_Beta(model, solver, provider, childAnalyzer, 0.05, 5, initialTemp);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            //return solver.LinearSystems[subdomainID].Solution;
            return parentAnalyzer.temperature[subdomainID];
        }
    }
}