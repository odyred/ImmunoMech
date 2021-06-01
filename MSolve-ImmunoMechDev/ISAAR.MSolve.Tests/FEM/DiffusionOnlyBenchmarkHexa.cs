﻿using ISAAR.MSolve.Analyzers;
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
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.Tests.FEM
{
    public class DiffusionOnlyBenchmarkHexa
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
            var comparer = new ValueComparer(1E-3);

            //                                               dofs:       4,       5,       6,       7,       8,       9,      13,      14,      15,      16,      17,      18,      22,      23,      24,      25,      26,  27
            var expectedSolution = Vector.CreateFromArray(new double[] { 135.054, 158.824, 135.054, 469.004, 147.059, 159.327, 178.178, 147.299, 139.469, 147.059, 191.717, 147.059, 135.054, 158.824, 135.054, 469.004, 147.059, 159.327 });
            int numFreeDofs = 18;
            if (solution.Length != 18) return false;
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
            double density = 1.0;
            double k = 1.0;
            double c = 1.0;
            double h = 0;

            // Nodes
            int numNodes = 27;
            var nodes = new Node[numNodes];
            nodes[0] = new Node(id: 0, x: 2.0, y: 2.0, z: 2.0);
            nodes[1] = new Node(id: 1, x: 2.0, y: 1.0, z: 2.0);
            nodes[2] = new Node(id: 2, x: 2.0, y: 0.0, z: 2.0);
            nodes[3] = new Node(id: 3, x: 2.0, y: 2.0, z: 1.0);
            nodes[4] = new Node(id: 4, x: 2.0, y: 1.0, z: 1.0);
            nodes[5] = new Node(id: 5, x: 2.0, y: 0.0, z: 1.0);
            nodes[6] = new Node(id: 6, x: 2.0, y: 2.0, z: 0.0);
            nodes[7] = new Node(id: 7, x: 2.0, y: 1.0, z: 0.0);
            nodes[8] = new Node(id: 8, x: 2.0, y: 0.0, z: 0.0);
            nodes[9] = new Node(id: 9, x: 1.0, y: 2.0, z: 2.0);
            nodes[10] = new Node(id: 10, x: 1.0, y: 1.0, z: 2.0);
            nodes[11] = new Node(id: 11, x: 1.0, y: 0.0, z: 2.0);
            nodes[12] = new Node(id: 12, x: 1.0, y: 2.0, z: 1.0);
            nodes[13] = new Node(id: 13, x: 1.0, y: 1.0, z: 1.0);
            nodes[14] = new Node(id: 14, x: 1.0, y: 0.0, z: 1.0);
            nodes[15] = new Node(id: 15, x: 1.0, y: 2.0, z: 0.0);
            nodes[16] = new Node(id: 16, x: 1.0, y: 1.0, z: 0.0);
            nodes[17] = new Node(id: 17, x: 1.0, y: 0.0, z: 0.0);
            nodes[18] = new Node(id: 18, x: 0.0, y: 2.0, z: 2.0);
            nodes[19] = new Node(id: 19, x: 0.0, y: 1.0, z: 2.0);
            nodes[20] = new Node(id: 20, x: 0.0, y: 0.0, z: 2.0);
            nodes[21] = new Node(id: 21, x: 0.0, y: 2.0, z: 1.0);
            nodes[22] = new Node(id: 22, x: 0.0, y: 1.0, z: 1.0);
            nodes[23] = new Node(id: 23, x: 0.0, y: 0.0, z: 1.0);
            nodes[24] = new Node(id: 24, x: 0.0, y: 2.0, z: 0.0);
            nodes[25] = new Node(id: 25, x: 0.0, y: 1.0, z: 0.0);
            nodes[26] = new Node(id: 26, x: 0.0, y: 0.0, z: 0.0);

            for (int i = 0; i < numNodes; ++i) model.NodesDictionary[i] = nodes[i];

            // Elements
            int numElements = 8;
            var elementFactory = new ConvectionDiffusionElement3DFactory(new ConvectionDiffusionMaterial(k, new double[] { 0, 0, 0 }, 0));
            var elements = new ConvectionDiffusionElement3D[8];
            elements[0] = elementFactory.CreateElement(CellType.Hexa8, new Node[] { nodes[13], nodes[4], nodes[3], nodes[12], nodes[10], nodes[1], nodes[0], nodes[9] });
            elements[1] = elementFactory.CreateElement(CellType.Hexa8, new Node[] { nodes[14], nodes[5], nodes[4], nodes[13], nodes[11], nodes[2], nodes[1], nodes[10] });
            elements[2] = elementFactory.CreateElement(CellType.Hexa8, new Node[] { nodes[16], nodes[7], nodes[6], nodes[15], nodes[13], nodes[4], nodes[3], nodes[12] });
            elements[3] = elementFactory.CreateElement(CellType.Hexa8, new Node[] { nodes[17], nodes[8], nodes[7], nodes[16], nodes[14], nodes[5], nodes[4], nodes[13] });
            elements[4] = elementFactory.CreateElement(CellType.Hexa8, new Node[] { nodes[22], nodes[13], nodes[12], nodes[21], nodes[19], nodes[10], nodes[9], nodes[18] });
            elements[5] = elementFactory.CreateElement(CellType.Hexa8, new Node[] { nodes[23], nodes[14], nodes[13], nodes[22], nodes[20], nodes[11], nodes[10], nodes[19] });
            elements[6] = elementFactory.CreateElement(CellType.Hexa8, new Node[] { nodes[25], nodes[16], nodes[15], nodes[24], nodes[22], nodes[13], nodes[12], nodes[21] });
            elements[7] = elementFactory.CreateElement(CellType.Hexa8, new Node[] { nodes[26], nodes[17], nodes[16], nodes[25], nodes[23], nodes[14], nodes[13], nodes[22] });

            for (int i = 0; i < numElements; ++i)
            {
                var elementWrapper = new Element() { ID = i, ElementType = elements[i] };
                foreach (var node in elements[i].Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary[i] = elementWrapper;
                model.SubdomainsDictionary[subdomainID].Elements.Add(elementWrapper);
            }

            // Dirichlet BC
            model.NodesDictionary[0].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
            model.NodesDictionary[1].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
            model.NodesDictionary[2].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
            model.NodesDictionary[9].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
            model.NodesDictionary[10].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
            model.NodesDictionary[11].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
            model.NodesDictionary[18].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
            model.NodesDictionary[19].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });
            model.NodesDictionary[20].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 100.0 });

            // Neumann BC
            double q = 100;
            model.Loads.Add(new Load() { Amount = q, Node = model.NodesDictionary[6], DOF = ThermalDof.Temperature });
            model.Loads.Add(new Load() { Amount = q, Node = model.NodesDictionary[24], DOF = ThermalDof.Temperature });

            return model;
        }

        private static IVectorView SolveModel(Model model)
        {
            SkylineSolver solver = (new SkylineSolver.Builder()).BuildSolver(model);
            var provider = new ProblemThermal(model, solver);

            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            return solver.LinearSystems[subdomainID].Solution;
        }
    }
}
