﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.SamplesConsole
{
    public class CNTExamples
    {
        public static void CNT_4_4_DisplacementControl()
        {
            double youngModulus = 16710.0;
            double poissonRatio = 0.034;
            double nodalDisplacement = 23.73;
            double area = 5.594673861218848e-003;
            double inertiaY = 2.490804749753243e-006;
            double inertiaZ = 2.490804749753243e-006;
            double torsionalInertia = inertiaY / 2.0;
            double effectiveAreaY = area;
            double effectiveAreaZ = area;

            string workingDirectory = @"..\..\..\Resources\Beam3DInputFiles";
            string geometryFileName = "CNT-4-4-L=25-Geometry.inp";
            string connectivityFileName = "CNT-4-4-L=25-ConnMatr.inp";
            int increments = 20;

            int nNodes = File.ReadLines(workingDirectory + '\\' + geometryFileName).Count();
            int nElems = File.ReadLines(workingDirectory + '\\' + connectivityFileName).Count();
            int monitorNode_1 = nNodes - 1;
            int monitorNode_2 = nNodes;

            // Create new 3D material
            var material_1 = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            var material_2 = new ElasticMaterial3D
            {
                YoungModulus = 100.0 * youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Model creation
            Model model = new Model();

            // Subdomains
            model.SubdomainsDictionary.Add(1, new Subdomain(1));

            // Node creation
            IList<Node> nodes = new List<Node>();
            using (TextReader reader = File.OpenText(workingDirectory + '\\' + geometryFileName))
            {
                for (int i = 0; i < nNodes; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int nodeID = int.Parse(bits[0]);
                    double nodeX = double.Parse(bits[1]);
                    double nodeY = double.Parse(bits[2]);
                    double nodeZ = double.Parse(bits[3]);
                    nodes.Add(new Node( id: nodeID, x: nodeX, y:  nodeY, z: nodeZ ));
                    model.NodesDictionary.Add(nodeID, nodes[i]);
                }
            }

            // Constraints
            IList<Node> constraintsNodes = new List<Node>();
            constraintsNodes.Add(nodes[1 - 1]);
            constraintsNodes.Add(nodes[2 - 1]);
            constraintsNodes.Add(nodes[617 - 1]);
            constraintsNodes.Add(nodes[618 - 1]);
            constraintsNodes.Add(nodes[1029 - 1]);
            constraintsNodes.Add(nodes[1030 - 1]);
            constraintsNodes.Add(nodes[1441 - 1]);
            constraintsNodes.Add(nodes[1442 - 1]);
            constraintsNodes.Add(nodes[1649 - 1]);

            for (int i = 0; i < 9; i++)
            {
                int iNode = constraintsNodes[i].ID;
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
            }

            // Applied displacement
            model.NodesDictionary[monitorNode_2].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY, Amount = nodalDisplacement });

            // Create new Beam3D section and element
            var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

            // Generate elements
            using (TextReader reader = File.OpenText(workingDirectory + '\\' + connectivityFileName))
            {
                for (int i = 0; i < (nElems - 16); i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
                    IList<Node> elementNodes = new List<Node>();
                    elementNodes.Add(model.NodesDictionary[node1]);
                    elementNodes.Add(model.NodesDictionary[node2]);
                    // create element
                    var beam_1 = new Beam3DCorotationalQuaternion(elementNodes, material_1, 7.85, beamSection);
                    var element = new Element { ID = elementID, ElementType = beam_1 };
                    // Add nodes to the created element
                    element.AddNode(model.NodesDictionary[node1]);
                    element.AddNode(model.NodesDictionary[node2]);
                    // beam stiffness matrix
                    var a = beam_1.StiffnessMatrix(element);
                    // Add beam element to the element and subdomains dictionary of the model
                    model.ElementsDictionary.Add(element.ID, element);
                    model.SubdomainsDictionary[1].Elements.Add(element);
                }
                for (int i = (nElems - 16); i < nElems; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
                    IList<Node> elementNodes = new List<Node>();
                    elementNodes.Add(model.NodesDictionary[node1]);
                    elementNodes.Add(model.NodesDictionary[node2]);
                    // create element
                    var beam_2 = new Beam3DCorotationalQuaternion(elementNodes, material_2, 7.85, beamSection);
                    var element = new Element { ID = elementID, ElementType = beam_2 };
                    // Add nodes to the created element
                    element.AddNode(model.NodesDictionary[node1]);
                    element.AddNode(model.NodesDictionary[node2]);
                    // beam stiffness matrix
                    var a = beam_2.StiffnessMatrix(element);
                    // Add beam element to the element and subdomains dictionary of the model
                    model.ElementsDictionary.Add(element.ID, element);
                    model.SubdomainsDictionary[1].Elements.Add(element);
                }
            }

            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Choose the provider of the problem -> here a structural problem
            var provider = new ProblemStructural(model, solver);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[1]) };
            var childAnalyzerBuilder = new DisplacementControlAnalyzer.Builder(model, solver, provider, increments);
            var childAnalyzer = childAnalyzerBuilder.Build();

            // Choose parent analyzer -> Parent: Static
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Request output
            //int monDOF = linearSystems[1].Solution.Length - 5;  //(6 * monitorNode_2 - 4);
            int monDOF = 9841;
            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { monDOF });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Check output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[1][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monDOF}, u = {log.DOFValues[monDOF]}");
            //Assert.Equal(-21.2328445476855, log.DOFValues[monDOF], 8);
        }

        public static void CNT_4_4_NewtonRaphson()
        {
            double youngModulus = 16710.0;
            double poissonRatio = 0.034;
            double nodalLoad = 1.2;
            double area = 5.594673861218848e-003;
            double inertiaY = 2.490804749753243e-006;
            double inertiaZ = 2.490804749753243e-006;
            double torsionalInertia = inertiaY / 2.0;
            double effectiveAreaY = area;
            double effectiveAreaZ = area;
            string workingDirectory = @"..\..\..\Resources\Beam3DInputFiles";
            string geometryFileName = "CNT-4-4-L=25-Geometry.inp";
            string connectivityFileName = "CNT-4-4-L=25-ConnMatr.inp";
            int increments = 100;

            //Read number of nodes and number of elements from input files
            int nNodes = File.ReadLines(workingDirectory + '\\' + geometryFileName).Count();
            int nElems = File.ReadLines(workingDirectory + '\\' + connectivityFileName).Count();
            int monitorNode_1 = nNodes - 1;
            int monitorNode_2 = nNodes;

            // Create new 3D material
            var material_1 = new ElasticMaterial3D
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio,
            };

            var material_2 = new ElasticMaterial3D
            {
                YoungModulus = 100.0 * youngModulus,
                PoissonRatio = poissonRatio,
            };

            // Model creation
            Model model = new Model();

            // Subdomains
            model.SubdomainsDictionary.Add(1, new Subdomain(1));

            // Node creation
            IList<Node> nodes = new List<Node>();
            using (TextReader reader = File.OpenText(workingDirectory + '\\' + geometryFileName))
            {
                for (int i = 0; i < nNodes; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int nodeID = int.Parse(bits[0]);
                    double nodeX = double.Parse(bits[1]);
                    double nodeY = double.Parse(bits[2]);
                    double nodeZ = double.Parse(bits[3]);
                    nodes.Add(new Node( id: nodeID, x: nodeX, y:  nodeY, z: nodeZ ));
                    model.NodesDictionary.Add(nodeID, nodes[i]);
                }
            }

            // Constraints
            IList<Node> constraintsNodes = new List<Node>();
            constraintsNodes.Add(nodes[1 - 1]);
            constraintsNodes.Add(nodes[2 - 1]);
            constraintsNodes.Add(nodes[617 - 1]);
            constraintsNodes.Add(nodes[618 - 1]);
            constraintsNodes.Add(nodes[1029 - 1]);
            constraintsNodes.Add(nodes[1030 - 1]);
            constraintsNodes.Add(nodes[1441 - 1]);
            constraintsNodes.Add(nodes[1442 - 1]);
            constraintsNodes.Add(nodes[1649 - 1]);

            for (int i = 0; i < 9; i++)
            {
                int iNode = constraintsNodes[i].ID;
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationZ });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.RotationX });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.RotationY });
                model.NodesDictionary[iNode].Constraints.Add(new Constraint { DOF = StructuralDof.RotationZ });
            }

            // Create new Beam3D section and element
            var beamSection = new BeamSection3D(area, inertiaY, inertiaZ, torsionalInertia, effectiveAreaY, effectiveAreaZ);

            // Generate elements
            using (TextReader reader = File.OpenText(workingDirectory + '\\' + connectivityFileName))
            {
                for (int i = 0; i < (nElems - 16); i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
                    IList<Node> elementNodes = new List<Node>();
                    elementNodes.Add(model.NodesDictionary[node1]);
                    elementNodes.Add(model.NodesDictionary[node2]);
                    // create element
                    var beam_1 = new Beam3DCorotationalQuaternion(elementNodes, material_1, 7.85, beamSection);
                    var element = new Element { ID = elementID, ElementType = beam_1 };
                    // Add nodes to the created element
                    element.AddNode(model.NodesDictionary[node1]);
                    element.AddNode(model.NodesDictionary[node2]);
                    // beam stiffness matrix
                    var a = beam_1.StiffnessMatrix(element);
                    // Add beam element to the element and subdomains dictionary of the model
                    model.ElementsDictionary.Add(element.ID, element);
                    model.SubdomainsDictionary[1].Elements.Add(element);
                }
                for (int i = (nElems - 16); i < nElems; i++)
                {
                    string text = reader.ReadLine();
                    string[] bits = text.Split(',');
                    int elementID = int.Parse(bits[0]);
                    int node1 = int.Parse(bits[1]);
                    int node2 = int.Parse(bits[2]);
                    // element nodes
                    IList<Node> elementNodes = new List<Node>();
                    elementNodes.Add(model.NodesDictionary[node1]);
                    elementNodes.Add(model.NodesDictionary[node2]);
                    // create element
                    var beam_2 = new Beam3DCorotationalQuaternion(elementNodes, material_2, 7.85, beamSection);
                    var element = new Element { ID = elementID, ElementType = beam_2 };
                    // Add nodes to the created element
                    element.AddNode(model.NodesDictionary[node1]);
                    element.AddNode(model.NodesDictionary[node2]);
                    // beam stiffness matrix
                    var a = beam_2.StiffnessMatrix(element);
                    // Add beam element to the element and subdomains dictionary of the model
                    model.ElementsDictionary.Add(element.ID, element);
                    model.SubdomainsDictionary[1].Elements.Add(element);
                }
            }

            // Add nodal load values at the top nodes of the model
            model.Loads.Add(new Load() { Amount = nodalLoad, Node = model.NodesDictionary[monitorNode_2], DOF = StructuralDof.TranslationY });

            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Choose the provider of the problem -> here a structural problem
            var provider = new ProblemStructural(model, solver);

            // Choose child analyzer -> Child: NewtonRaphsonNonLinearAnalyzer
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[1]) };
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(model, solver, provider, increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-3;
            LoadControlAnalyzer childAnalyzer = childAnalyzerBuilder.Build();

            // Choose parent analyzer -> Parent: Static
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Request output
            //int monDOF = linearSystems[1].Solution.Length - 5;  //(6 * monitorNode_2 - 4);
            int monDOF = 9841;
            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { monDOF });

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Check output
            DOFSLog log = (DOFSLog)childAnalyzer.Logs[1][0]; //There is a list of logs for each subdomain and we want the first one
            Console.WriteLine($"dof = {monDOF}, u = {log.DOFValues[monDOF]}");
            //Assert.Equal(23.7373042863345, log.DOFValues[monDOF], 8);
        }
    }
}
