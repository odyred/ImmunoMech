﻿using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Problems.Structural.Elements;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.SamplesConsole
{
    public class TrussExample
    {
        private const int subdomainID = 0;

        public static IList<Node> CreateNodes()
        {
            IList<Node> nodes = new List<Node>();
            Node node1 = new Node( id: 1, x: 0, y:  0 );
            Node node2 = new Node( id: 2, x: 0, y:  40 );
            Node node3 = new Node( id: 3, x: 40, y:  40 );

            nodes.Add(node1);
            nodes.Add(node2);
            nodes.Add(node3);

            return nodes;
        }

        public static void Run()
        {
            double youngMod = 10e6;
            //double poisson = 0.3;
            double loadX = 500;
            double loadY = 300;
            double sectionArea = 1.5;

            IList<Node> nodes = TrussExample.CreateNodes();

            var model = new Model();

            model.SubdomainsDictionary.Add(0, new Subdomain(subdomainID));

            for (int i = 0; i < nodes.Count; i++)
            {
                model.NodesDictionary.Add(i + 1, nodes[i]);
            }

            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            model.NodesDictionary[1].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            model.NodesDictionary[2].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            model.NodesDictionary[2].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });


            var element1 = new Element()
            {
                ID = 1,
                ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea }
            };

            var element2 = new Element()
            {
                ID = 2,
                ElementType = new Rod2D(youngMod) { Density = 1, SectionArea = sectionArea }
            };

            element1.AddNode(model.NodesDictionary[1]);
            element1.AddNode(model.NodesDictionary[3]);

            element2.AddNode(model.NodesDictionary[2]);
            element2.AddNode(model.NodesDictionary[3]);

            model.ElementsDictionary.Add(element1.ID, element1);
            model.ElementsDictionary.Add(element2.ID, element2);

            model.SubdomainsDictionary[subdomainID].Elements.Add(element1);
            model.SubdomainsDictionary[subdomainID].Elements.Add(element2);

            model.Loads.Add(new Load() { Amount = loadX, Node = model.NodesDictionary[3], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = loadY, Node = model.NodesDictionary[3], DOF = StructuralDof.TranslationY });

            // Solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Output requests
            var logFactory = new TotalDisplacementsLog.Factory(model.SubdomainsDictionary[subdomainID]);
            logFactory.WatchDof(model.NodesDictionary[3], StructuralDof.TranslationX);
            logFactory.WatchDof(model.NodesDictionary[3], StructuralDof.TranslationY);
            childAnalyzer.LogFactories[subdomainID] = logFactory;

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Print output
            var logger = (TotalDisplacementsLog)(childAnalyzer.Logs[subdomainID][0]); //There is a list of logs for each subdomain and we want the first one
            double ux = logger.GetDisplacementAt(model.NodesDictionary[3], StructuralDof.TranslationX);
            double uy = logger.GetDisplacementAt(model.NodesDictionary[3], StructuralDof.TranslationY);
            Console.WriteLine($"Displacements of Node 3: Ux = {ux}, Uy = {uy}");
        }
    }
}
