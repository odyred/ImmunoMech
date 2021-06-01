﻿using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using ISAAR.MSolve.FEM; 
using ISAAR.MSolve.FEM.Elements; 
using ISAAR.MSolve.FEM.Entities; 
using ISAAR.MSolve.FEM.Materials; 
using ISAAR.MSolve.Materials; 
using ISAAR.MSolve.SamplesConsole; 
using ISAAR.MSolve.Solvers.Interfaces; 

namespace ISAAR.MSolve.SamplesConsole
{
    public class ProgramElegxoiDdm
    {
        public static void SolveRVEExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(1, new Subdomain() { ID = 1 });

            // EPILOGH MONTELOU
            int model__builder_choice;
            model__builder_choice =4;   // 9 einai to megalo me to renumbering pou tsekaretai

            
            if (model__builder_choice == 1) // 
            { DddmExamplesBuilder.Reference1RVEExample10000(model); }
            if (model__builder_choice == 2) // 
            { DddmExamplesBuilder.Reference1RVEExample10000_Hexaonly(model); }
            if (model__builder_choice == 8) // 
            { RVEExamplesBuilder.Reference2RVEExample10000withRenumbering(model); }
            //if (model__builder_choice == 4) // 
            //{ DddmExamplesBuilder.Reference1RVEExample10000_Hexaonly(model); }
            if (model__builder_choice == 4) // 
            { RVEExamplesBuilder_Debug_ody.FewElementsRVECheckExample2GrapheneSheets(model); }
                //RVEExamplesBuilder.FewElementsRVECheckExample2GrapheneSheets(model); }
            if (model__builder_choice == 5) // 
            { RVEExamplesBuilder.Reference2RVEExample100_000withRenumbering_mono_hexa(model); }
            if (model__builder_choice == 6) // 
            { RVEExamplesBuilder.Reference2RVEExample500_000withRenumbering_mono_hexa(model); }

            // i)
            //DddmExamplesBuilder.MakeModelDictionariesZeroBasedForDecomposer(model);


            model.ConnectDataStructures();

            // ii)
            //AutomaticDomainDecomposer domainDecomposer = new AutomaticDomainDecomposer(model, 2); //2o orisma arithmoos subdomains
            //domainDecomposer.UpdateModel();


            //comment section 1 palaia version
            //SolverSkyline solver = new SolverSkyline(model);
            //ProblemStructural provider = new ProblemStructural(model, solver.SubdomainsDictionary);

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            linearSystems[1] = new SkylineLinearSystem(1, model.Subdomains[0].Forces);
            
            ProblemStructural provider = new ProblemStructural(model, linearSystems);


            // PARADEIGMA A: LinearAnalyzer analyzer = new LinearAnalyzer(solver, solver.SubdomainsDictionary);
            //SolverSkyline2 solver = new SolverSkyline2(linearSystems[1]); //H MARIA XRHSIMOPOIEI TON sklinesolver 
            //LinearAnalyzer childAnalyzer = new LinearAnalyzer(solver, linearSystems);
            //---------------------------------------------------------------------------------------------------------------------------------

            // PARADEIGMA B: Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, 17, model.TotalDOFs);//1. increments einai to 17 (arxika eixame thesei2 26 incr)
            //PALIA DIATUPWSH: NewtonRaphsonNonLinearAnalyzer analyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystems, provider, 10, 48); 
            // NEA DIATUPWSH:
            var solver = new SolverSkyline(linearSystems[1]);
            var linearSystemsArray = new[] { linearSystems[1] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            var increments = 2;
            var childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);
            childAnalyzer.SetMaxIterations = 100;
            childAnalyzer.SetIterationsForMatrixRebuild = 1;
            //---------------------------------------------------------------------------------------------------------------------------------


            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);

            childAnalyzer.LogFactories[1] = new LinearAnalyzerLogFactory(new int[] { 47 });





            //comment section 2 palaia version
            //int increments = 1;
            //Analyzers.NewtonRaphsonNonLinearAnalyzer3 analyzer = new NewtonRaphsonNonLinearAnalyzer3(solver, solver.SubdomainsDictionary, provider, increments, model.TotalDOFs);//1. increments einai to 1 (arxika eixame thesei2 26 incr)
            //StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, analyzer, solver.SubdomainsDictionary);
            //analyzer.SetMaxIterations = 100;
            //analyzer.SetIterationsForMatrixRebuild = 1;



            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


        }

        //static void Main(string[] args)
        //{
        //    SolveRVEExample(); //|
        //}

    }
}
