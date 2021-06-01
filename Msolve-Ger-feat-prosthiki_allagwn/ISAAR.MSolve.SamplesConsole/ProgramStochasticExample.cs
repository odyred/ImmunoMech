using System;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Skyline;
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
using System.Runtime.CompilerServices;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Embedding;


namespace ISAAR.MSolve.SamplesConsole
{
    public class ProgramStochasticExample
    {
        public static void SolveStochasticExample()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Tuple<Analyzers.rveMatrixParameters, Analyzers.grapheneSheetParameters> gpmp = RandomGrapheneModelBuilder.Get2GrapheneSheetsRveExampleParameters();
            int graphene_sheets_number=2;
            double b1 = 10; double b2 = 10; double sigma_f = 0.2;
            IList<IStochasticCoefficientsProvider2D> coefficientsProviders=new List<IStochasticCoefficientsProvider2D> { new SpectralRepresentation2DRandomField(b1, b2, sigma_f, 0.01) } ;
            for (int j = 0; j < graphene_sheets_number-1; j++)
            { coefficientsProviders.Add(new SpectralRepresentation2DRandomField(b1, b2, sigma_f, 0.01)); }
            
            Model model = new Model();

            RandomGrapheneModelBuilder randomGrapheneModelBuilder1 = new RandomGrapheneModelBuilder(model,coefficientsProviders,graphene_sheets_number,gpmp);
            randomGrapheneModelBuilder1.BuildStochasticModel();

            //model.ConnectDataStructures();
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

            var increments = 1;
            var childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);
            childAnalyzer.SetMaxIterations = 100;
            childAnalyzer.SetIterationsForMatrixRebuild = 1;
            //---------------------------------------------------------------------------------------------------------------------------------

            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);


            MonteCarloAnalyzerTemp stohasticAnalyzer = new MonteCarloAnalyzerTemp(model, provider, parentAnalyzer, linearSystems, randomGrapheneModelBuilder1, 50, 2);
            stohasticAnalyzer.Initialize();
            stohasticAnalyzer.BuildMatrices();
            stohasticAnalyzer.Solve();
        }





        //private static Tuple<rveMatrixParameters, grapheneSheetParameters> Get2GrapheneSheetsRveExampleParameters()
        //{
        //    rveMatrixParameters mp;
        //    mp = new rveMatrixParameters()
        //    {
        //        E_disp = 3.5, //Gpa
        //        ni_disp = 0.4, // stather Poisson
        //        L01 = 150, // diastaseis
        //        L02 = 150,
        //        L03 = 40,
        //        hexa1 = 2,// diakritopoihsh
        //        hexa2 = 2,
        //        hexa3 = 2
        //    };

        //    grapheneSheetParameters gp;
        //    gp = new grapheneSheetParameters()
        //    {
        //         //parametroi shell
        //        E_shell = 27196.4146610211, // GPa = 1000Mpa = 1000N / mm2
        //        ni_shell = 0.0607, // stathera poisson
        //        elem1 = 1, // DIORTHOSI 2 graohene sheets
        //        elem2 = 2,
        //        L1 = 0.5 * 105,// nm  // DIORTHOSI 2 graphene sheets
        //        L2 = 110,// nm
        //        L3 = 112.5096153846, // nm
        //        a1_shell = 0, // nm
        //        tk = 0.0125016478913782,  // 0.0125016478913782nm //0.125*40,

        //        //parametroi cohesive epifaneias
        //        T_o_3 = 0.05,// Gpa = 1000Mpa = 1000N / mm2
        //        D_o_3 = 0.5, // nm
        //        D_f_3 = 4, // nm
        //        T_o_1 = 0.05,// Gpa
        //        D_o_1 = 0.5, // nm
        //        D_f_1 = 4, // nm
        //        n_curve = 1.4
        //    };

        //    Tuple<rveMatrixParameters, grapheneSheetParameters> gpmp = new Tuple<rveMatrixParameters, grapheneSheetParameters>(mp, gp);
        //    return gpmp;
        //}



    }
}

