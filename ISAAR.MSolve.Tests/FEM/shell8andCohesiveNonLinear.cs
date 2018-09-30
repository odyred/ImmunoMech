﻿using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
//using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.Problems;
//using ISAAR.MSolve.SamplesConsole;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.Solvers.Skyline;
using System.Collections.Generic;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    public static class shell8andCohesiveNonLinear
    {
        private const int subdomainID = 1;

        [Fact]
        private static void RunTest()
        {
            IReadOnlyList<Dictionary<int, double>> expectedDisplacements = GetExpectedDisplacements();
            IncrementalDisplacementsLog computedDisplacements = SolveModel();
            Assert.True(AreDisplacementsSame(expectedDisplacements, computedDisplacements));
        }

        private static bool AreDisplacementsSame(IReadOnlyList<Dictionary<int, double>> expectedDisplacements, IncrementalDisplacementsLog computedDisplacements)
        {
            var comparer = new ValueComparer(1E-13);
            for (int iter = 0; iter < expectedDisplacements.Count; ++iter)
            {
                foreach (int dof in expectedDisplacements[iter].Keys)
                {
                    if (!comparer.AreEqual(expectedDisplacements[iter][dof], computedDisplacements.GetTotalDisplacement(iter, subdomainID, dof)))
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        private static IReadOnlyList<Dictionary<int, double>> GetExpectedDisplacements()
        {
            var expectedDisplacements = new Dictionary<int, double>[5]; //TODO: this should be 11 EINAI ARRAY APO DICTIONARIES

            expectedDisplacements[0] = new Dictionary<int, double> {
    { 0,-1.501306714739379200e-05 }, {11,4.963733738129762700e-06 }, {23,-1.780945407868164000e-05 }, {35,-1.499214801866597700e-05 }, {39,-5.822833969672407000e-05 }};
            expectedDisplacements[1] = new Dictionary<int, double> {
    { 0,-1.500991892603258200e-05 }, {11,4.962619842284933700e-06 }, {23,-1.780557361554878800e-05 }, {35,-1.498958552760043300e-05 }, {39,-5.821676140518532700e-05 }};
            expectedDisplacements[2] = new Dictionary<int, double> {
    { 0,-3.001954880280654900e-05 }, {11,9.925100656459454300e-06 }, {23,-3.561116405105287500e-05 }, {35,-2.997946837567278200e-05 }, {39,-1.164336113147122200e-04 }};
            expectedDisplacements[3] = new Dictionary<int, double> {
    { 0,-3.074327250557275400e-05 }, {11,1.064972618934132200e-05 }, {23,-3.846410374899734500e-05 }, {35,-3.069783728673289500e-05 }, {39,-1.191612724602051200e-04 }};
            expectedDisplacements[4] = new Dictionary<int, double> {
    { 0,-3.074281618470940200e-05 }, {11,1.064926767847933800e-05 }, {23,-3.846254167902190800e-05 }, {35,-3.069737876105337200e-05 }, {39,-1.191596225032377100e-04 }};


            return expectedDisplacements;
        }

        private static IncrementalDisplacementsLog SolveModel()
        {
            VectorExtensions.AssignTotalAffinityCount();
            Model model = new Model();
            model.SubdomainsDictionary.Add(subdomainID, new Subdomain() { ID = subdomainID });

            ShellAndCohesiveRAM_11tlkShellPaktwsh(model);


            model.ConnectDataStructures();            

            var linearSystems = new Dictionary<int, ILinearSystem>(); //I think this should be done automatically 
            linearSystems[subdomainID] = new SkylineLinearSystem(subdomainID, model.Subdomains[0].Forces);

            ProblemStructural provider = new ProblemStructural(model, linearSystems);

            var solver = new SolverSkyline(linearSystems[subdomainID]);
            var linearSystemsArray = new[] { linearSystems[subdomainID] };
            var subdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.Subdomains[0]) };
            var subdomainMappers = new[] { new SubdomainGlobalMapping(model.Subdomains[0]) };

            var increments = 2;
            var childAnalyzer = new NewtonRaphsonNonLinearAnalyzer(solver, linearSystemsArray, subdomainUpdaters, subdomainMappers, provider, increments, model.TotalDOFs);

            var watchDofs = new Dictionary<int, int[]>();
            watchDofs.Add(subdomainID, new int[5] { 0, 11, 23, 35, 39 });
            var log1 = new IncrementalDisplacementsLog(watchDofs);
            childAnalyzer.IncrementalDisplacementsLog = log1;


            childAnalyzer.SetMaxIterations = 100;
            childAnalyzer.SetIterationsForMatrixRebuild = 1;

            StaticAnalyzer parentAnalyzer = new StaticAnalyzer(provider, childAnalyzer, linearSystems);
            
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();


            return log1;
        }

        public static void ShellAndCohesiveRAM_11tlkShellPaktwsh(Model model)
        {
            //PROELEFSI: dhmiourgithike kata to ParadeigmataElegxwnBuilder.ShellAndCohesiveRAM_11ShellPaktwsh(model);
            // allaxame to cohesive element
            // gewmetria
            double Tk = 0.5;

            int nodeID = 1;

            double startX = 0;
            double startY = 0;
            double startZ = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ });
                nodeID++;
            }

            // katw strwsh pou tha paktwthei

            startX = 0;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
                nodeID++;
            }

            startX = 0.25;
            for (int l = 0; l < 2; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.5, Z = startZ - 0.5 * Tk });
                nodeID++;
            }

            startX = 0.5;
            for (int l = 0; l < 3; l++)
            {
                model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + l * 0.25, Z = startZ - 0.5 * Tk });
                nodeID++;
            }

            double[][] VH = new double[8][];

            for (int j = 0; j < 8; j++)
            {
                VH[j] = new double[3];
                VH[j][0] = 0;
                VH[j][1] = 0;
                VH[j][2] = 1;
            }
            // perioxh gewmetrias ews edw

            // constraints

            nodeID = 9;
            for (int j = 0; j < 8; j++)
            {
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);
                nodeID++;
            }
            //perioxh constraints ews edw

            // perioxh materials 
            BenzeggaghKenaneCohMat material1 = new BenzeggaghKenaneCohMat()
            {
                T_o_3 = 57, // New load case argurhs NR_shell_coh.m
                D_o_3 = 5.7e-5,
                D_f_3 = 0.0098245610,
                T_o_1 = 57,
                D_o_1 = 5.7e-5,
                D_f_1 = 0.0098245610,
                n_curve = 1.4,
            };

            //ElasticMaterial3D material2 = new ElasticMaterial3D()
            //{
            //    YoungModulus = 1353000,
            //    PoissonRatio = 0.3,
            //};
            ShellElasticMaterial material2 = new ShellElasticMaterial()
            {
                YoungModulus = 1353000,
                PoissonRatio = 0.3,
                ShearCorrectionCoefficientK = 5 / 6,
            };
            // perioxh materials ews edw




            //eisagwgh tou shell element
            double[] Tk_vec = new double[8];
            for (int j = 0; j < 8; j++)
            {
                Tk_vec[j] = Tk;
            }

            Element e1;
            e1 = new Element()
            {
                ID = 1,
                ElementType = new Shell8dispCopyGetRAM_1(material2, 3, 3, 3)
                {
                    oVn_i = VH,
                    tk = Tk_vec,
                }
            };
            e1.NodesDictionary.Add(8, model.NodesDictionary[8]);
            e1.NodesDictionary.Add(3, model.NodesDictionary[3]);
            e1.NodesDictionary.Add(1, model.NodesDictionary[1]);
            e1.NodesDictionary.Add(6, model.NodesDictionary[6]);
            e1.NodesDictionary.Add(5, model.NodesDictionary[5]);
            e1.NodesDictionary.Add(2, model.NodesDictionary[2]);
            e1.NodesDictionary.Add(4, model.NodesDictionary[4]);
            e1.NodesDictionary.Add(7, model.NodesDictionary[7]);

            int subdomainID = 1; // tha mporei kai na dinetai sto hexabuilder opws sto MakeBeamBuilding
            model.ElementsDictionary.Add(e1.ID, e1);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);
            //eisagwgh shell ews edw

            // eisagwgh tou cohesive element
            int[] coh_global_nodes;
            coh_global_nodes = new int[] { 8, 3, 1, 6, 5, 2, 4, 7, 16, 11, 9, 14, 13, 10, 12, 15 };

            Element e2;
            e2 = new Element()
            {
                ID = 2,
                ElementType = new cohesive_shell_to_hexaCopyGetEmbeRAM_11_tlk(material1, 3, 3)
                {
                    oVn_i = VH,
                    tk = Tk_vec,
                    ShellElementSide = 0,
                }
            };

            for (int j = 0; j < 16; j++)
            {
                e2.NodesDictionary.Add(coh_global_nodes[j], model.NodesDictionary[coh_global_nodes[j]]);
            }

            model.ElementsDictionary.Add(e2.ID, e2);
            model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e2.ID, e2);
            // eisagwgh cohesive ews edw

            // perioxh loads
            double value_ext;
            value_ext = 2 * 2.5 * 0.5;

            int[] points_with_negative_load;
            points_with_negative_load = new int[] { 1, 3, 6, 8 };
            int[] points_with_positive_load;
            points_with_positive_load = new int[] { 2, 4, 5, 7 };

            Load load1;
            Load load2;

            // LOADCASE '' orthi ''
            //for (int j = 0; j < 4; j++)
            //{
            //    load1 = new Load()
            //    {
            //        Node = model.NodesDictionary[points_with_negative_load[j]],
            //        DOF = DOFType.Z,
            //        Amount = -0.3333333 * value_ext,
            //    };
            //    model.Loads.Add(load1);

            //    load2 = new Load()
            //    {
            //        Node = model.NodesDictionary[points_with_positive_load[j]],
            //        DOF = DOFType.Z,
            //        Amount = 1.3333333 * value_ext,
            //    };
            //    model.Loads.Add(load2);
            //}

            // LOADCASE '' orthi '' dixws ta duo prwta fortia  (-0.3333) kai (1.3333)
            for (int j = 0; j < 3; j++)
            {
                load1 = new Load()
                {
                    Node = model.NodesDictionary[points_with_negative_load[j + 1]],
                    DOF = DOFType.Z,
                    Amount = -0.3333333 * value_ext,
                };
                model.Loads.Add(load1);

                load2 = new Load()
                {
                    Node = model.NodesDictionary[points_with_positive_load[j + 1]],
                    DOF = DOFType.Z,
                    Amount = 1.3333333 * value_ext,
                };
                model.Loads.Add(load2);
            }


            // perioxh loads ews edw
        }

    }

}
