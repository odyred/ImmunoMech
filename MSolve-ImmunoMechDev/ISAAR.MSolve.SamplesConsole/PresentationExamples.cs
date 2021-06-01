﻿using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Mesh.Generation;
using ISAAR.MSolve.Discretization.Mesh.Generation.GMSH;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Logging.VTK;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.Direct;
using ISSAR.MSolve.Discretization.Loads;

namespace ISAAR.MSolve.SamplesConsole
{
    public class PresentationTests
    {
        private static void SolveStaticQuadRetainingWall()
        {
            #region Model
            double youngModulus = 2.1e09;
            double poissonRatio = 0.3;

            var material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
            Model model = new Model();

            model.SubdomainsDictionary.Add(0, new Subdomain(0));

            #region Nodes
            model.NodesDictionary.Add(0, new Node( id: 0, x: 0.0, y:  0.0, z: 0.0 ));
            model.NodesDictionary.Add(1, new Node( id: 1, x: 0.6, y:  0.0, z: 0.0 ));

            model.NodesDictionary.Add(2, new Node( id: 2, x: 0.0, y:  0.5, z: 0.0 ));
            model.NodesDictionary.Add(3, new Node( id: 3, x: 0.6, y:  0.5, z: 0.0 ));
            model.NodesDictionary.Add(4, new Node( id: 4, x: 1.2, y:  0.5, z: 0.0 ));
            model.NodesDictionary.Add(5, new Node( id: 5, x: 1.8, y:  0.5, z: 0.0 ));
            model.NodesDictionary.Add(6, new Node( id: 6, x: 2.4, y:  0.5, z: 0.0 ));
            model.NodesDictionary.Add(7, new Node( id: 7, x: 3.0, y:  0.5, z: 0.0 ));
            model.NodesDictionary.Add(8, new Node( id: 8, x: 3.6, y:  0.5, z: 0.0 ));
            model.NodesDictionary.Add(9, new Node( id: 9, x: 4.3, y:  0.5, z: 0.0 ));

            model.NodesDictionary.Add(10, new Node( id: 10, x: 0.0, y:  1.2, z: 0.0 ));
            model.NodesDictionary.Add(11, new Node( id: 11, x: 0.6, y:  1.2, z: 0.0 ));
            model.NodesDictionary.Add(12, new Node( id: 12, x: 1.2, y:  1.2, z: 0.0 ));
            model.NodesDictionary.Add(13, new Node( id: 13, x: 1.8, y:  1.2, z: 0.0 ));
            model.NodesDictionary.Add(14, new Node( id: 14, x: 2.4, y:  1.2, z: 0.0 ));
            model.NodesDictionary.Add(15, new Node( id: 15, x: 3.0, y:  1.2, z: 0.0 ));
            model.NodesDictionary.Add(16, new Node( id: 16, x: 3.6, y:  1.2, z: 0.0 ));
            model.NodesDictionary.Add(17, new Node( id: 17, x: 4.3, y:  1.2, z: 0.0 ));

            model.NodesDictionary.Add(18, new Node( id: 18, x: 3, y:  1.9, z: 0.0 ));
            model.NodesDictionary.Add(19, new Node( id: 19, x: 3.5756, y:  1.9, z: 0.0 ));

            model.NodesDictionary.Add(20, new Node( id: 20, x: 3, y:  2.6, z: 0.0 ));
            model.NodesDictionary.Add(21, new Node( id: 21, x: 3.5512, y:  2.6, z: 0.0 ));

            model.NodesDictionary.Add(22, new Node( id: 22, x: 3, y:  3.3, z: 0.0 ));
            model.NodesDictionary.Add(23, new Node( id: 23, x: 3.5267, y:  3.3, z: 0.0 ));

            model.NodesDictionary.Add(24, new Node( id: 24, x: 3, y:  4.0, z: 0.0 ));
            model.NodesDictionary.Add(25, new Node( id: 25, x: 3.5023, y:  4.0, z: 0.0 ));

            model.NodesDictionary.Add(26, new Node( id: 26, x: 3, y:  4.7, z: 0.0 ));
            model.NodesDictionary.Add(27, new Node( id: 27, x: 3.4779, y:  4.7, z: 0.0 ));

            model.NodesDictionary.Add(28, new Node( id: 28, x: 3, y:  5.5, z: 0.0 ));
            model.NodesDictionary.Add(29, new Node( id: 29, x: 3.45, y:  5.5, z: 0.0 ));
            #endregion

            #region Elements

            #region element0
            var element0 = new Element()
            {
                ID = 0,
                ElementType = new Quad4(material)
            };
            element0.AddNode(model.NodesDictionary[0]);
            element0.AddNode(model.NodesDictionary[1]);
            element0.AddNode(model.NodesDictionary[3]);
            element0.AddNode(model.NodesDictionary[2]);
            model.ElementsDictionary.Add(element0.ID, element0);
            model.SubdomainsDictionary[0].Elements.Add(element0);
            #endregion

            #region element1
            var element1 = new Element()
            {
                ID = 1,
                ElementType = new Quad4(material)
            };
            element1.AddNode(model.NodesDictionary[2]);
            element1.AddNode(model.NodesDictionary[3]);
            element1.AddNode(model.NodesDictionary[11]);
            element1.AddNode(model.NodesDictionary[10]);
            model.ElementsDictionary.Add(element1.ID, element1);
            model.SubdomainsDictionary[0].Elements.Add(element1);
            #endregion

            #region element2
            var element2 = new Element()
            {
                ID = 2,
                ElementType = new Quad4(material)
            };
            element2.AddNode(model.NodesDictionary[3]);
            element2.AddNode(model.NodesDictionary[4]);
            element2.AddNode(model.NodesDictionary[12]);
            element2.AddNode(model.NodesDictionary[11]);
            model.ElementsDictionary.Add(element2.ID, element2);
            model.SubdomainsDictionary[0].Elements.Add(element2);
            #endregion

            #region element3
            var element3 = new Element()
            {
                ID = 3,
                ElementType = new Quad4(material)
            };
            element3.AddNode(model.NodesDictionary[4]);
            element3.AddNode(model.NodesDictionary[5]);
            element3.AddNode(model.NodesDictionary[13]);
            element3.AddNode(model.NodesDictionary[12]);
            model.ElementsDictionary.Add(element3.ID, element3);
            model.SubdomainsDictionary[0].Elements.Add(element3);
            #endregion

            #region element4
            var element4 = new Element()
            {
                ID = 4,
                ElementType = new Quad4(material)
            };
            element4.AddNode(model.NodesDictionary[5]);
            element4.AddNode(model.NodesDictionary[6]);
            element4.AddNode(model.NodesDictionary[14]);
            element4.AddNode(model.NodesDictionary[13]);
            model.ElementsDictionary.Add(element4.ID, element4);
            model.SubdomainsDictionary[0].Elements.Add(element4);
            #endregion

            #region element5
            var element5 = new Element()
            {
                ID = 5,
                ElementType = new Quad4(material)
            };
            element5.AddNode(model.NodesDictionary[6]);
            element5.AddNode(model.NodesDictionary[7]);
            element5.AddNode(model.NodesDictionary[15]);
            element5.AddNode(model.NodesDictionary[14]);
            model.ElementsDictionary.Add(element5.ID, element5);
            model.SubdomainsDictionary[0].Elements.Add(element5);
            #endregion

            #region element6
            var element6 = new Element()
            {
                ID = 6,
                ElementType = new Quad4(material)
            };
            element6.AddNode(model.NodesDictionary[7]);
            element6.AddNode(model.NodesDictionary[8]);
            element6.AddNode(model.NodesDictionary[16]);
            element6.AddNode(model.NodesDictionary[15]);
            model.ElementsDictionary.Add(element6.ID, element6);
            model.SubdomainsDictionary[0].Elements.Add(element6);
            #endregion

            #region element7
            var element7 = new Element()
            {
                ID = 7,
                ElementType = new Quad4(material)
            };
            element7.AddNode(model.NodesDictionary[8]);
            element7.AddNode(model.NodesDictionary[9]);
            element7.AddNode(model.NodesDictionary[17]);
            element7.AddNode(model.NodesDictionary[16]);
            model.ElementsDictionary.Add(element7.ID, element7);
            model.SubdomainsDictionary[0].Elements.Add(element7);
            #endregion

            #region element8
            var element8 = new Element()
            {
                ID = 8,
                ElementType = new Quad4(material)
            };
            element8.AddNode(model.NodesDictionary[15]);
            element8.AddNode(model.NodesDictionary[16]);
            element8.AddNode(model.NodesDictionary[19]);
            element8.AddNode(model.NodesDictionary[18]);
            model.ElementsDictionary.Add(element8.ID, element8);
            model.SubdomainsDictionary[0].Elements.Add(element8);
            #endregion

            #region element9
            var element9 = new Element()
            {
                ID = 9,
                ElementType = new Quad4(material)
            };
            element9.AddNode(model.NodesDictionary[18]);
            element9.AddNode(model.NodesDictionary[19]);
            element9.AddNode(model.NodesDictionary[21]);
            element9.AddNode(model.NodesDictionary[20]);
            model.ElementsDictionary.Add(element9.ID, element9);
            model.SubdomainsDictionary[0].Elements.Add(element9);
            #endregion

            #region element10
            var element10 = new Element()
            {
                ID = 10,
                ElementType = new Quad4(material)
            };
            element10.AddNode(model.NodesDictionary[20]);
            element10.AddNode(model.NodesDictionary[21]);
            element10.AddNode(model.NodesDictionary[23]);
            element10.AddNode(model.NodesDictionary[22]);
            model.ElementsDictionary.Add(element10.ID, element10);
            model.SubdomainsDictionary[0].Elements.Add(element10);
            #endregion

            #region element11
            var element11 = new Element()
            {
                ID = 11,
                ElementType = new Quad4(material)
            };
            element11.AddNode(model.NodesDictionary[22]);
            element11.AddNode(model.NodesDictionary[23]);
            element11.AddNode(model.NodesDictionary[25]);
            element11.AddNode(model.NodesDictionary[24]);
            model.ElementsDictionary.Add(element11.ID, element11);
            model.SubdomainsDictionary[0].Elements.Add(element11);
            #endregion

            #region element12
            var element12 = new Element()
            {
                ID = 12,
                ElementType = new Quad4(material)
            };
            element12.AddNode(model.NodesDictionary[24]);
            element12.AddNode(model.NodesDictionary[25]);
            element12.AddNode(model.NodesDictionary[27]);
            element12.AddNode(model.NodesDictionary[26]);
            model.ElementsDictionary.Add(element12.ID, element12);
            model.SubdomainsDictionary[0].Elements.Add(element12);
            #endregion

            #region element13
            var element13 = new Element()
            {
                ID = 13,
                ElementType = new Quad4(material)
            };
            element13.AddNode(model.NodesDictionary[26]);
            element13.AddNode(model.NodesDictionary[27]);
            element13.AddNode(model.NodesDictionary[29]);
            element13.AddNode(model.NodesDictionary[28]);
            model.ElementsDictionary.Add(element13.ID, element13);
            model.SubdomainsDictionary[0].Elements.Add(element13);
            #endregion
            #endregion

            #region Constrains
            var constrainedNodes = new int[] { 0, 1, 3, 4, 5, 6, 7, 8, 9 };
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            }

            #endregion

            #region Loads
            //ground vertical loads
            model.Loads.Add(new Load() { Amount = -25800, Node = model.NodesDictionary[10], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[11], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[12], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[13], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[14], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[15], DOF = StructuralDof.TranslationY });

            //ground horizontal loads
            model.Loads.Add(new Load() { Amount = -2130, Node = model.NodesDictionary[28], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -11490, Node = model.NodesDictionary[26], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -20990, Node = model.NodesDictionary[24], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -30790, Node = model.NodesDictionary[22], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -40600, Node = model.NodesDictionary[22], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -25800, Node = model.NodesDictionary[20], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -50390, Node = model.NodesDictionary[18], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -28460, Node = model.NodesDictionary[15], DOF = StructuralDof.TranslationX });
            #endregion

            #endregion

            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        private static void SolveDynamicQuadRetainingWall()
        {
            #region Model
            double youngModulus = 2.1e09;
            double poissonRatio = 0.3;
            double density = 20;

            var material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };
            Model model = new Model();

            model.SubdomainsDictionary.Add(0, new Subdomain(0));

            #region Nodes
            model.NodesDictionary.Add(0, new Node(id: 0, x: 0.0, y: 0.0, z: 0.0));
            model.NodesDictionary.Add(1, new Node(id: 1, x: 0.6, y: 0.0, z: 0.0));

            model.NodesDictionary.Add(2, new Node(id: 2, x: 0.0, y: 0.5, z: 0.0));
            model.NodesDictionary.Add(3, new Node(id: 3, x: 0.6, y: 0.5, z: 0.0));
            model.NodesDictionary.Add(4, new Node(id: 4, x: 1.2, y: 0.5, z: 0.0));
            model.NodesDictionary.Add(5, new Node(id: 5, x: 1.8, y: 0.5, z: 0.0));
            model.NodesDictionary.Add(6, new Node(id: 6, x: 2.4, y: 0.5, z: 0.0));
            model.NodesDictionary.Add(7, new Node(id: 7, x: 3.0, y: 0.5, z: 0.0));
            model.NodesDictionary.Add(8, new Node(id: 8, x: 3.6, y: 0.5, z: 0.0));
            model.NodesDictionary.Add(9, new Node(id: 9, x: 4.3, y: 0.5, z: 0.0));

            model.NodesDictionary.Add(10, new Node(id: 10, x: 0.0, y: 1.2, z: 0.0));
            model.NodesDictionary.Add(11, new Node(id: 11, x: 0.6, y: 1.2, z: 0.0));
            model.NodesDictionary.Add(12, new Node(id: 12, x: 1.2, y: 1.2, z: 0.0));
            model.NodesDictionary.Add(13, new Node(id: 13, x: 1.8, y: 1.2, z: 0.0));
            model.NodesDictionary.Add(14, new Node(id: 14, x: 2.4, y: 1.2, z: 0.0));
            model.NodesDictionary.Add(15, new Node(id: 15, x: 3.0, y: 1.2, z: 0.0));
            model.NodesDictionary.Add(16, new Node(id: 16, x: 3.6, y: 1.2, z: 0.0));
            model.NodesDictionary.Add(17, new Node(id: 17, x: 4.3, y: 1.2, z: 0.0));

            model.NodesDictionary.Add(18, new Node(id: 18, x: 3, y: 1.9, z: 0.0));
            model.NodesDictionary.Add(19, new Node(id: 19, x: 3.5756, y: 1.9, z: 0.0));

            model.NodesDictionary.Add(20, new Node(id: 20, x: 3, y: 2.6, z: 0.0));
            model.NodesDictionary.Add(21, new Node(id: 21, x: 3.5512, y: 2.6, z: 0.0));

            model.NodesDictionary.Add(22, new Node(id: 22, x: 3, y: 3.3, z: 0.0));
            model.NodesDictionary.Add(23, new Node(id: 23, x: 3.5267, y: 3.3, z: 0.0));

            model.NodesDictionary.Add(24, new Node(id: 24, x: 3, y: 4.0, z: 0.0));
            model.NodesDictionary.Add(25, new Node(id: 25, x: 3.5023, y: 4.0, z: 0.0));

            model.NodesDictionary.Add(26, new Node(id: 26, x: 3, y: 4.7, z: 0.0));
            model.NodesDictionary.Add(27, new Node(id: 27, x: 3.4779, y: 4.7, z: 0.0));

            model.NodesDictionary.Add(28, new Node(id: 28, x: 3, y: 5.5, z: 0.0));
            model.NodesDictionary.Add(29, new Node(id: 29, x: 3.45, y: 5.5, z: 0.0));
            #endregion

            #region Elements

            #region element0
            var element0 = new Element()
            {
                ID = 0,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element0.AddNode(model.NodesDictionary[0]);
            element0.AddNode(model.NodesDictionary[1]);
            element0.AddNode(model.NodesDictionary[3]);
            element0.AddNode(model.NodesDictionary[2]);
            model.ElementsDictionary.Add(element0.ID, element0);
            model.SubdomainsDictionary[0].Elements.Add(element0);
            #endregion

            #region element1
            var element1 = new Element()
            {
                ID = 1,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element1.AddNode(model.NodesDictionary[2]);
            element1.AddNode(model.NodesDictionary[3]);
            element1.AddNode(model.NodesDictionary[11]);
            element1.AddNode(model.NodesDictionary[10]);
            model.ElementsDictionary.Add(element1.ID, element1);
            model.SubdomainsDictionary[0].Elements.Add(element1);
            #endregion

            #region element2
            var element2 = new Element()
            {
                ID = 2,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element2.AddNode(model.NodesDictionary[3]);
            element2.AddNode(model.NodesDictionary[4]);
            element2.AddNode(model.NodesDictionary[12]);
            element2.AddNode(model.NodesDictionary[11]);
            model.ElementsDictionary.Add(element2.ID, element2);
            model.SubdomainsDictionary[0].Elements.Add(element2);
            #endregion

            #region element3
            var element3 = new Element()
            {
                ID = 3,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element3.AddNode(model.NodesDictionary[4]);
            element3.AddNode(model.NodesDictionary[5]);
            element3.AddNode(model.NodesDictionary[13]);
            element3.AddNode(model.NodesDictionary[12]);
            model.ElementsDictionary.Add(element3.ID, element3);
            model.SubdomainsDictionary[0].Elements.Add(element3);
            #endregion

            #region element4
            var element4 = new Element()
            {
                ID = 4,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element4.AddNode(model.NodesDictionary[5]);
            element4.AddNode(model.NodesDictionary[6]);
            element4.AddNode(model.NodesDictionary[14]);
            element4.AddNode(model.NodesDictionary[13]);
            model.ElementsDictionary.Add(element4.ID, element4);
            model.SubdomainsDictionary[0].Elements.Add(element4);
            #endregion

            #region element5
            var element5 = new Element()
            {
                ID = 5,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element5.AddNode(model.NodesDictionary[6]);
            element5.AddNode(model.NodesDictionary[7]);
            element5.AddNode(model.NodesDictionary[15]);
            element5.AddNode(model.NodesDictionary[14]);
            model.ElementsDictionary.Add(element5.ID, element5);
            model.SubdomainsDictionary[0].Elements.Add(element5);
            #endregion

            #region element6
            var element6 = new Element()
            {
                ID = 6,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element6.AddNode(model.NodesDictionary[7]);
            element6.AddNode(model.NodesDictionary[8]);
            element6.AddNode(model.NodesDictionary[16]);
            element6.AddNode(model.NodesDictionary[15]);
            model.ElementsDictionary.Add(element6.ID, element6);
            model.SubdomainsDictionary[0].Elements.Add(element6);
            #endregion

            #region element7
            var element7 = new Element()
            {
                ID = 7,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element7.AddNode(model.NodesDictionary[8]);
            element7.AddNode(model.NodesDictionary[9]);
            element7.AddNode(model.NodesDictionary[17]);
            element7.AddNode(model.NodesDictionary[16]);
            model.ElementsDictionary.Add(element7.ID, element7);
            model.SubdomainsDictionary[0].Elements.Add(element7);
            #endregion

            #region element8
            var element8 = new Element()
            {
                ID = 8,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element8.AddNode(model.NodesDictionary[15]);
            element8.AddNode(model.NodesDictionary[16]);
            element8.AddNode(model.NodesDictionary[19]);
            element8.AddNode(model.NodesDictionary[18]);
            model.ElementsDictionary.Add(element8.ID, element8);
            model.SubdomainsDictionary[0].Elements.Add(element8);
            #endregion

            #region element9
            var element9 = new Element()
            {
                ID = 9,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element9.AddNode(model.NodesDictionary[18]);
            element9.AddNode(model.NodesDictionary[19]);
            element9.AddNode(model.NodesDictionary[21]);
            element9.AddNode(model.NodesDictionary[20]);
            model.ElementsDictionary.Add(element9.ID, element9);
            model.SubdomainsDictionary[0].Elements.Add(element9);
            #endregion

            #region element10
            var element10 = new Element()
            {
                ID = 10,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element10.AddNode(model.NodesDictionary[20]);
            element10.AddNode(model.NodesDictionary[21]);
            element10.AddNode(model.NodesDictionary[23]);
            element10.AddNode(model.NodesDictionary[22]);
            model.ElementsDictionary.Add(element10.ID, element10);
            model.SubdomainsDictionary[0].Elements.Add(element10);
            #endregion

            #region element11
            var element11 = new Element()
            {
                ID = 11,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element11.AddNode(model.NodesDictionary[22]);
            element11.AddNode(model.NodesDictionary[23]);
            element11.AddNode(model.NodesDictionary[25]);
            element11.AddNode(model.NodesDictionary[24]);
            model.ElementsDictionary.Add(element11.ID, element11);
            model.SubdomainsDictionary[0].Elements.Add(element11);
            #endregion

            #region element12
            var element12 = new Element()
            {
                ID = 12,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element12.AddNode(model.NodesDictionary[24]);
            element12.AddNode(model.NodesDictionary[25]);
            element12.AddNode(model.NodesDictionary[27]);
            element12.AddNode(model.NodesDictionary[26]);
            model.ElementsDictionary.Add(element12.ID, element12);
            model.SubdomainsDictionary[0].Elements.Add(element12);
            #endregion

            #region element13
            var element13 = new Element()
            {
                ID = 13,
                ElementType = new Quad4(material)
                {
                    Density = density,
                    RayleighAlpha = 0.05,
                    RayleighBeta = 0.05
                }
            };
            element13.AddNode(model.NodesDictionary[26]);
            element13.AddNode(model.NodesDictionary[27]);
            element13.AddNode(model.NodesDictionary[29]);
            element13.AddNode(model.NodesDictionary[28]);
            model.ElementsDictionary.Add(element13.ID, element13);
            model.SubdomainsDictionary[0].Elements.Add(element13);
            #endregion
            #endregion

            #region Constrains
            var constrainedNodes = new int[] { 0, 1, 3, 4, 5, 6, 7, 8, 9 };
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                model.NodesDictionary[constrainedNodes[i]].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            }

            #endregion

            #region Loads
            //ground vertical loads
            model.Loads.Add(new Load() { Amount = -25800, Node = model.NodesDictionary[10], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[11], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[12], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[13], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[14], DOF = StructuralDof.TranslationY });
            model.Loads.Add(new Load() { Amount = -25800 * 2, Node = model.NodesDictionary[15], DOF = StructuralDof.TranslationY });

            //ground horizontal loads
            model.Loads.Add(new Load() { Amount = -2130, Node = model.NodesDictionary[28], DOF = StructuralDof.TranslationX  });
            model.Loads.Add(new Load() { Amount = -11490, Node = model.NodesDictionary[26], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -20990, Node = model.NodesDictionary[24], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -30790, Node = model.NodesDictionary[22], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -40600, Node = model.NodesDictionary[22], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -25800, Node = model.NodesDictionary[20], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -50390, Node = model.NodesDictionary[18], DOF = StructuralDof.TranslationX });
            model.Loads.Add(new Load() { Amount = -28460, Node = model.NodesDictionary[15], DOF = StructuralDof.TranslationX });

            model.MassAccelerationHistoryLoads.Add(new MassAccelerationHistoryLoad("..\\..\\..\\elcentro_NS.dat", 1) { DOF = StructuralDof.TranslationX });
            #endregion

            #endregion

            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear dynamic analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.02, 53.74);
            //parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
            parentAnalyzerBuilder.SetNewmarkParameters(0.6, 1.0);
            NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }

        private static void SolveStaticLinearWall()
        {
            #region Read Data
            string workingDirectory = @"C:\Users\Dimitris\Desktop\Presentation";
            string meshFileName = "wall.msh";

            (IReadOnlyList<Node> nodes, IReadOnlyList<CellConnectivity<Node>> elements) = 
                GenerateMeshFromGmsh(workingDirectory, meshFileName);
            #endregion

            #region CreateModel

            const double height = 3.5;
            const double thickness = 0.1;
            const double youngModulus = 2E6;
            const double poissonRatio = 0.3;
            const double maxLoad = 1000.0;
            // Initialize
            int numberOfNodes = nodes.Count;
            int numberOfElements = elements.Count;

            // Materials
            var material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };

            // Subdomains
            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain(0));

            // Nodes
            for (int i = 0; i < numberOfNodes; ++i) model.NodesDictionary.Add(i, nodes[i]);

            // Elements
            var factory = new ContinuumElement2DFactory(thickness, material, null);
            for (int i = 0; i < numberOfElements; ++i)
            {
                ContinuumElement2D element = factory.CreateElement(elements[i].CellType, elements[i].Vertices);
                var elementWrapper = new Element() { ID = i, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(i, elementWrapper);
                model.SubdomainsDictionary[0].Elements.Add(elementWrapper);
            }

            // Constraints
            double tol = 1E-10;
            Node[] constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol).ToArray();
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            }

            // Loads
            Node[] loadedNodes = nodes.Where(
                node => (Math.Abs(node.Y - height) <= tol) && ((Math.Abs(node.X) <= tol))).ToArray();
            if (loadedNodes.Length != 1) throw new Exception("Only 1 node was expected at the top left corner");
            model.Loads.Add(new Load() { Amount = maxLoad, Node = loadedNodes[0], DOF = StructuralDof.TranslationX });


            #endregion

            #region Analysis
            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Logging displacement, strain, and stress fields.
            string outputDirectory = workingDirectory + "\\Plots";
            childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = true,
                LogStrains = true,
                LogStresses = true
            };

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion
        }

        private static void SolveStaticNonLinearWall()
        {
            #region Read Data
            string workingDirectory = @"C:\Users\Dimitris\Desktop\Presentation";
            string meshFileName = "wall.msh";

            (IReadOnlyList<Node> nodes, IReadOnlyList<CellConnectivity<Node>> elements) = GenerateMeshFromGmsh(workingDirectory, meshFileName);
            #endregion

            #region CreateModel

            const double height = 3.5;
            const double thickness = 0.1;
            const double youngModulus = 2E6;
            const double poissonRatio = 0.3;
            const double maxLoad = 1000.0;
            // Initialize
            int numberOfNodes = nodes.Count;
            int numberOfElements = elements.Count;

            // Materials
            var material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };



            // Subdomains
            Model model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain(0));

            // Nodes
            for (int i = 0; i < numberOfNodes; ++i) model.NodesDictionary.Add(i, nodes[i]);

            // Elements
            var factory = new ContinuumElement2DFactory(thickness, material, null);
            for (int i = 0; i < numberOfElements; ++i)
            {
                ContinuumElement2D element = factory.CreateElement(elements[i].CellType, elements[i].Vertices);
                var elementWrapper = new Element() { ID = i, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(i, elementWrapper);
                model.SubdomainsDictionary[0].Elements.Add(elementWrapper);
            }

            // Constraints
            double tol = 1E-10;
            Node[] constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol).ToArray();
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            }

            // Loads
            Node[] loadedNodes = nodes.Where(
                node => (Math.Abs(node.Y - height) <= tol) && ((Math.Abs(node.X) <= tol))).ToArray();
            if (loadedNodes.Length != 1) throw new Exception("Only 1 node was expected at the top left corner");
            model.Loads.Add(new Load() { Amount = maxLoad, Node = loadedNodes[0], DOF = StructuralDof.TranslationX });


            #endregion

            #region Analysis
            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Logging displacement, strain, and stress fields.
            string outputDirectory = workingDirectory + "\\Plots";
            childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = true,
                LogStrains = true,
                LogStresses = true
            };

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion
        }

        private static void SolveDynamicWall()
        {
            #region Read Data
            string workingDirectory = @"C:\Users\Dimitris\Desktop\Presentation";
            string meshFileName = "wall.msh";

            (IReadOnlyList<Node> nodes, IReadOnlyList<CellConnectivity<Node>> elements) = 
                GenerateMeshFromGmsh(workingDirectory, meshFileName);
            #endregion

            #region CreateModel

            const double height = 3.5;
            const double thickness = 0.1;
            const double youngModulus = 2E6;
            const double poissonRatio = 0.3;
            const double maxLoad = 1000.0;
            // Initialize
            int numberOfNodes = nodes.Count;
            int numberOfElements = elements.Count;

            // Materials
            var material = new ElasticMaterial2D(StressState2D.PlaneStress)
            {
                YoungModulus = youngModulus,
                PoissonRatio = poissonRatio
            };

            var dynamicMaterial = new DynamicMaterial(25, 0.05, 0.05, true);

            // Subdomains
            var model = new Model();
            model.SubdomainsDictionary.Add(0, new Subdomain(0));

            // Nodes
            for (int i = 0; i < numberOfNodes; ++i) model.NodesDictionary.Add(i, nodes[i]);

            // Elements
            var factory = new ContinuumElement2DFactory(thickness, material, dynamicMaterial);
            for (int i = 0; i < numberOfElements; ++i)
            {
                ContinuumElement2D element = factory.CreateElement(elements[i].CellType, elements[i].Vertices);
                var elementWrapper = new Element() { ID = i, ElementType = element };
                foreach (Node node in element.Nodes) elementWrapper.AddNode(node);
                model.ElementsDictionary.Add(i, elementWrapper);
                model.SubdomainsDictionary[0].Elements.Add(elementWrapper);
            }

            // Constraints
            double tol = 1E-10;
            Node[] constrainedNodes = nodes.Where(node => Math.Abs(node.Y) <= tol).ToArray();
            for (int i = 0; i < constrainedNodes.Length; i++)
            {
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
                constrainedNodes[i].Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            }

            // Loads
            Node[] loadedNodes = nodes.Where(
                node => (Math.Abs(node.Y - height) <= tol) && ((Math.Abs(node.X) <= tol))).ToArray();
            if (loadedNodes.Length != 1) throw new Exception("Only 1 node was expected at the top left corner");
            model.Loads.Add(new Load() { Amount = maxLoad, Node = loadedNodes[0], DOF = StructuralDof.TranslationX });

            model.MassAccelerationHistoryLoads.Add(new MassAccelerationHistoryLoad("..\\..\\..\\elcentro_NS.dat", 1) { DOF = StructuralDof.TranslationX });


            #endregion

            #region Analysis
            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear dynamic analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzerBuilder = new NewmarkDynamicAnalyzer.Builder(model, solver, provider, childAnalyzer, 0.02, 53.74);
            //parentAnalyzerBuilder.SetNewmarkParametersForConstantAcceleration(); // Not necessary. This is the default
            parentAnalyzerBuilder.SetNewmarkParameters(0.6, 1.0);
            NewmarkDynamicAnalyzer parentAnalyzer = parentAnalyzerBuilder.Build();

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();

            // Logging displacement, strain, and stress fields.
            string outputDirectory = workingDirectory + "\\Plots";
            childAnalyzer.LogFactories[0] = new VtkLogFactory(model, outputDirectory)
            {
                LogDisplacements = true,
                LogStrains = true,
                LogStresses = true
            };

            // Run the analysis
            parentAnalyzer.BuildMatrices();
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
            #endregion
        }

        private static (IReadOnlyList<Node> nodes, IReadOnlyList<CellConnectivity<Node>> elements) GenerateMeshFromGmsh(
            string workingDirectory, string meshFileName)
        {
            using (var reader = new GmshReader<Node>(workingDirectory + "\\" + meshFileName))
            {
                return reader.CreateMesh((id, x, y, z) => new Node(id: id, x: x, y:  y, z: z ));
            }
        }

        public void Solve4QuadProblem()
        {

            #region Nodes
            var node0 = new Node( id: 0, x: 0.0, y:  0.0 );
            var node1 = new Node( id: 1, x: 2.0, y:  0.0 );
            var node2 = new Node( id: 2, x: 4.0, y:  0.0 );
            var node3 = new Node( id: 3, x: 0.0, y:  2.0 );
            var node4 = new Node( id: 4, x: 2.0, y:  2.0 );
            var node5 = new Node( id: 5, x: 4.0, y:  2.0 );
            var node6 = new Node( id: 6, x: 0.0, y:  4.0 );
            var node7 = new Node( id: 7, x: 2.0, y:  4.0 );
            var node8 = new Node( id: 8, x: 4.0, y:  4.0 );
            #endregion

            #region material
            var material = new ElasticMaterial2D(StressState2D.PlaneStress);
            material.YoungModulus = 3e7;
            material.PoissonRatio = 0.2;

            double thickness = 0.3;
            #endregion

            #region elements
            var element0 = new Element { ID = 0, ElementType = new Quad4(material) { Thickness = thickness }};
            element0.AddNodes(new[] { node0, node1, node4, node3 });

            var element1 = new Element { ID = 1, ElementType = new Quad4(material) { Thickness = thickness }};
            element1.AddNodes(new[] { node1, node2, node5, node4 });

            var element2 = new Element { ID = 2, ElementType = new Quad4(material) { Thickness = thickness }};
            element2.AddNodes(new[] { node3, node4, node7, node6 });

            var element3 = new Element { ID = 3, ElementType = new Quad4(material) { Thickness = thickness }};
            element3.AddNodes(new[] { node4, node5, node8, node7 });
            #endregion

            #region loads
            var load0 = new Load { Amount = 50, DOF = StructuralDof.TranslationX, Node = node8 };
            #endregion

            #region constraints
            node0.Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            node0.Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });

            node2.Constraints.Add(new Constraint { DOF = StructuralDof.TranslationX });
            node2.Constraints.Add(new Constraint { DOF = StructuralDof.TranslationY });
            #endregion

            #region subdomains
            var subdomain0 = new Subdomain(0);
            subdomain0.Elements.Add(element0);
            subdomain0.Elements.Add(element1);
            subdomain0.Elements.Add(element2);
            subdomain0.Elements.Add(element3);
            #endregion

            #region model
            var model = new Model();

            model.NodesDictionary.Add(0, node0);
            model.NodesDictionary.Add(1, node1);
            model.NodesDictionary.Add(2, node2);
            model.NodesDictionary.Add(3, node3);
            model.NodesDictionary.Add(4, node4);
            model.NodesDictionary.Add(5, node5);
            model.NodesDictionary.Add(6, node6);
            model.NodesDictionary.Add(7, node7);
            model.NodesDictionary.Add(8, node8);

            model.ElementsDictionary.Add(0, element0);
            model.ElementsDictionary.Add(1, element1);
            model.ElementsDictionary.Add(2, element2);
            model.ElementsDictionary.Add(3, element3);

            model.Loads.Add(load0);

            model.SubdomainsDictionary.Add(0, subdomain0);

            #endregion

            // Choose linear equation system solver
            var solverBuilder = new SkylineSolver.Builder();
            ISolver solver = solverBuilder.BuildSolver(model);

            // Structural problem provider
            var provider = new ProblemStructural(model, solver);

            // Linear static analysis
            var childAnalyzer = new LinearAnalyzer(model, solver, provider);
            var parentAnalyzer = new StaticAnalyzer(model, solver, provider, childAnalyzer);

            // Run the analysis
            parentAnalyzer.Initialize();
            parentAnalyzer.Solve();
        }
    }
}
