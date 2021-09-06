using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Analyzers.Dynamic;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Problems;
using ISAAR.MSolve.Solvers.Direct;
using Xunit;
using ISAAR.MSolve.Discretization.Interfaces;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using System.Linq;
using ISAAR.MSolve.Solvers;
using System.IO;
using ISAAR.MSolve.FEM.Readers;
using ISAAR.MSolve.FEM.Readers.Interfaces;
using ISAAR.MSolve.FEM.Loading.SurfaceLoads;
using static ISAAR.MSolve.FEM.Loading.SurfaceLoads.WeakDirichlet;
using ISAAR.MSolve.FEM.Loading;
using ISAAR.MSolve.FEM.Elements.BoundaryConditionElements;
using System;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISSAR.MSolve.Discretization.Loads;
using ISAAR.MSolve.FEM.Loading.BodyLoads;
using System.Reflection;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.FEM.Loading.Interfaces;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Tests
{
	public class tumModel_u_lg_c_ox_Tumc_p_Csdf_endoc_Cvegf_a1_a2_phis_ΝΚ_CD8_Treg_
	{
		private const double timestep = .2;
		private const double time = 14;
		private const int subdomainID = 0;
		private static readonly double[] loxc = new double[] { 2.1 / 24d / 3600d, 1.54 / 24d / 3600d }; //1/s
		private static readonly double[] Aox = new double[] { 2200.0d / 24d / 3600d, 2200d / 24d / 3600d }; //mol/(m^3*s)
		private static readonly double[] Dox = new double[] { 1.78e-9, 1.79e-9 }; //m^2/s
		private static readonly double[] kox = new double[] { .00464, .00464 }; //mol/m^3
		private static readonly double[] Koxc = new double[] { 0.0083, 0.0083 }; //mol/m^3
		private static readonly double[] Dcell = new double[] { 5.4e-3, 1.8e-4 }; //m^2/s
		private static readonly double[] Dvegf = new double[] { 3.1e-11, 3.1e-11 }; //m^2/s
		private static readonly double[] khy = new double[] { 6.5e-11, 6.5e-11 };  //m^2/(Pa*d);
		private static readonly double[] lp = new double[] { 2.7e-12, 2.7e-12 };  //m^2*s/kg;
		private static readonly double lpSv = 3.75e-4;  //m*s/kg;
		private static double pc0 = 3.32e-11; //mol/m^3
		private static double[] pc = new double[] { 2 * 3.32e-11, 3.32e-11 }; //mol/m^3
		private static double T0 = 1; //kg/m^3
		private static double Cv0 = 1; //kg/m^3
		private static double Sfn = 1; //kg/m^3
		private static double ptc = 0.55 / 24d / 3600d; //1/s
		private static double pti = 0.21 / 24d / 3600d; //kg/m^3
		private static double lm1tum = 3d / 24d / 3600d; //1/s
		private static double Cs0 = 1; //kg/m^3
		private static double l2 = 1e-5; //m^3/kg/s
		private static double l4 = 1e-11; //m/s
		private static double l10 = 6.8e-9; //1/s
		private static double l11 = 4e-8; //m/s
		private static double l13 = 4e-9; //1/s
		private static double cvox = 0.2; //mol/m^3
		private static double Lwv = 5e-6; //m
		private static double mtox = 8e-3 * 1.1 * 1e-6 / 3600d; //m^2/s
		private static double pv = 4000; //kg/m/s^2
		private static double r0 = 4e-7; //m
		private static double Cvegf0 = 1; //kg/m^3
		private static double[] Svin = new double[] { 7000 - (7000 / 250e-9) * (r0 - 400e-9), 7000 }; //1/m
		private static double WvsTc = 1;
		private static double WvsSv = 1;
		private static double xn = 1e-12;//m^5/kg/s
		private static double Dsv = 1e-7 * 1e-8;//m^2/s
		private static double s1 = 1e6;//m^3/mol
		private static double s2 = 1e6;//m^3/mol
		private static double a10 = 1e-6;//mol/m^3
		private static double a20 = 1e-6;//mol/m^3
		private static double aD = 1;
		private static double bD = 1;
		private static double m1 = 4.56 * 100d / 3600d; //1/s
		private static double m2 = 4.56 * 100d / 3600d; //1/s
		private static double b1 = 2280d / 3600d; //1/s
		private static double b2 = 18240d / 3600d; //1/s
		private static double T0in = 500;
		public static double kec = 0;
		public static double canti = 0;
		private static double kantivegf = 0;
		private static double snkcon = 1.3e4 / T0in / 24d / 3600d; // 1/s
		private static double fdrnk0 = 0.0412 / 24d / 3600d; // 1/s
		private static double grrnk0 = 0.025 / 24d / 3600d; //1/s
		private static double hscnk = 2.02e7 / T0in / T0in;
		private static double pirnk = 1e-7 * T0in / 24d / 3600d; //1/s
		private static double lreg = 100 / 24d / 3600d; //1/s
		private static double lm2Nk = 7 / 24d / 3600d; //1/s
		private static double cnk = 3.23e-7 * T0in / 24d / 3600d; // 1/s
		private static double scd8 = 0.01 / 24d / 3600d; // 1/s
		private static double rst8 = (1.1e-7 / 24d / 3600d) * T0in; // 1/s
		private static double mdrt8 = 0.02 / 24d / 3600d; // 1/s
		private static double jrrt80 = 0.0375 / 24d / 3600d; // 1/s
		private static double dt8 = 1.43 / 24d / 3600d; // 1/s
		private static double limt8 = 1.36;
		private static double st8 = 2.73; 
		private static double qirt8= 3.42e-10 * T0in / 24d / 3600d; // 1/s
		private static double lm2CD8 = (-0.14286 * 0.01 + 1.3214) / 24d / 3600d; // 1/s
		private static double ksct8 = (2e7) / (T0in * T0in);
		private static double mtreg = 0.02 / 24d / 3600d; // 1/s
		private static double gtreg = 0.0375 / 24d / 3600d; // 1/s
		private static double[][] conv0 = new double[][] { new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 } };
		private static int solverSymmetric = 0, solverNonSymmetric = 0;
		private static bool reordering = false;
		private static ISolverBuilder builder, structuralBuilder;
		//private static ISolverBuilder builder, asymBuilder, structuralBuilder;
		private static double[] lgNode;
		private static double[] lgElement;
		private static double[] CsdfNode;
		private static double[] CsdfElement;
		private static double[] CvegfNode;
		private static double[] CvegfElement;
		private static double[] coxNode;
		private static double[] coxElement;
		private static double[] tumcNode;
		private static double[] pSolution;
		private static double[] tumcElement;
		private static double[] pNode;
		private static double[] pElement;
		private static double[] endocNode;
		private static double[] endocElement;
		private static double[] a1Node;
		private static double[] a1Element;
		private static double[] a2Node;
		private static double[] a2Element;
		private static double[] phisNode;
		private static double[] phisElement;
		private static double[] NKcellsNode;
		private static double[] NKcellsElement;
		private static double[] T8imcellsNode;
		private static double[] T8imcellsElement;
		private static double[] TregNode;
		private static double[] TregElement;
		private static double[] dd0;
		private static double[][] vNode;
		private static Dictionary<int, double[]> vElement;/*= new Dictionary<int, double[]>();*/
		private static double[][] uNode;
		private static double[][] prevUNode;
		private static Dictionary<int, double[]> uElement = new Dictionary<int, double[]>();
		private static double[][] dcoxdx; /*= new Dictionary<int, double[]>();*/
		private static double[][] dpdx; /*= new Dictionary<int, double[]>();*/
		private static double[][] dphisdx; /*= new Dictionary<int, double[]>();*/
		private static Dictionary<int, double> OxygenTransportK;
		private static Dictionary<int, double[]> OxygenTransportU;
		private static Dictionary<int, double> OxygenTransportL;
		private static Dictionary<int, double[]> CancerTransportU;
		private static Dictionary<int, double> CancerTransportL;
		private static Dictionary<int, double[]> lgU;
		private static Dictionary<int, double> lgL;
		private static Dictionary<int, double[]> PressureU;
		private static Dictionary<int, double> PressureL;
		private static Dictionary<int, double> CvegfK;
		private static Dictionary<int, double[]> CvegfU;
		private static Dictionary<int, double> CvegfL;
		private static Dictionary<int, double[]> phisU;
		private static Dictionary<int, double> phisL;
		private static Dictionary<int, double[]> NKcellsU;
		private static Dictionary<int, double> NKcellsL;
		private static Dictionary<int, double[]> CD8U;
		private static Dictionary<int, double> CD8L;
		private static Dictionary<int, double[]> TregU;
		private static double TregL = (mtreg - gtreg) * 24d * 3600d;
		private static double[][] PreviousStrains;
		private static double[][] Strains;
		private static double[][] Stresses;
		private static Dictionary<int, double[]> uXt;
		private static Dictionary<int, IVector> Displacements;
		private static Tuple<Model, IModelReader> oxModel, gModel, ctModel, prModel, csdfModel,
			endocModel, cvegfModel, a1Model, a2Model, phisModel, NKcellsModel, CD8Model, TregModel, structModel;
		private static int pressureModelFreeDOFs = 0;
		private static string inputFile = "mesh446elem.mphtxt";
		private static int NewtonRaphsonIncrements = 5;
		private static int NewtonRaphosnIterations = 10;
		private static double NewtonRaphsonTolerarance = 1e-3;
		private static int NewtonRaphsonIterForMatrixRebuild = 5;

		static tumModel_u_lg_c_ox_Tumc_p_Csdf_endoc_Cvegf_a1_a2_phis_ΝΚ_CD8_Treg_()
		{
			if (solverSymmetric == 0)
			{
				builder = new SkylineSolver.Builder();
				structuralBuilder = new SkylineSolver.Builder();
			}
			else if (solverSymmetric == 1)
			{
				IDofReorderingStrategy reorderingStrategy;
				if (reordering)
				{
					reorderingStrategy = AmdReordering.CreateWithCSparseAmd();
				}
				else
				{
					reorderingStrategy = new NullReordering();
				}
				builder = new CSparseCholeskySolver.Builder()
				{
					DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), reorderingStrategy)
				};
				structuralBuilder = new CSparseCholeskySolver.Builder()
				{
					DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), reorderingStrategy)
				};
			}
			else
			{
				IDofReorderingStrategy reorderingStrategy;
				if (reordering)
				{
					reorderingStrategy = AmdReordering.CreateWithSuiteSparseAmd();
				}
				else
				{
					reorderingStrategy = new NullReordering();
				}
				builder = new SuiteSparseSolver.Builder()
				{
					DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), reorderingStrategy)
				};
				structuralBuilder = new SuiteSparseSolver.Builder()
				{
					DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), reorderingStrategy)
				};
			}
		}
		[Fact]
		private static void RunTest()
		{
			var path1 = Path.Combine(Directory.GetCurrentDirectory(), $"solutionNorms");
			if (!Directory.Exists(path1))
			{
				Directory.CreateDirectory(path1);
			}
			var path2 = Path.Combine(path1, $"solutionNorm.txt");
			ISAAR.MSolve.Discretization.Logging.GlobalLogger.OpenOutputFile(path2);
			var DoxDays = new double[Dox.Length];
			for (int i = 0; i < Dox.Length; i++)
			{
				DoxDays[i] = 24d * 3600d * Dox[i];
			}

			var DcellDays = new double[Dcell.Length];
			for (int i = 0; i < Dox.Length; i++)
			{
				DcellDays[i] = 24d * 3600d * Dcell[i];
			}

			endocModel = CreateEndocModel();
			oxModel = CreateOxygenTransportModel(DoxDays);
			ctModel = CreateCancerTransportModel(DcellDays[0]);
			gModel = CreateGrowthModel();
			prModel = CreatePressureModel(khy);
			csdfModel = CreateCsdfModel();
			cvegfModel = CreateCvegfModel();
			a1Model = CreateAng1Model();
			a2Model = CreateAng2Model();
			phisModel = CreatePhisModel();
			NKcellsModel = CreateNKcellsModel();
			CD8Model = CreateCD8Model();
			TregModel = CreateTregModel();

			var models = new[] { endocModel.Item1, oxModel.Item1, ctModel.Item1, gModel.Item1,
				prModel.Item1, csdfModel.Item1, csdfModel.Item1, a1Model.Item1, a2Model.Item1, 
				phisModel.Item1, NKcellsModel.Item1, CD8Model.Item1, TregModel.Item1};
			var modelReaders = new[] { endocModel.Item2, oxModel.Item2, ctModel.Item2, gModel.Item2,
				prModel.Item2, csdfModel.Item2, cvegfModel.Item2, a1Model.Item2, a2Model.Item2,
				phisModel.Item2, NKcellsModel.Item2, CD8Model.Item2, TregModel.Item2};
			IVectorView[] solutions = SolveModelsWithStructuralAnalyzer(models, modelReaders);

			Assert.True(CompareResults(solutions[0]));
		}
		private static void Paraview(int timeStep)
		{
			var path0 = Path.Combine(Directory.GetCurrentDirectory(), $"paraviewOutput");
			var path3 = Path.Combine(path0, $"results{timeStep}.vtu");
			var numberOfPoints = structModel.Item1.Nodes.Count;
			var numberOfCells = structModel.Item1.Elements.Count;
			using (StreamWriter outputFile = new StreamWriter(path3))
			{
				outputFile.WriteLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
				outputFile.WriteLine("  <UnstructuredGrid>");
				outputFile.WriteLine($"     <Piece NumberOfPoints=\"{numberOfPoints}\" NumberOfCells=\"{numberOfCells}\">");
				outputFile.WriteLine("          <Points>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"position\" NumberOfComponents=\"3\" format =\"ascii\">");
				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{structModel.Item1.Nodes[i].X} {structModel.Item1.Nodes[i].Y} {structModel.Item1.Nodes[i].Z} ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("          </Points>");
				outputFile.WriteLine("          <PointData>");

				outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"node_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{i + 1}");
				outputFile.WriteLine("              </DataArray>");


				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"totalDisplacement\" NumberOfComponents=\"1\" format =\"ascii\">");
				for (int i = 0; i < numberOfPoints; i++)
				{
					double dist = Math.Sqrt(Math.Pow(oxModel.Item1.Nodes[i].X - structModel.Item1.Nodes[i].X, 2) +
						Math.Pow(oxModel.Item1.Nodes[i].Y - structModel.Item1.Nodes[i].Y, 2) +
						Math.Pow(oxModel.Item1.Nodes[i].Z - structModel.Item1.Nodes[i].Z, 2));
					outputFile.WriteLine($"{dist} ");
				}
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"oxygen\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{coxNode[i]} ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"growth factor\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < lgNode.Length; i++)
					outputFile.WriteLine($"{lgNode[i]} ");
				for (int i = lgNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"1 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"cancer cells\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < tumcNode.Length; i++)
					outputFile.WriteLine($"{tumcNode[i]} ");
				for (int i = tumcNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"pressure\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{pNode[i]} ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Csdf\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < CsdfNode.Length; i++)
					outputFile.WriteLine($"{CsdfNode[i]} ");
				for (int i = CsdfNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"endoc\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < endocNode.Length; i++)
					outputFile.WriteLine($"{endocNode[i]} ");
				for (int i = endocNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Cvegf\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < CvegfNode.Length; i++)
					outputFile.WriteLine($"{CvegfNode[i]} ");
				for (int i = CvegfNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Ang1\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < a1Node.Length; i++)
					outputFile.WriteLine($"{a1Node[i]} ");
				for (int i = a1Node.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Ang2\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < a2Node.Length; i++)
					outputFile.WriteLine($"{a2Node[i]} ");
				for (int i = a2Node.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"phis\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < phisNode.Length; i++)
					outputFile.WriteLine($"{phisNode[i]} ");
				for (int i = phisNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0.3 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"NKcells\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < NKcellsNode.Length; i++)
					outputFile.WriteLine($"{NKcellsNode[i]} ");
				for (int i = NKcellsNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"CD8cells\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < T8imcellsNode.Length; i++)
					outputFile.WriteLine($"{T8imcellsNode[i]} ");
				for (int i = T8imcellsNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Tregs\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < TregNode.Length; i++)
					outputFile.WriteLine($"{TregNode[i]} ");
				for (int i = TregNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("          </PointData>");
				outputFile.WriteLine("          <CellData>");
				outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"element_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < numberOfCells; i++)
				{
					outputFile.WriteLine($"{i + 1}");
				}
				outputFile.WriteLine("              </DataArray>");
				outputFile.WriteLine("          </CellData>");
				outputFile.WriteLine("          <Cells>");

				outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"connectivity\">");
				for (int i = 0; i < numberOfCells; i++)
				{
					for (int j = 0; j < structModel.Item1.Elements[i].Nodes.Count; j++)
						outputFile.Write($"{structModel.Item1.Elements[i].Nodes[j].ID} ");
					outputFile.WriteLine("");
				}
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">");
				var offset = 0;
				for (int i = 0; i < numberOfCells; i++)
				{
					offset += structModel.Item1.Elements[i].Nodes.Count;
					outputFile.WriteLine($"{offset} ");
				}
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Int32\" Name =\"types\" NumberOfComponents =\"1\" format=\"ascii\">");
				for (int i = 0; i < numberOfCells; i++)
				{
					if (structModel.Item1.Elements[i].Nodes.Count == 4)
						outputFile.WriteLine($"{10} ");
					else outputFile.WriteLine($"{5} ");
				}
				outputFile.WriteLine("              </DataArray>");
				outputFile.WriteLine("          </Cells>");
				outputFile.WriteLine("      </Piece>");
				outputFile.WriteLine("  </UnstructuredGrid>");
				outputFile.WriteLine("</VTKFile>");
			}
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
		private static Tuple<double[][], double[][]> GetStrainsStresses(int elementsNo)
		{
			if (structModel == null)
			{
				double[][] strains = new double[elementsNo][];
				double[][] stresses = new double[elementsNo][];
				for (int i = 0; i < elementsNo; i++)
				{
					strains[i] = new double[6];
					stresses[i] = new double[6];
				}
				return new Tuple<double[][], double[][]>(strains, stresses);
			}
			else
			{
				IList<Element> elements = structModel.Item1.Elements;
				double[][] strains = new double[elements.Count][];
				double[][] stresses = new double[elements.Count][];
				if (Displacements == null)
				{
					Displacements = new Dictionary<int, IVector>();
					Displacements.Add(0, Vector.CreateZero(structModel.Item1.GlobalDofOrdering.NumGlobalFreeDofs));
				}
				foreach (Element e in elements)
				{
					double[] localVector = e.Subdomain.FreeDofOrdering.ExtractVectorElementFromSubdomain(e, Displacements[0]);
					var strainStresses = e.ElementType.CalculateStresses(e, localVector,
						new double[e.ElementType.GetElementDofTypes(e).SelectMany(x => x).Count()]);
					strains[e.ID] = new double[strainStresses.Item1.Length];
					stresses[e.ID] = new double[strainStresses.Item2.Length];
					Array.Copy(strainStresses.Item1, strains[e.ID], strains[e.ID].Length);
					Array.Copy(strainStresses.Item2, stresses[e.ID], stresses[e.ID].Length);
				}
				return new Tuple<double[][], double[][]>(strains, stresses);
			}
		}
		private static double[][] GetCDStrains(Model model, IVector solution)
		{
			IList<Element> elements = model.Elements;
			double[][] strains = new double[elements.Count][];
			if (solution == null)
			{
				solution = Vector.CreateZero(model.GlobalDofOrdering.NumGlobalFreeDofs);
			}
			foreach (Element e in elements)
			{
				double[] localVector = e.Subdomain.FreeDofOrdering.ExtractVectorElementFromSubdomain(e, solution);
				var strainStresses = e.ElementType.CalculateStresses(e, localVector,
					new double[e.ElementType.GetElementDofTypes(e).SelectMany(x => x).Count()]);
				strains[e.ID] = new double[strainStresses.Item1.Length];
				Array.Copy(strainStresses.Item1, strains[e.ID], strains[e.ID].Length);
			}
			return strains;
		}
		private static Dictionary<int, double[]> StructuralSpaceTimeDerivatives(double[][] current, double[][] previous)
		{
			Dictionary<int, double[]> spaceTimeDerivatives = new Dictionary<int, double[]>();
			for (int i = 0; i < current.Length; i++)
			{
				spaceTimeDerivatives[i] = new double[3];
				for (int j = 0; j < 3; j++)
				{
					spaceTimeDerivatives[i][j] = (current[i][j] - previous[i][j]) / timestep;
				}
			}
			return spaceTimeDerivatives;
		}
		private static void PressureCoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			var modelReader = oxModel.Item2;
			foreach (Element element in modelReader.elementDomains[0])
			{

				l[element.ID] = dd0[element.ID] >= 0d ? lp[0] * endocElement[element.ID] * Svin[0] * dd0[element.ID] * 24d * 3600d : 0;
				u[element.ID] = conv0[0];
			}
			foreach (Element element in modelReader.elementDomains[1])
			{
				l[element.ID] = (lp[1] * Svin[1] + lpSv) * 24d * 3600d;
				u[element.ID] = conv0[1];
			}
			PressureL = l;
			PressureU = u;
		}
		private static void OxygenTransportCoefficientsCalculation(IList<IList<int>> elementIDsPerDomain, Dictionary<int, double> k, Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			if (endocElement == null)
			{
				endocElement = new double[elementIDsPerDomain[0].Count];
				foreach (var id in elementIDsPerDomain[0])
				{
					endocElement[id] = 0.5;
				}
			}
			if (tumcElement == null)
			{
				tumcElement = new double[elementIDsPerDomain[0].Count];
				for (int i = 0; i < elementIDsPerDomain[0].Count; i++)
				{
					tumcElement[i] = 0.96;
				}
			}
			if (coxElement == null)
			{
				coxElement = new double[k.Count];
				for (int i = 0; i < k.Count; i++)
				{
					coxElement[i] = 0;
				}
			}
			if (phisElement == null)
			{
				phisElement = new double[elementIDsPerDomain[0].Count];
				for (int i = 0; i < elementIDsPerDomain[0].Count; i++)
				{
					phisElement[i] = 0.3;
				}
			}
			if (dd0 == null) dd0 = new double[elementIDsPerDomain[0].Count];
			if (Stresses == null)
			{
				Stresses = new double[elementIDsPerDomain[0].Count][];
				for (int i = 0; i < elementIDsPerDomain[0].Count; i++)
				{
					Stresses[i] = new double[6];
				}
			}
			foreach (var element in endocModel.Item1.Elements)
			{
				double sumX = 0;
				double sumY = 0;
				double sumZ = 0;
				foreach (var node in element.Nodes)
				{
					sumX += node.X;
					sumY += node.Y;
					sumZ += node.Z;
				}
				double ri = Math.Sqrt(Math.Pow(sumX / element.Nodes.Count, 2)) +
					Math.Sqrt(Math.Pow(sumY / element.Nodes.Count, 2)) +
					Math.Sqrt(Math.Pow(sumZ / element.Nodes.Count, 2));
				double phi = Math.Atan2(sumY / element.Nodes.Count, sumX / element.Nodes.Count);
				double th = Math.Acos(sumZ / element.Nodes.Count / ri);
				double srr = Stresses[element.ID][2] * Math.Cos(th) * Math.Cos(th) + Stresses[element.ID][0] * Math.Cos(phi) * Math.Cos(phi) * Math.Sin(th) * Math.Sin(phi) +
					Stresses[element.ID][5] * Math.Cos(phi) * Math.Sin(2 * th) + 2 * Stresses[element.ID][3] * Math.Cos(phi) * Math.Sin(th) * Math.Sin(th) * Math.Sin(phi) +
					Stresses[element.ID][4] * Math.Sin(2 * th) * Math.Sin(phi) + Stresses[element.ID][1] * Math.Sin(th) * Math.Sin(th) * Math.Sin(phi) * Math.Sin(phi);
				double stt = Stresses[element.ID][0] * Math.Cos(th) * Math.Cos(th) * Math.Cos(phi) * Math.Cos(phi) + Stresses[element.ID][2] * Math.Sin(th) * Math.Sin(th) -
					Stresses[element.ID][5] * Math.Sin(phi) * Math.Sin(2 * th) - Stresses[element.ID][4] * Math.Sin(phi) * Math.Sin(2 * th) +
					Stresses[element.ID][1] * Math.Cos(th) * Math.Cos(th) * Math.Sin(phi) * Math.Sin(phi) + Stresses[element.ID][3] * Math.Cos(th) * Math.Cos(th) * Math.Sin(2 * phi);
				double spp = Stresses[element.ID][1] * Math.Cos(phi) * Math.Cos(phi) - 2 * Stresses[element.ID][3] * Math.Cos(phi) * Math.Sin(phi) +
					Stresses[element.ID][0] * Math.Sin(phi) * Math.Sin(phi);
				double sbulk = (srr + stt + spp) / 3d;

				//var dd = (-2e-4) * tumcElement[element.ID] / (4 * Math.PI * Math.Pow(ri, 2)) + 1;
				//dd0[element.ID] = dd >= 0 ? dd : 0;
				dd0[element.ID] = (0.0002 * Math.Pow(sbulk / 1000, 2) + 0.02618 * (sbulk / 1000) + 1);
			}
			if (uXt == null)
			{
				uXt = new Dictionary<int, double[]>();
				foreach (var domain in elementIDsPerDomain)
				{
					foreach (var id in domain)
					{
						uXt[id] = new double[3];
					}
				}
			}
			if (dpdx == null)
			{
				dpdx = new double[elementIDsPerDomain[0].Count + elementIDsPerDomain[1].Count][];
				foreach (var domain in elementIDsPerDomain)
				{
					foreach (var id in domain)
					{
						dpdx[id] = new double[3];
					}
				}
			}

			if (vElement == null)
			{
				vElement = new Dictionary<int, double[]>();
				foreach (var domain in elementIDsPerDomain)
				{
					foreach (var id in domain)
					{
						vElement[id] = new double[3];
					}
				}
			}
			foreach (var id in elementIDsPerDomain[0])
			{
				//l[e.ID] = 24 * 3600 * Dox[0] / Lwv * Svin * dd0[e.ID] * SvDElement[e.ID];
				k[id] = Dox[0] * 24d * 3600d;
				u[id][0] = (-khy[0] * dpdx[id][0] / (1 - phisElement[id])) + vElement[id][0];
				u[id][1] = (-khy[0] * dpdx[id][1] / (1 - phisElement[id])) + vElement[id][1];
				u[id][2] = (-khy[0] * dpdx[id][2] / (1 - phisElement[id])) + vElement[id][2];
				l[id] = 24d * 3600d * (Dox[0] / Lwv * Svin[0] * dd0[id] * endocElement[id]
				+ ((Aox[0]) / (kox[0] + coxElement[id])) * tumcElement[id]);
			}
			foreach (var id in elementIDsPerDomain[1])
			{
				//l[e.ID] = 24 * 3600 * Dox[1] / Lwv * Svin * SvDElement[e.ID];
				k[id] = Dox[1] * 24d * 3600d;
				u[id][0] = (-khy[1] * dpdx[id][0] / 0.7) + vElement[id][0];
				u[id][1] = (-khy[1] * dpdx[id][1] / 0.7) + vElement[id][1];
				u[id][2] = (-khy[1] * dpdx[id][2] / 0.7) + vElement[id][2];
				l[id] = 24d * 3600d * (Dox[1] / Lwv * Svin[1]);
			}
			OxygenTransportK = k;
			OxygenTransportU = u;
			OxygenTransportL = l;
		}
		private static void CvegfCoefficientsCalculation(IList<IList<int>> elementIDsPerDomain, Dictionary<int, double> k, Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			foreach (IList<int> domain in elementIDsPerDomain)
			{
				foreach (var id in domain)
				{
					k[id] = Dvegf[elementIDsPerDomain.IndexOf(domain)] * 24d * 3600d;
				}
			}
			CvegfK = k;
			CvegfU = u;
			foreach (var id in elementIDsPerDomain[0])
			{
				l[id] = 24d * 3600d * (l11 * Svin[0] * endocElement[id] + l13);
			}
			CvegfL = l;
		}
		private static void lgCoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			if (lgElement == null)
			{
				lgElement = new double[oxModel.Item1.Elements.Count];
				for (int i = 0; i < lgElement.Length; i++)
				{
					lgElement[i] = 1;
				}
			}

			for (int i = 0; i < l.Count; i++)
			{
				var Grox = (loxc[0] * coxElement[i]) / (coxElement[i] + Koxc[0]);
				var cnknew = (1d + 20d * coxElement[i] / 0.2) * cnk;
				var Rtumc = Grox * Sfn - (cnknew * NKcellsElement[i] + ptc + pti + lm1tum * 0.01) * tumcElement[i];
				l[i] = -24d * 3600d * Rtumc / 3d;
			}
			lgU = u;
			lgL = l;
		}
		private static void phisCoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			if (phisElement == null)
			{
				phisElement = new double[endocModel.Item1.Elements.Count];
				for (int i = 0; i < phisElement.Length; i++)
				{
					phisElement[i] = 0.3;
				}
			}

			if (uXt == null)
			{
				uXt = new Dictionary<int, double[]>();
				for (int i = 0; i < oxModel.Item1.Elements.Count; i++)
				{
					uXt[i] = new double[3];
				}
			}
			for (int i = 0; i < l.Count; i++)
			{
				l[i] = uXt[i].Sum();
			}
			phisU = u;
			phisL = l;
		}
		private static void CD8CoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			if (T8imcellsElement == null)
			{
				T8imcellsElement = new double[l.Count];
				for (int i = 0; i < l.Count; i++)
				{
					T8imcellsElement[i] = 0d;
				}
			}

			if (TregElement == null)
			{
				TregElement = new double[l.Count];
				for (int i = 0; i < l.Count; i++)
				{
					TregElement[i] = 0.01;
				}
			}

			for (int i = 0; i < l.Count; i++)
			{
				var dt8new = dt8 + 8 * coxElement[i] / 0.2 * dt8;
				var Dfck = (dt8new * Math.Pow(T8imcellsElement[i] / tumcElement[i], limt8) * tumcElement[i]
					/ (st8 + (Math.Pow(T8imcellsElement[i] / tumcElement[i], limt8))));
				l[i] = (mdrt8 - jrrt80 * (Dfck * Dfck) / (ksct8 + (Dfck * Dfck)) + qirt8 * tumcElement[i]
					+ lreg * TregElement[i] - lm2CD8 * 0.01) * 24d * 3600d;
			}
			CD8U = u;
			CD8L = l;
		}
		private static void ctCoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			if (NKcellsElement == null)
			{
				NKcellsElement = new double[oxModel.Item2.elementDomains[0].Count];
				for (int i = 0; i < NKcellsElement.Length; i++)
				{
					NKcellsElement[i] = 0d;
				}
			}

			for (int i = 0; i < l.Count; i++)
			{
				var cnknew = (1d + 20d * coxElement[i] / 0.2) * cnk;
				l[i] = (cnknew * NKcellsElement[i] + ptc + pti + lm1tum * 0.01) * 24d * 3600d;
			}
			CancerTransportU = u;
			CancerTransportL = l;
		}
		private static void NKcellsCoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			if (NKcellsElement == null)
			{
				NKcellsElement = new double[NKcellsModel.Item1.Elements.Count];
				for (int i = 0; i < NKcellsElement.Length; i++)
				{
					NKcellsElement[i] = 0d;
				}
			}

			if (TregElement == null)
			{
				TregElement = new double[l.Count];
				for (int i = 0; i < l.Count; i++)
				{
					TregElement[i] = 0.01;
				}
			}

			for (int i = 0; i < l.Count; i++)
			{
				var fdrnk = fdrnk0 + (1 - coxElement[i] / 0.2) * 1.025 * fdrnk0;
				l[i] = -(grrnk0 * (Math.Pow(tumcElement[i], 2) - fdrnk) / (hscnk + Math.Pow(tumcElement[i], 2))
					- pirnk * tumcElement[i] - lreg * TregElement[i] + lm2Nk * 0.01) * 24d * 3600d;
			}
			NKcellsU = u;
			NKcellsL = l;
		}
		private static void UpdateModels(Dictionary<int, IVector>[] prevStepSolutions, IStructuralModel[] modelsToReplace, ISolver[] solversToReplace,
			IConvectionDiffusionIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
		{

			//
			modelsToReplace[0] = CreateEndocModel().Item1;
			modelsToReplace[1] = CreateOxygenTransportModel(new double[] { Dox[0] * 24d * 3600d, Dox[1] * 24d * 3600d }).Item1;
			modelsToReplace[2] = CreateCancerTransportModel(Dcell[0]).Item1;
			modelsToReplace[3] = CreateGrowthModel().Item1;
			modelsToReplace[4] = CreatePressureModel(khy).Item1;
			modelsToReplace[5] = CreateCsdfModel().Item1;
			modelsToReplace[6] = CreateCvegfModel().Item1;
			modelsToReplace[7] = CreateAng1Model().Item1;
			modelsToReplace[8] = CreateAng2Model().Item1;
			modelsToReplace[9] = CreatePhisModel().Item1;
			modelsToReplace[10] = CreateNKcellsModel().Item1;
			modelsToReplace[11] = CreateCD8Model().Item1;
			modelsToReplace[12] = CreateTregModel().Item1;

			for (int i = 0; i < modelsToReplace.Length; i++)
			{
				if (i == 1)
				{
					var asymBuilder = new DenseMatrixSolver.Builder();
					asymBuilder.IsMatrixPositiveDefinite = false;
					solversToReplace[i] = asymBuilder.BuildSolver(modelsToReplace[i]);
				}
				else
					solversToReplace[i] = builder.BuildSolver(modelsToReplace[i]);
				providersToReplace[i] = new ProblemConvectionDiffusion2((Model)modelsToReplace[i], solversToReplace[i]);
				childAnalyzersToReplace[i] = new LinearAnalyzer(modelsToReplace[i], solversToReplace[i], providersToReplace[i]);
			}

		}
		private static void UpdateSolutions(ISolver[] solversToReplace)
		{
			endocNode = solversToReplace[0].LinearSystems[0].Solution.CopyToArray();
			coxNode = solversToReplace[1].LinearSystems[0].Solution.CopyToArray();
			tumcNode = solversToReplace[2].LinearSystems[0].Solution.CopyToArray();
			lgNode = solversToReplace[3].LinearSystems[0].Solution.CopyToArray();
			pSolution = solversToReplace[4].LinearSystems[0].Solution.CopyToArray();
			if (pNode == null) pNode = new double[prModel.Item1.Nodes.Count];
			int freeDofNo = 0;
			for (int i = 0; i < prModel.Item1.Nodes.Count; i++)
			{
				if (prModel.Item1.GlobalDofOrdering.GlobalFreeDofs.Contains(prModel.Item1.Nodes[i], ThermalDof.Temperature))
				{
					pNode[i] = pSolution[freeDofNo];
					freeDofNo++;
				}
			}
			CsdfNode = solversToReplace[5].LinearSystems[0].Solution.CopyToArray();
			CvegfNode = solversToReplace[6].LinearSystems[0].Solution.CopyToArray();
			a1Node = solversToReplace[7].LinearSystems[0].Solution.CopyToArray();
			a2Node = solversToReplace[8].LinearSystems[0].Solution.CopyToArray();
			phisNode = solversToReplace[9].LinearSystems[0].Solution.CopyToArray();
			NKcellsNode = solversToReplace[10].LinearSystems[0].Solution.CopyToArray();
			T8imcellsNode = solversToReplace[11].LinearSystems[0].Solution.CopyToArray();
			TregNode = solversToReplace[12].LinearSystems[0].Solution.CopyToArray();
			//endoc
			if (endocElement == null) endocElement = new double[endocModel.Item1.Elements.Count];
			foreach (var e in endocModel.Item1.Elements)
			{
				endocElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					endocElement[e.ID] += endocNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//c_ox
			dcoxdx = GetCDStrains(oxModel.Item1, Vector.CreateFromArray(coxNode));
			foreach (var e in oxModel.Item1.Elements)
			{
				coxElement[e.ID] = 0d;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					coxElement[e.ID] += coxNode[i] / (e.Nodes.Count);
				}
			}

			//Tumc
			if (tumcElement == null) tumcElement = new double[ctModel.Item1.Elements.Count];
			foreach (var e in ctModel.Item1.Elements)
			{
				tumcElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					tumcElement[e.ID] += tumcNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//lGrowth
			if (lgElement == null) lgElement = new double[gModel.Item1.Elements.Count];
			foreach (var e in gModel.Item1.Elements)
			{
				lgElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					lgElement[e.ID] += lgNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//pressure
			dpdx = GetCDStrains(prModel.Item1, Vector.CreateFromArray(pNode));
			if (pElement == null) pElement = new double[prModel.Item1.Elements.Count];
			//if (dpdx == null) dpdx = new Dictionary<int, double[]>();
			foreach (var e in prModel.Item1.Elements)
			{
				pElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					pElement[e.ID] += pNode[i] / (e.Nodes.Count);
				}
			}

			//Csdf
			if (CsdfElement == null) CsdfElement = new double[csdfModel.Item1.Elements.Count];
			foreach (var e in csdfModel.Item1.Elements)
			{
				CsdfElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					CsdfElement[e.ID] += CsdfNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//Cvegf
			if (CvegfElement == null) CvegfElement = new double[cvegfModel.Item1.Elements.Count];
			foreach (var e in cvegfModel.Item1.Elements)
			{
				CvegfElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					CvegfElement[e.ID] += CvegfNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//Ang1
			if (a1Element == null) a1Element = new double[a1Model.Item1.Elements.Count];
			foreach (var e in a1Model.Item1.Elements)
			{
				a1Element[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					a1Element[e.ID] += a1Node[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//Ang2
			if (a2Element == null) a2Element = new double[a2Model.Item1.Elements.Count];
			foreach (var e in a2Model.Item1.Elements)
			{
				a2Element[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					a2Element[e.ID] += a2Node[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//phis
			//dphisdx = GetCDStrains(phisModel.Item1, Vector.CreateFromArray(phisNode));
			if (phisElement == null) phisElement = new double[phisModel.Item1.Elements.Count];
			foreach (var e in phisModel.Item1.Elements)
			{
				phisElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					phisElement[e.ID] += phisNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//NKcell
			if (NKcellsElement == null) NKcellsElement = new double[NKcellsModel.Item1.Elements.Count];
			foreach (var e in NKcellsModel.Item1.Elements)
			{
				NKcellsElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					NKcellsElement[e.ID] += NKcellsNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//T8imcell
			if (T8imcellsElement == null) T8imcellsElement = new double[CD8Model.Item1.Elements.Count];
			foreach (var e in CD8Model.Item1.Elements)
			{
				T8imcellsElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					T8imcellsElement[e.ID] += T8imcellsNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

			//Tregs
			if (TregElement == null) TregElement = new double[TregModel.Item1.Elements.Count];
			foreach (var e in TregModel.Item1.Elements)
			{
				TregElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					TregElement[e.ID] += TregNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}

		}
		private static void ReplaceLambdaGInModel(IStructuralModel model, double[] lg)
		{
			foreach (var e in model.Elements)
			{
				var et = (ContinuumElement3DNonLinearDefGrad)e.ElementType;
				var bindFlags = BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.Static;
				FieldInfo field = typeof(ContinuumElement3DNonLinearDefGrad).GetField("lambdag", bindFlags);
				field.SetValue(et, lg[e.ID]);
			}
		}
		private static void UpdateStructuralModel(IStructuralModel[] modelsToReplace,
			ISolver[] solversToReplace, IStaticProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
		{
			if (lgElement == null)
			{
				lgElement = new double[modelsToReplace[0].Elements.Count];
				foreach (Element e in structModel.Item2.elementDomains[1])
				{
					lgElement[e.ID] = 1d;
				}
			}
			ReplaceLambdaGInModel(modelsToReplace[0], lgElement);
			UpdateStructuralLoads();
			modelsToReplace[0] = structModel.Item1;
			solversToReplace[0] = structuralBuilder.BuildSolver(modelsToReplace[0]);
			providersToReplace[0] = new ProblemStructural(modelsToReplace[0], solversToReplace[0]);
			var increments = 2;
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(modelsToReplace[0], solversToReplace[0], (INonLinearProvider)providersToReplace[0], increments);
			childAnalyzerBuilder.ResidualTolerance = 1E-5;
			childAnalyzerBuilder.MaxIterationsPerIncrement = 50;
			childAnalyzerBuilder.NumIterationsForMatrixRebuild = 5;
			childAnalyzersToReplace[0] = childAnalyzerBuilder.Build();
		}
		private static void UpdateStructuralModelSolution(IChildAnalyzer[] childAnalyzersToReplace)
		{
			if (prevUNode == null)
			{
				prevUNode = new double[structModel.Item1.Nodes.Count][];
				for (int i = 0; i < structModel.Item1.Nodes.Count; i++)
				{
					prevUNode[i] = new double[3];
				}
			}
			if (uNode == null)
			{
				uNode = new double[structModel.Item1.Nodes.Count][];
				vNode = new double[structModel.Item1.Nodes.Count][];
				for (int i = 0; i < structModel.Item1.Nodes.Count; i++)
				{
					uNode[i] = new double[3];
					vNode[i] = new double[3];
				}
			}
			if (Displacements == null)
			{
				Displacements = new Dictionary<int, IVector>();
				Displacements.Add(0, Vector.CreateZero(structModel.Item1.GlobalDofOrdering.NumGlobalFreeDofs));
			}
			Displacements = childAnalyzersToReplace[0].Responses;

			int freeDofNo = 0;
			for (int i = 0; i < structModel.Item1.Nodes.Count; i++)
			{
				if (structModel.Item1.GlobalDofOrdering.GlobalFreeDofs.Contains(structModel.Item1.Nodes[i], StructuralDof.TranslationX))
				{
					uNode[i][0] = Displacements[0][freeDofNo];
					vNode[i][0] = (uNode[i][0] - prevUNode[i][0]) / timestep;
					freeDofNo++;
				}
				if (structModel.Item1.GlobalDofOrdering.GlobalFreeDofs.Contains(structModel.Item1.Nodes[i], StructuralDof.TranslationY))
				{
					uNode[i][1] = Displacements[0][freeDofNo];
					vNode[i][1] = (uNode[i][1] - prevUNode[i][1]) / timestep;
					freeDofNo++;
				}
				if (structModel.Item1.GlobalDofOrdering.GlobalFreeDofs.Contains(structModel.Item1.Nodes[i], StructuralDof.TranslationZ))
				{
					uNode[i][2] = Displacements[0][freeDofNo];
					vNode[i][2] = (uNode[i][1] - prevUNode[i][1]) / timestep;
					freeDofNo++;
				}
			}

			foreach (Element e in structModel.Item1.Elements)
			{
				vElement[e.ID] = new double[3];
				uElement[e.ID] = new double[3];
				foreach (Node n in e.Nodes)
				{
					for (int i = 0; i < 3; i++)
					{
						vElement[e.ID][i] += vNode[n.ID][i] / e.Nodes.Count;
						uElement[e.ID][i] += uNode[n.ID][i] / e.Nodes.Count;
					}
				}
			}

			Strains = GetStrainsStresses(structModel.Item1.Elements.Count).Item1;
			Stresses = GetStrainsStresses(structModel.Item1.Elements.Count).Item2;
			uXt = StructuralSpaceTimeDerivatives(Strains, PreviousStrains);
		}
		private static Tuple<Model, IModelReader> CreateEndocModel()
		{
			ComsolMeshReader3 modelReader;
			Model model;
			if (endocModel == null)
			{
				Console.WriteLine("Creating endoc Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader3(filename, new double[] { 1 }, new double[] { 1 }, conv0, new double[] { 0 });
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating endoc Model...");
				modelReader = (ComsolMeshReader3)endocModel.Item2;
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 1 }, conv0, new double[] { 0 });
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (endocElement == null)
			{
				endocElement = new double[model.Elements.Count];
				foreach (var e in model.Elements)
				{
					endocElement[e.ID] = 0.5;
				}
			}

			if (CvegfElement == null) CvegfElement = new double[model.Elements.Count];

			int[] domainIDs = new int[] { 0 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					double fendoc = 0.021429 + kec * canti * endocElement[element.ID] * 24d * 3600d;
					var endocMaterial = new ConvectionDiffusionMaterial(1, 1, conv0[0], 0);
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(endocMaterial, fendoc, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateOxygenTransportModel(double[] k)
		{
			ComsolMeshReader5 modelReader;
			Model model;

			if (oxModel == null)
			{
				Console.WriteLine("Creating Oxygen Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0, 1 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
				modelReader = new ComsolMeshReader5(filename, new double[] { 1, 1 }, OxygenTransportCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Oxygen Model...");
				modelReader = (ComsolMeshReader5)oxModel.Item2;
				OxygenTransportCoefficientsCalculation(modelReader.elementIDsPerDomain, OxygenTransportK, OxygenTransportU, OxygenTransportL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1, 1 }, OxygenTransportK, OxygenTransportU, OxygenTransportL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}
			if (coxElement == null)
			{
				coxElement = new double[model.Elements.Count];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					coxElement[i] = 0d;/* 0.9673;*/
				}
			}
			if (tumcElement == null)
			{
				tumcElement = new double[modelReader.elementDomains[0].Count];
				for (int i = 0; i < modelReader.elementDomains[0].Count; i++)
				{
					tumcElement[i] = 0.96;/* 0.9673;*/
				}
			}
			if (endocElement == null)
			{
				endocElement = new double[modelReader.elementDomains[0].Count];
				for (int i = 0; i < modelReader.elementDomains[0].Count; i++)
				{
					endocElement[i] = 0.5;/* 0.9673;*/
				}
			}

			int[] domainIDs = new int[] { 0, 1 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					double Rox;
					var material = new ConvectionDiffusionMaterial
						(1, OxygenTransportK[element.ID], OxygenTransportU[element.ID], OxygenTransportL[element.ID]);
					if (domainID == 0)
					{
						//Rox = (Dox[domainID] / Lwv * Svin[domainID] * dd0[element.ID] * endocElement[element.ID] * (cvox - coxElement[element.ID]) 
						//-((Aox[0] * coxElement[element.ID]) / (kox[0] + coxElement[element.ID])) * tumcElement[element.ID]) * (24d * 3600d);
						Rox = (Dox[domainID] / Lwv * Svin[domainID] * dd0[element.ID] * endocElement[element.ID] * cvox) * (24d * 3600d);
					}
					else
					{
						//Rox = (Dox[domainID] / Lwv * Svin[domainID] * (cvox - coxElement[element.ID])) * (24d * 3600d);
						Rox = (Dox[domainID] / Lwv * Svin[domainID] * cvox) * (24d * 3600d);
					}
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(material, Rox, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateGrowthModel()
		{
			ComsolMeshReader4 modelReader;
			Model model;

			if (gModel == null)
			{
				Console.WriteLine("Creating Growth Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 1 }, new double[] { 0 }, lgCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Growth Model...");
				modelReader = (ComsolMeshReader4)gModel.Item2;
				lgCoefficientsCalculation(lgU, lgL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 0 }, lgU, lgL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (lgElement == null)
			{
				lgElement = new double[oxModel.Item1.Elements.Count];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					lgElement[i] = 1;
				}
			}

			//var materialODE = new ConvectionDiffusionMaterial(1, 0, conv0[0], 0);
			////double[] Grox = new double[model.Elements.Count];
			//int[] domainIDs = new int[] { 0, };
			//foreach (int domainID in domainIDs)
			//{
			//	foreach (Element element in modelReader.elementDomains[domainID])
			//	{
			//		var Grox = (loxc[domainID] * coxElement[element.ID]) / (coxElement[element.ID] + Koxc[domainID]);
			//		var Rtumc = Grox * Sfn - (ptc + pti + lm1tum * 0.01) * tumcElement[element.ID];
			//		var fg = 24d * 3600d * Rtumc * lgElement[element.ID] / 3d;
			//		var nodes = (IReadOnlyList<Node>)element.Nodes;
			//		var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, fg, ThermalDof.Temperature);
			//		var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
			//		var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
			//		model.BodyLoads.Add(bodyLoadElement);
			//	}
			//}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateCancerTransportModel(double k)
		{
			ComsolMeshReader4 modelReader;
			Model model;
			if (ctModel == null)
			{
				Console.WriteLine("Creating Cancer Transport Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 1 }, new double[] { 0 }, ctCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Cancer Transport Model...");
				modelReader = (ComsolMeshReader4)ctModel.Item2;
				ctCoefficientsCalculation(CancerTransportU, CancerTransportL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 0 }, CancerTransportU, CancerTransportL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (tumcElement == null)
			{
				tumcElement = new double[oxModel.Item1.Elements.Count];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					tumcElement[i] = 0.96;
				}
			}

			foreach (Element element in modelReader.elementDomains[0])
			{
				var materialODE = new ConvectionDiffusionMaterial(1, 0, conv0[0], CancerTransportL[element.ID]);
				var Grox = (loxc[0] * coxElement[element.ID]) / (coxElement[element.ID] + Koxc[0]);
				var Rtumc = (Grox * Sfn) * 24d * 3600d;
				var nodes = (IReadOnlyList<Node>)element.Nodes;
				var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, Rtumc, ThermalDof.Temperature);
				var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
				var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
				model.BodyLoads.Add(bodyLoadElement);
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateCD8Model()
		{
			ComsolMeshReader4 modelReader;
			Model model;
			if (CD8Model == null)
			{
				Console.WriteLine("Creating T8 immune cells Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 1 }, new double[] { 0 }, CD8CoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating T8 immune cells Model...");
				modelReader = (ComsolMeshReader4)CD8Model.Item2;
				CD8CoefficientsCalculation(CD8U, CD8L);
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 0 }, CD8U, CD8L);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (T8imcellsElement == null)
			{
				T8imcellsElement = new double[model.Elements.Count];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					T8imcellsElement[i] = 0d;
				}
			}

			foreach (Element element in modelReader.elementDomains[0])
			{
				var materialODE = new ConvectionDiffusionMaterial(1, 0, CD8U[element.ID], CD8L[element.ID]);
				var dt8new = dt8 + 8 * coxElement[element.ID] / 0.2 * dt8;
				var Dfck = (dt8new * Math.Pow(T8imcellsElement[element.ID] / tumcElement[element.ID], limt8) * tumcElement[element.ID]
					/ (st8 + (Math.Pow(T8imcellsElement[element.ID] / tumcElement[element.ID], limt8))));
				var Rt8im = (scd8 + rst8 * NKcellsElement[element.ID] * tumcElement[element.ID]) * 24d * 3600d;
				var nodes = (IReadOnlyList<Node>)element.Nodes;
				var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, Rt8im, ThermalDof.Temperature);
				var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
				var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
				model.BodyLoads.Add(bodyLoadElement);
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateTregModel()
		{
			ComsolMeshReader3 modelReader;
			Model model;
			if (TregModel == null)
			{
				Console.WriteLine("Creating Treg cells Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader3(filename, new double[] { 1 }, new double[] { 0 }, conv0, new double[] { TregL });
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Treg cells Model...");
				modelReader = (ComsolMeshReader3)TregModel.Item2;
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 0 }, conv0, new double[] { TregL });
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (TregElement == null)
			{
				TregElement = new double[model.Elements.Count];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					TregElement[i] = 0.01;
				}
			}

			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreatePressureModel(double[] k)
		{
			ComsolMeshReader4 modelReader;
			Model model;
			if (prModel == null)
			{
				Console.WriteLine("Creating Pressure Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0, 1 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 0, 0 },
					k, PressureCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);

			}
			else
			{
				Console.WriteLine("Updating Pressure Model...");
				modelReader = (ComsolMeshReader4)prModel.Item2;
				PressureCoefficientsCalculation(PressureU, PressureL);
				modelReader = modelReader.UpdateModelReader(new double[] { 0, 0 },
					k, PressureU, PressureL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			//Dirichlet BC
			int[] boundaryIDs = new int[] { 6, 8, 9 };
			foreach (int boundaryID in boundaryIDs)
			{
				foreach (Node node in modelReader.nodeBoundaries[boundaryID])
				{
					model.NodesDictionary[node.ID].Constraints.Add(new Constraint()
					{
						Amount = 0,
						DOF = ThermalDof.Temperature
					});
				}
			}
			foreach (Node n in model.Nodes)
			{
				if (n.Constraints.Count == 0)
					pressureModelFreeDOFs++;
			}
			//

			if (pElement == null)
			{
				pElement = new double[model.Elements.Count];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					pElement[i] = 0;/* 0.9673;*/
				}
			}
			if (PreviousStrains == null)
			{
				PreviousStrains = new double[model.Elements.Count][];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					PreviousStrains[i] = new double[6];
				}
			}
			if (Strains == null)
			{
				Strains = new double[model.Elements.Count][];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					Strains[i] = new double[6];
				}
			}
			if (uXt == null)
			{
				uXt = new Dictionary<int, double[]>();
				for (int i = 0; i < model.Elements.Count; i++)
				{
					uXt[i] = new double[3];
				}
			}

			if (PressureL == null)
			{
				for (int i = 0; i < model.Elements.Count; i++)
				{
					PressureL.Add(i, 0d);
				}
			}

			int[] domainIDs = new int[] { 0, 1 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					var prMaterial = new ConvectionDiffusionMaterial(0, k[domainID], conv0[0], PressureL[element.ID]);
					var Grox = (loxc[domainID] * coxElement[element.ID]) / (coxElement[element.ID] + Koxc[domainID]);
					double Rtumc;
					double fp;//= Rtumc + pv * PressureL[element.ID] - uXt[element.ID].Sum();
					var cnknew = (1d + 20d * coxElement[element.ID] / 0.2) * cnk;
					if (domainID == 0)
					{
						Rtumc = (Grox * Sfn - (cnknew * NKcellsElement[element.ID] + ptc + pti + lm1tum * 0.01) * tumcElement[element.ID]) * 24d * 3600d;
						fp = Rtumc + pv * PressureL[element.ID] - uXt[element.ID].Sum();
					}
					else
						fp = (lp[1] * Svin[1] * pv + lpSv * Math.PI) * 24d * 3600d - uXt[element.ID].Sum();
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(prMaterial, fp, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateCsdfModel()
		{
			ComsolMeshReader3 modelReader;
			Model model;

			if (csdfModel == null)
			{
				Console.WriteLine("Creating Csdf Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader3(filename, new double[] { 1 }, new double[] { Dvegf[0] }, conv0, new double[] { l13 * 24 * 3600 });
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Csdf Model...");
				modelReader = (ComsolMeshReader3)csdfModel.Item2;
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { Dvegf[0] }, conv0, new double[] { l13 * 24d * 3600d });
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (CvegfElement == null) CvegfElement = new double[model.Elements.Count];

			var materialODE = new ConvectionDiffusionMaterial(1, Dvegf[0], conv0[0], l13 * 24 * 3600);
			int[] domainIDs = new int[] { 0, };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					double Ga = 2 - coxElement[element.ID] / 0.2;
					var Cs_stD = 24 * 3600 * ((l10 * Ga * tumcElement[element.ID] * T0 * 100 / Cs0) / 2 +
						(l10 * endocElement[element.ID] * 100 * Cv0 * CvegfElement[element.ID] / Cs0) / 2);
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, Cs_stD, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateCvegfModel()
		{
			ComsolMeshReader5 modelReader;
			Model model;

			if (cvegfModel == null)
			{
				Console.WriteLine("Creating Cvegf Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader5(filename, new double[] { 1 }, CvegfCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Cvegf Model...");
				modelReader = (ComsolMeshReader5)cvegfModel.Item2;
				CvegfCoefficientsCalculation(modelReader.elementIDsPerDomain, CvegfK, CvegfU, CvegfL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, CvegfK, CvegfU, CvegfL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			int[] domainIDs = new int[] { 0, };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					var materialODE = new ConvectionDiffusionMaterial(1, CvegfK[element.ID], CvegfU[element.ID], CvegfL[element.ID]);
					double Ga = 2 - coxElement[element.ID] / 0.2;
					var Rvegf = 24 * 3600 * ((l10 * Ga * tumcElement[element.ID] * T0 * T0in * coxElement[element.ID] / 0.2 / Cvegf0)
						 + kantivegf * canti);
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, Rvegf, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateAng1Model()
		{
			ComsolMeshReader3 modelReader;
			Model model;

			if (a1Model == null)
			{
				Console.WriteLine("Creating Ang1 Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader3(filename, new double[] { 1 }, new double[] { 0 },
					conv0, new double[] { m1 / 2d * 24d * 3600d });
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Ang1 Model...");
				modelReader = (ComsolMeshReader3)a1Model.Item2;
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 0 },
					conv0, new double[] { m1 / 2d * 24d * 3600d });
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (a1Element == null) a1Element = new double[model.Elements.Count];

			var materialODE = new ConvectionDiffusionMaterial(1, 0, conv0[0], m1 / 2d * 24d * 3600d);
			//double[] Grox = new double[model.Elements.Count];
			int[] domainIDs = new int[] { 0 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					double Ga = 2 - coxElement[element.ID] / 0.2;
					double Sv;
					Sv = Svin[0] * dd0[element.ID] * endocElement[element.ID];
					var Ra1 = 24d * 3600d * ((b1 * 1e-10 * Ga * pc[0] * 0.5 * Sv / (pc0 * 200) / a10) + m1) / 2;
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, Ra1, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateAng2Model()
		{
			ComsolMeshReader3 modelReader;
			Model model;

			if (a2Model == null)
			{
				Console.WriteLine("Creating Ang2 Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader3(filename, new double[] { 1 }, new double[] { 0 },
					conv0, new double[] { m2 * 24d * 3600d });
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Ang2 Model...");
				modelReader = (ComsolMeshReader3)a2Model.Item2;
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 0 },
					conv0, new double[] { m2 * 24d * 3600d });
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (a2Element == null) a2Element = new double[model.Elements.Count];

			var materialODE = new ConvectionDiffusionMaterial(1, 0, conv0[0], m2 * 24d * 3600d);
			//double[] Grox = new double[model.Elements.Count];
			int[] domainIDs = new int[] { 0 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					double Ga = 2 - coxElement[element.ID] / 0.2;
					double Sv;
					Sv = Svin[0] * dd0[element.ID] * endocElement[element.ID];
					var Ra2 = 24d * 3600d * (b2 * 1e-11 * endocElement[element.ID] * Sv * Ga / a20) / 2;
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, Ra2, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreatePhisModel()
		{
			ComsolMeshReader4 modelReader;
			Model model;

			if (phisModel == null)
			{
				Console.WriteLine("Creating phis Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 1 }, new double[] { 0 }, phisCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating phis Model...");
				modelReader = (ComsolMeshReader4)phisModel.Item2;
				phisCoefficientsCalculation(phisU, phisL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 0 }, phisU, phisL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (phisElement == null)
			{
				phisElement = new double[endocModel.Item1.Elements.Count];
				for (int i = 0; i < phisElement.Length; i++)
				{
					phisElement[i] = 0.3;
				}
			}
			if (dphisdx == null)
			{
				dphisdx = new double[phisElement.Length][];
				for (int i = 0; i < phisElement.Length; i++)
				{
					dphisdx[i] = new double[3];
				}
			}

			int[] domainIDs = new int[] { 0, };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					var Grox = (loxc[domainID] * coxElement[element.ID]) / (coxElement[element.ID] + Koxc[domainID]);
					var cnknew = (1d + 20d * coxElement[element.ID] / 0.2) * cnk;
					var Rtumc = Grox * Sfn - (cnknew * NKcellsElement[element.ID] + ptc + pti + lm1tum * 0.01) * tumcElement[element.ID];
					var materialODE = new ConvectionDiffusionMaterial(1, 0, phisU[element.ID], phisL[element.ID]);
					var fphis = 24d * 3600d * (Rtumc - (vElement[element.ID][0] * dphisdx[element.ID][0] +
						vElement[element.ID][1] * dphisdx[element.ID][1] + vElement[element.ID][2] * dphisdx[element.ID][2]));
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, fphis, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateNKcellsModel()
		{
			ComsolMeshReader4 modelReader;
			Model model;

			if (NKcellsModel == null)
			{
				Console.WriteLine("Creating NKcells Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 1 }, new double[] { 0 }, NKcellsCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating NKcells Model...");
				modelReader = (ComsolMeshReader4)NKcellsModel.Item2;
				NKcellsCoefficientsCalculation(NKcellsU, NKcellsL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 0 }, NKcellsU, NKcellsL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			if (NKcellsElement == null)
			{
				NKcellsElement = new double[model.Elements.Count];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					NKcellsElement[i] = 0d;
				}
			}

			foreach (Element element in modelReader.elementDomains[0])
			{
				var materialODE = new ConvectionDiffusionMaterial(1, 0, conv0[0], NKcellsL[element.ID]);
				var Rnk = snkcon * 24d * 3600d;
				var nodes = (IReadOnlyList<Node>)element.Nodes;
				var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, Rnk, ThermalDof.Temperature);
				var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
				var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
				model.BodyLoads.Add(bodyLoadElement);
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateStructuralModel(double[] MuLame, double[] PoissonV, IDynamicMaterial[] commonDynamicMaterialProperties,
			double b, double[] l, double[] lambdag)
		{
			double[] C1 = new double[MuLame.Length];
			double[] C2 = new double[MuLame.Length];
			double[] bulkModulus = new double[MuLame.Length];
			for (int i = 0; i < MuLame.Length; i++)
			{
				//poissonV[i] = 0.2;
				C1[i] = MuLame[i] / 2;
				C2[i] = 0;
				bulkModulus[i] = 2 * MuLame[i] * (1 + PoissonV[i]) / (3 * (1 - 2 * PoissonV[i]));
			}
			string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", inputFile);
			ComsolMeshReader1 modelReader;
			if (lambdag == null)
			{
				modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties);
			}
			else
			{
				modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties, lambdag);
			}
			Model model = modelReader.CreateModelFromFile();
			//Boundary Conditions
			var lx = l[0];
			var ly = l[1];
			var lz = l[2];
			var distributedLoad = new DistributedLoad(lx, ly, lz);


			int[] boundaryIDs = new int[] { 0, 3 };
			foreach (int boundaryID in boundaryIDs)
			{
				foreach (Node node in modelReader.nodeBoundaries[boundaryID])
				{
					node.Constraints.Add(new Constraint()
					{
						Amount = b,
						DOF = StructuralDof.TranslationX
					});
				}
			}
			boundaryIDs = new int[] { 1, 4 };
			foreach (int boundaryID in boundaryIDs)
			{
				foreach (Node node in modelReader.nodeBoundaries[boundaryID])
				{
					node.Constraints.Add(new Constraint()
					{
						Amount = b,
						DOF = StructuralDof.TranslationY
					});
				}
			}
			boundaryIDs = new int[] { 2, 7 };
			foreach (int boundaryID in boundaryIDs)
			{
				foreach (Node node in modelReader.nodeBoundaries[boundaryID])
				{
					node.Constraints.Add(new Constraint()
					{
						Amount = b,
						DOF = StructuralDof.TranslationZ
					});
				}
			}
			if (dpdx == null)
			{
				dpdx = new double[model.Elements.Count][];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					dpdx[i] = new double[3];
				}
			}
			int[] domainIDs = new int[] { 0, 1 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var bodyLoadX = new GravityLoad(1d, -dpdx[element.ID][0], StructuralDof.TranslationX);
					var bodyLoadElementFactoryX = new BodyLoadElementFactory(bodyLoadX, model);
					var bodyLoadElementX = bodyLoadElementFactoryX.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElementX);
					var bodyLoadY = new GravityLoad(1d, -dpdx[element.ID][1], StructuralDof.TranslationY);
					var bodyLoadElementFactoryY = new BodyLoadElementFactory(bodyLoadY, model);
					var bodyLoadElementY = bodyLoadElementFactoryY.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElementY);
					var bodyLoadZ = new GravityLoad(1d, -dpdx[element.ID][2], StructuralDof.TranslationZ);
					var bodyLoadElementFactoryZ = new BodyLoadElementFactory(bodyLoadZ, model);
					var bodyLoadElementZ = bodyLoadElementFactoryZ.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElementZ);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static void UpdateStructuralLoads()
		{
			structModel.Item1.BodyLoads.Clear();
			int[] domainIDs = new int[] { 0, 1 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in structModel.Item2.elementDomains[domainID])
				{
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var bodyLoadX = new GravityLoad(1d, -dpdx[element.ID][0], StructuralDof.TranslationX);
					var bodyLoadElementFactoryX = new BodyLoadElementFactory(bodyLoadX, structModel.Item1);
					var bodyLoadElementX = bodyLoadElementFactoryX.CreateElement(CellType.Tet4, nodes);
					structModel.Item1.BodyLoads.Add(bodyLoadElementX);
					var bodyLoadY = new GravityLoad(1d, -dpdx[element.ID][1], StructuralDof.TranslationY);
					var bodyLoadElementFactoryY = new BodyLoadElementFactory(bodyLoadY, structModel.Item1);
					var bodyLoadElementY = bodyLoadElementFactoryY.CreateElement(CellType.Tet4, nodes);
					structModel.Item1.BodyLoads.Add(bodyLoadElementY);
					var bodyLoadZ = new GravityLoad(1d, -dpdx[element.ID][2], StructuralDof.TranslationZ);
					var bodyLoadElementFactoryZ = new BodyLoadElementFactory(bodyLoadZ, structModel.Item1);
					var bodyLoadElementZ = bodyLoadElementFactoryZ.CreateElement(CellType.Tet4, nodes);
					structModel.Item1.BodyLoads.Add(bodyLoadElementZ);
				}
			}
		}
		private static IVectorView[] SolveModelsWithStructuralAnalyzer(Model[] models, IModelReader[] modelReaders)
		{
			double[] E = new double[] { 7e4, 2.1e4 };
			double[] poissonV = new double[] { .45, .2 };
			double[] muLame = new double[] { E[0] / 2d / (1 + poissonV[0]), E[1] / 2d / (1 + poissonV[1]) };
			IDynamicMaterial[] dynamicMaterials = new DynamicMaterial[] { new DynamicMaterial(0, 0, 0, true), new DynamicMaterial(0, 0, 0, true) };
			structModel = CreateStructuralModel(muLame, poissonV, dynamicMaterials, 0, new double[] { 0, 0, 0 }, lgElement);//.Item1; // new Model();
			var structuralModel = structModel.Item1;
			var structuralSolver = structuralBuilder.BuildSolver(structuralModel);
			var structuralProvider = new ProblemStructural(structuralModel, structuralSolver);
			//var structuralChildAnalyzer = new LinearAnalyzer(structuralModel, structuralSolver, structuralProvider);
			var increments = NewtonRaphsonIncrements;
			var structuralChildAnalyzerBuilder = new LoadControlAnalyzer.Builder(structuralModel, structuralSolver, structuralProvider, increments);
			structuralChildAnalyzerBuilder.ResidualTolerance = NewtonRaphsonTolerarance;
			structuralChildAnalyzerBuilder.MaxIterationsPerIncrement = NewtonRaphosnIterations;
			structuralChildAnalyzerBuilder.NumIterationsForMatrixRebuild = NewtonRaphsonIterForMatrixRebuild;
			//childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
			LoadControlAnalyzer structuralChildAnalyzer = structuralChildAnalyzerBuilder.Build();
			var structuralParentAnalyzer = new StaticAnalyzer(UpdateStructuralModel, UpdateStructuralModelSolution, structuralModel, structuralSolver,
				structuralProvider, structuralChildAnalyzer);
			structuralParentAnalyzer.Initialize();

			Vector[] initialValues = new Vector[models.Length];
			var value0 = new Dictionary<int, double[]>();
			for (int i = 0; i < models.Length; i++)
			{
				double[] v0;
				if (i == 4)
					v0 = new double[pressureModelFreeDOFs];
				else
					v0 = new double[models[i].Nodes.Count];
				value0.Add(i, v0);
			}
			foreach (Node node in models[0].Nodes)
			{
				value0[0][node.ID] = 0.5;
			}
			foreach (Node node in models[1].Nodes)
			{
				value0[1][node.ID] = 0d; /* 0.96733;*/
			}
			foreach (Node node in models[2].Nodes)
			{
				value0[2][node.ID] = 0.96;
			}
			foreach (Node node in models[3].Nodes)
			{
				value0[3][node.ID] = 1;
			}
			for (int i = 0; i < pressureModelFreeDOFs; i++)
			{
				value0[4][i] = 0; /* 0.96733;*/
			}
			foreach (Node node in models[5].Nodes)
			{
				value0[5][node.ID] = 0;
			}
			foreach (Node node in models[6].Nodes)
			{
				value0[6][node.ID] = 0;
			}
			foreach (Node node in models[7].Nodes)
			{
				value0[7][node.ID] = 0;
			}
			foreach (Node node in models[8].Nodes)
			{
				value0[8][node.ID] = 0;
			}
			foreach (Node node in models[9].Nodes)
			{
				value0[9][node.ID] = 0.3;
			}
			foreach (Node node in models[10].Nodes)
			{
				value0[10][node.ID] = 0d;
			}
			foreach (Node node in models[11].Nodes)
			{
				value0[11][node.ID] = 0d;
			}
			foreach (Node node in models[12].Nodes)
			{
				value0[12][node.ID] = 0.01;
			}

			IConvectionDiffusionIntegrationProvider[] providers = new IConvectionDiffusionIntegrationProvider[models.Length];
			IChildAnalyzer[] childAnalyzers = new IChildAnalyzer[models.Length];
			var solvers = new ISolver[models.Length];
			for (int i = 0; i < models.Length; i++)
			{
				initialValues[i] = Vector.CreateFromArray(value0[i]);
				//var builder = new DenseMatrixSolver.Builder();
				//builder.IsMatrixPositiveDefinite = false;
				if (i == 0)
				{
					var asymBuilder = new DenseMatrixSolver.Builder();
					asymBuilder.IsMatrixPositiveDefinite = false;
					solvers[i] = asymBuilder.BuildSolver(models[i]);
				}
				else
					solvers[i] = builder.BuildSolver(models[i]);
				providers[i] = new ProblemConvectionDiffusion2(models[i], solvers[i]);
				childAnalyzers[i] = new LinearAnalyzer(models[i], solvers[i], providers[i]);
			}

			var parentAnalyzer = new ConvectionDiffusionImplicitDynamicAnalyzerMultiModel2(UpdateModels, UpdateSolutions, models, solvers,
				providers, childAnalyzers, timestep, time, structuralParentAnalyzer: structuralParentAnalyzer, initialTemperature: initialValues);
			//var parentAnalyzer = new ConvectionDiffusionImplicitDynamicAnalyzerMultiModel(UpdateModels, models, solvers,
			//	providers, childAnalyzers, timestep, time, initialTemperature: initialValues);
			parentAnalyzer.Initialize();

			for (int i = 0; i < time / timestep; i++)
			{
				parentAnalyzer.SolveTimestep(i);
				prevUNode = uNode.Clone() as double[][];
				PreviousStrains = Strains.Clone() as double[][];
				//structuralParentAnalyzer.SolveTimestep(i);
				Paraview(i);
			}

			return solvers.Select(x => x.LinearSystems[subdomainID].Solution).ToArray();
		}
	}
}