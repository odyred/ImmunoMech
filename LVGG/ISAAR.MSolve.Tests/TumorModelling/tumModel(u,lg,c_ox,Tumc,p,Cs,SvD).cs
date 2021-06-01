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

namespace ISAAR.MSolve.Tests
{
	public class tumModel_u_lg_c_ox_Tumc_p_Cs_SvD_
	{
		private const double timestep = 1;
		private const double time = 30;
		private const int subdomainID = 0;
		private static readonly double[] loxc = new double[] { .07 / 24 / 3600, 1.0 / 24 / 3600 }; //1/s
		private static readonly double[] Aox = new double[] { 2200.0 / 24 / 3600, 2200.0 / 24 / 3600 }; //mol/(m^3*s)
		private static readonly double[] Dox = new double[] { 1.78e-9, 1.79e-9 }; //m^2/s
		private static readonly double[] kox = new double[] { .00464, .00464 }; //mol/m^3
		private static readonly double[] Koxc = new double[] { 0.0083, 0.0083 }; //mol/m^3
		private static readonly double[] Dcell = new double[] { 5.4e-3, 1.8e-4 }; //m^2/s
		private static readonly double[] khy = new double[] { 6.5e-11, 6.5e-11  };  //m^3/(Pa*d);
		private static readonly double[] lp = new double[] { 2.7e-12, 2.7e-12 };  //m^2*s/kg;
		private static double T0 = 1; //kg/m^3
		private static double Cv0 = 1; //kg/m^3
		private static double Cs0 = 1; //kg/m^3
		private static double l2 = 1e-5; //m^3/kg/s
		private static double l4 = 1e-11; //m/s
		private static double l10 = 6.8e-9; //1/s
		private static double l13 = 4e-9; //1/s
		private static double cvox = 0.2; //mol/m^3
		private static double Lwv = 5e-6; //m
		private static double mtox = 8e-3 * 1.1 * 1e-6 / 3600; //m^2/s
		private static double pv = 4000; //kg/m/s^2
		private static double Svin = 7000; //1/m
		private static double WvsTc = 1;
		private static double WvsSv = 1;
		private static double xn = 1e-12;//m^5/kg/s
		private static double Dsv = 1e-7 * 1e-8;//m^2/s
		private static double s1 = 1e6;//m^3/mol
		private static double s2 = 1e6;//m^3/mol
		private static double a10 = 1;//m^3/mol
		private static double a20 = 1;//m^3/mol
		private static double aD = 1;
		private static double bD = 1;
		private static double[][] conv0 = new double[][] { new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 } };
		//private static double fox = -((Aox * c_ox) / (kox + c_ox * cvox)) * 0.3;
		private static SkylineSolver.Builder builder = new SkylineSolver.Builder();
		private static DenseMatrixSolver.Builder ctBuilder = new DenseMatrixSolver.Builder();
		private static SkylineSolver.Builder structuralBuilder = new SkylineSolver.Builder();
		private static double[] lgNode;
		private static double[] lgElement;
		private static double[] CsNode;
		private static double[] CsElement;
		private static double[] c_oxNode;
		private static double[] c_oxElement;
		private static double[] tumcNode;
		private static double[] pSolution;
		private static double[] tumcElement;
		private static double[] pNode;
		private static double[] pElement;
		private static double[] SvDNode;
		private static double[] SvDElement;
		private static Dictionary<int, double> uof = new Dictionary<int, double>();
		private static Dictionary<int, double> vof = new Dictionary<int, double>();
		private static Dictionary<int, double> wof = new Dictionary<int, double>();
		private static Dictionary<int, double> duof = new Dictionary<int, double>();
		private static Dictionary<int, double> dvof = new Dictionary<int, double>();
		private static Dictionary<int, double> dwof = new Dictionary<int, double>();
		private static double[] dcoxdvf;
		private static double[] dd0;
		private static Dictionary<int, double[]> aNode = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> aElement = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> vNode = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> vElement = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> uNode = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> uElement = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> dcoxNode = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> dcoxElement; /*= new Dictionary<int, double[]>();*/
		private static Dictionary<int, double[]> dpNode = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> dpSolution = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> dpElement; /*= new Dictionary<int, double[]>();*/
		private static Dictionary<int, double[]> ddpNode = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> ddpSolution = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> ddpElement; /*= new Dictionary<int, double[]>();*/
		private static Dictionary<int, double[]> ddcoxNode = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> ddcoxElement; /*= new Dictionary<int, double[]>();*/
		private static Dictionary<int, double[]> dCsNode = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> dCsElement; /*= new Dictionary<int, double[]>();*/
		private static Dictionary<int, double[]> ddCsNode = new Dictionary<int, double[]>();
		private static Dictionary<int, double[]> ddCsElement; /*= new Dictionary<int, double[]>();*/
		//private static Dictionary<int, double[]> ddcoxNode = new Dictionary<int, double[]>();
		//private static Dictionary<int, double[]> ddcoxElement = new Dictionary<int, double[]>();
		private static Dictionary<int, double> Grox = new Dictionary<int, double>();
		private static Dictionary<int, double> RTumc = new Dictionary<int, double>();
		private static Dictionary<int, double> fsv = new Dictionary<int, double>();
		private static Dictionary<int, double[]> OxygenTransportU;
		private static Dictionary<int, double> OxygenTransportL;
		private static Dictionary<int, double[]> CancerTransportU;
		private static Dictionary<int, double> CancerTransportL;
		private static Dictionary<int, double[]> PressureU;
		private static Dictionary<int, double> PressureL;
		private static Dictionary<int, double> SvDK;
		private static Dictionary<int, double[]> SvDU;
		private static Dictionary<int, double> SvDL;
		private static Dictionary<int, IVector> Accelerations;
		private static Dictionary<int, IVector> Velocities;
		private static double[][] PreviousSpaceDerivatives;
		private static double[][] SpaceDerivatives;
		private static Dictionary<int, double[]> uXt;
		private static Dictionary<int, IVector> Displacements;
		private static Tuple<Model, IModelReader> oxModel, gModel, ctModel, prModel, csModel, SvDModel, structModel;
		private static int pressureModelFreeDOFs = 0;

		[Fact]
		private static void RunTest()
		{
			var DoxDays = new double[Dox.Length];
			for (int i = 0; i < Dox.Length; i++)
			{
				DoxDays[i] = 24 * 3600 * Dox[i];
			}

			var DcellDays = new double[Dcell.Length];
			for (int i = 0; i < Dox.Length; i++)
			{
				DcellDays[i] = 24 * 3600 * Dcell[i];
			}

			SvDModel = CreateSvDModel();
			oxModel = CreateOxygenTransportModel(DoxDays);
			ctModel = CreateCancerTransportModel(DcellDays[0]);
			gModel = CreateGrowthModel();
			prModel = CreatePressureModel(khy);
			csModel = CreateCsModel();
			var models = new[] { SvDModel.Item1, oxModel.Item1, ctModel.Item1, gModel.Item1, prModel.Item1, csModel.Item1};
			var modelReaders = new[] { SvDModel.Item2, oxModel.Item2, ctModel.Item2, gModel.Item2, prModel.Item2, csModel.Item2};
			IVectorView[] solutions = SolveModelsWithNewmark(models, modelReaders);

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

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"oxygen\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < numberOfPoints; i++)
					outputFile.WriteLine($"{c_oxNode[i]} ");
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

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"Cs\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < CsNode.Length; i++)
					outputFile.WriteLine($"{CsNode[i]} ");
				for (int i = CsNode.Length; i < numberOfPoints; i++)
					outputFile.WriteLine($"0 ");
				outputFile.WriteLine("              </DataArray>");

				outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"SvD\" NumberOfComponents=\"1\" format=\"ascii\">");
				for (int i = 0; i < CsNode.Length; i++)
					outputFile.WriteLine($"{SvDNode[i]} ");
				for (int i = CsNode.Length; i < numberOfPoints; i++)
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
		private static double[][] GetStrains(int elementsNo)
		{
			if (structModel == null)
			{
				double[][] strains = new double[elementsNo][];
				for (int i = 0; i < elementsNo; i++)
				{
					strains[i] = new double[6];
				}
				return strains;
			}
			else
			{
				IList<Element> elements = structModel.Item1.Elements;
				double[][] strains = new double[elements.Count][];
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
				}
				return strains;
			}
		}
		private static Dictionary<int, double[]> StructuralSpaceTimeDerivatives(double[][] current, double[][] previous)
		{
			Dictionary<int, double[]> spaceTimeDerivatives = new Dictionary<int, double[]>();
			for (int i = 0; i < current.Length; i++)
			{
				var temp = current[i];
				double[] temp1 = new double[] { temp[0], temp[1], temp[2] };
				var previousTemp = previous[i];
				double[] temp2 = new double[] { previousTemp[0], previousTemp[1], previousTemp[2] };
				temp1.Zip(temp2, (x, y) => (x - y) / timestep);
				spaceTimeDerivatives.Add(i, temp1);
			}
			return spaceTimeDerivatives;
		}
		private static void PressureCoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			var modelReader = SvDModel.Item2;
			foreach (Element element in modelReader.elementDomains[0])
			{

				l[element.ID] = dd0[element.ID] >= 0d ? lp[0] * SvDElement[element.ID] * Svin * dd0[element.ID] * 24d * 3600d : 0;
				u[element.ID] = conv0[0];
			}
			foreach (Element element in modelReader.elementDomains[1])
			{
				l[element.ID] = lp[1] * Svin * SvDElement[element.ID] * 24d * 3600d;
				u[element.ID] = conv0[1];
			}
			PressureL = l;
			PressureU = u;
		}
		private static void OxygenTransportCoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			OxygenTransportU = new Dictionary<int,double[]>(u);
			if (tumcElement == null)
			{
				tumcElement = new double[u.Count];
				for (int i = 0; i < u.Count; i++)
				{
					tumcElement[i] = 0.96;
				}
			}
			if (c_oxElement == null)
			{
				c_oxElement = new double[u.Count];
				for (int i = 0; i < u.Count; i++)
				{
					c_oxElement[i] = 0d;
				}
			}
			if (dd0 == null) dd0 = new double[SvDModel.Item2.elementDomains[0].Count];
			foreach (var element in SvDModel.Item2.elementDomains[0])
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
				var dd = (-2e-10) * tumcElement[element.ID] / (4 * Math.PI * Math.Pow(ri, 2)) + 1;
				dd0[element.ID] = dd >= 0 ? dd : 0;
			}
            foreach (var e in SvDModel.Item2.elementDomains[0])
            {
                //l[e.ID] = 24 * 3600 * Dox[0] / Lwv * Svin * dd0[e.ID] * SvDElement[e.ID];
                l[e.ID] = 24 * 3600 * (Dox[0] / Lwv * Svin * dd0[e.ID] * SvDElement[e.ID]
				+ ((Aox[0]) / (kox[0] + c_oxElement[e.ID] * cvox)) * (tumcElement[e.ID] + 0.3));
            }
			foreach (var e in SvDModel.Item2.elementDomains[1])
            {
                //l[e.ID] = 24 * 3600 * Dox[1] / Lwv * Svin * SvDElement[e.ID];
                l[e.ID] = 24 * 3600 * (Dox[1] / Lwv * Svin * SvDElement[e.ID]
					+ ((Aox[1]) / (kox[1] + c_oxElement[e.ID] * cvox)));
            }
			OxygenTransportL = l;
		}
		private static void TumorCellsCoefficientsCalculation(Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			if (ddcoxElement == null | dCsElement == null)
			{
				ddcoxElement = u;
				dcoxElement = u;
				dCsElement = u;
			}
			else
				for (int i = 0; i < u.Count; i++)
				{
					u[i] = new double[] { 24 * 3600 * mtox * (cvox * dcoxElement[i][0] + WvsTc * dCsElement[i][0] * Cs0),
						24 * 3600 * mtox * (cvox * dcoxElement[i][1] + WvsTc * dCsElement[i][1] * Cs0),
						24 * 3600 * mtox * (cvox * dcoxElement[i][2] + WvsTc * dCsElement[i][2] * Cs0)};
				}
			CancerTransportU = u;
			if (c_oxElement == null | CsElement == null)
			{
				c_oxElement = new double[u.Count];
				CsElement = new double[u.Count];
				for (int i = 0; i < u.Count; i++)
				{
					c_oxElement[i] = 0;
					CsElement[i] = 0;
				}
			}
			else
				for (int i = 0; i < u.Count; i++)
				{
					l[i] = 24 * 3600 * mtox * (cvox * ddcoxElement[i].Sum() + WvsTc * Cs0 * ddCsElement[i].Sum());
				}
			CancerTransportL = l;
		}
		private static void SvDCoefficientsCalculation(IList<IList<Element>> elementDomains, Dictionary<int, double> k, Dictionary<int, double[]> u, Dictionary<int, double> l)
		{
			if (SvDElement == null)
			{
				SvDElement = new double[u.Count];
                foreach (var e in elementDomains[0])
                {
					SvDElement[e.ID] = 0.5;
                }
				foreach (var e in elementDomains[1])
				{
					SvDElement[e.ID] = 1;
				}
			}
			for (int i = 0; i < u.Count; i++)
			{
				k[i] = 24 * 3600 * Dsv * Math.Pow((1+s1*a10),-aD) * Math.Pow((1+s2*a20),bD);
			}
			SvDK = k;
			if (dCsElement == null)
			{
				dCsElement = new Dictionary<int, double[]>(u);
				for (int i = 0; i < u.Count; i++)
				{
					u[i] = new double[] { 24 * 3600 * xn * Cv0,
						24 * 3600 * xn * Cv0,
						24 * 3600 * xn * Cv0};
				}
			}
			else
                for (int i = 0; i < u.Count; i++)
                {
                    u[i] = new double[] { 24 * 3600 * xn * (Cv0 + WvsSv * dCsElement[i][0] * Cs0),
                        24 * 3600 * xn * (Cv0 + WvsSv * dCsElement[i][1] * Cs0),
                        24 * 3600 * xn * (Cv0 + WvsSv * dCsElement[i][2] * Cs0)};
                }
            SvDU = u;
			if (ddCsElement == null)
				ddCsElement = new Dictionary<int, double[]>(u);
            else
                for (int i = 0; i < u.Count; i++)
                {
                    l[i] = 24 * 3600 * xn * WvsSv * Cs0 * ddCsElement[i].Sum();
                }
            SvDL = l;
		}
		private static void UpdateModels(Dictionary<int, IVector>[] prevStepSolutions, IStructuralModel[] modelsToReplace, ISolver[] solversToReplace,
			IConvectionDiffusionIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
		{
			SvDNode = solversToReplace[0].LinearSystems[0].Solution.CopyToArray();
			c_oxNode = solversToReplace[1].LinearSystems[0].Solution.CopyToArray();
			tumcNode = solversToReplace[2].LinearSystems[0].Solution.CopyToArray();
			lgNode = solversToReplace[3].LinearSystems[0].Solution.CopyToArray();
			pSolution = solversToReplace[4].LinearSystems[0].Solution.CopyToArray();
			if (pNode == null) pNode = new double[modelsToReplace[4].Nodes.Count];
			int freeDofNo = 0;
			for (int i = 0; i < modelsToReplace[4].Nodes.Count; i++)
			{
				if (modelsToReplace[4].GlobalDofOrdering.GlobalFreeDofs.Contains(modelsToReplace[4].Nodes[i], ThermalDof.Temperature))
				{
					pNode[i] = pSolution[freeDofNo];
					freeDofNo++;
				}
			}
			CsNode = solversToReplace[5].LinearSystems[0].Solution.CopyToArray();

			//SvD
			if (SvDElement == null) SvDElement = new double[modelsToReplace[0].Elements.Count];
			foreach (var e in modelsToReplace[0].Elements)
			{
				SvDElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					SvDElement[e.ID] += SvDNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}
			//c_ox
			var c_oxSubdomain = solversToReplace[1].LinearSystems[0].Subdomain;
			var c_oxFirstDerivatives = providersToReplace[1].GetFirstSpaceDerivatives(c_oxSubdomain, Vector.CreateFromArray(c_oxNode));
			var c_oxSecondDerivatives = providersToReplace[1].GetSecondSpaceDerivatives(c_oxSubdomain, Vector.CreateFromArray(c_oxNode));
			for (int i = 0; i < c_oxSecondDerivatives.NumRows; i++)
			{
				dcoxNode[i] = c_oxFirstDerivatives.GetRow(i).CopyToArray();
				ddcoxNode[i] = c_oxSecondDerivatives.GetRow(i).CopyToArray();
			}
			if (c_oxElement == null) c_oxElement = new double[modelsToReplace[1].Elements.Count];
			foreach (var e in modelsToReplace[1].Elements)
			{
				c_oxElement[e.ID] = 0;
				dcoxElement[e.ID] = new double[] { 0, 0, 0 };
				ddcoxElement[e.ID] = new double[] { 0, 0, 0 };
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					c_oxElement[e.ID] += c_oxNode[i] / (e.Nodes.Count);
					for (int j = 0; j < 3; j++)
					{
						dcoxElement[e.ID][j] += dcoxNode[i][j] / (e.Nodes.Count);
						ddcoxElement[e.ID][j] += ddcoxNode[i][j] / (e.Nodes.Count);
					}
				}
			}
			//Tumc
			if (tumcElement == null) tumcElement = new double[modelsToReplace[2].Elements.Count];
			foreach (var e in modelsToReplace[2].Elements)
			{
				tumcElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					tumcElement[e.ID] += tumcNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}
			//lGrowth
			if (lgElement == null) lgElement = new double[modelsToReplace[3].Elements.Count];
			foreach (var e in modelsToReplace[3].Elements)
			{
				lgElement[e.ID] = 0;
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					lgElement[e.ID] += lgNode[e.Nodes[i].ID] / (e.Nodes.Count);
				}
			}
			//pressure
			var pSubdomain = solversToReplace[4].LinearSystems[0].Subdomain;
			var pFirstDerivatives = providersToReplace[4].GetFirstSpaceDerivatives(pSubdomain, Vector.CreateFromArray(pSolution));
			var pSecondDerivatives = providersToReplace[4].GetSecondSpaceDerivatives(pSubdomain, Vector.CreateFromArray(pSolution));
			for (int i = 0; i < pFirstDerivatives.NumRows; i++)
			{
				dpSolution[i] = pFirstDerivatives.GetRow(i).CopyToArray();
				ddpSolution[i] = pSecondDerivatives.GetRow(i).CopyToArray();
			}
			freeDofNo = 0;
			for (int i = 0; i < modelsToReplace[4].Nodes.Count; i++)
			{
				if (modelsToReplace[4].GlobalDofOrdering.GlobalFreeDofs.Contains(modelsToReplace[4].Nodes[i], ThermalDof.Temperature))
				{
					dpNode[i] = dpSolution[freeDofNo];
					ddpNode[i] = ddpSolution[freeDofNo];
					freeDofNo++;
				}
			}
			if (pElement == null) pElement = new double[modelsToReplace[4].Elements.Count];
			if (dpElement == null) dpElement = new Dictionary<int, double[]>();
			if (ddpElement == null) ddpElement = new Dictionary<int, double[]>();
			foreach (var e in modelsToReplace[4].Elements)
			{
				pElement[e.ID] = 0;
				dpElement[e.ID] = new double[] { 0, 0, 0 };
				ddpElement[e.ID] = new double[] { 0, 0, 0 };
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					pElement[e.ID] += pNode[i] / (e.Nodes.Count);
					for (int j = 0; j < 3; j++)
					{
						dpElement[e.ID][j] += dpNode[i][j] / (e.Nodes.Count);
						ddpElement[e.ID][j] += ddpNode[i][j] / (e.Nodes.Count);
					}
				}
			}
			//Cs
			var CsSubdomain = solversToReplace[5].LinearSystems[0].Subdomain;
			var CsFirstDerivatives = providersToReplace[5].GetFirstSpaceDerivatives(CsSubdomain, Vector.CreateFromArray(CsNode));
			var CsSecondDerivatives = providersToReplace[5].GetSecondSpaceDerivatives(CsSubdomain, Vector.CreateFromArray(CsNode));
			for (int i = 0; i < CsSecondDerivatives.NumRows; i++)
			{
				dCsNode[i] = CsFirstDerivatives.GetRow(i).CopyToArray();
				ddCsNode[i] = CsSecondDerivatives.GetRow(i).CopyToArray();
			}
			if (CsElement == null) CsElement = new double[modelsToReplace[5].Elements.Count];
			if (dCsElement == null) dCsElement = new Dictionary<int, double[]>();
			if (ddCsElement == null) ddCsElement = new Dictionary<int, double[]>();
			foreach (var e in modelsToReplace[5].Elements)
			{
				CsElement[e.ID] = 0;
				dCsElement[e.ID] = new double[] { 0, 0, 0 };
				ddCsElement[e.ID] = new double[] { 0, 0, 0 };
				for (int i = 0; i < e.Nodes.Count; i++)
				{
					CsElement[e.ID] += CsNode[e.Nodes[i].ID] / (e.Nodes.Count);
					for (int j = 0; j < 3; j++)
					{
						dCsElement[e.ID][j] += dCsNode[i][j] / (e.Nodes.Count);
						ddCsElement[e.ID][j] += ddCsNode[i][j] / (e.Nodes.Count);
					}
				}
			}

			//
			modelsToReplace[0] = CreateSvDModel().Item1;
			modelsToReplace[1] = CreateOxygenTransportModel(new double[] { Dox[0] * 24d * 3600d, Dox[1] * 24d * 3600d }).Item1;
			modelsToReplace[2] = CreateCancerTransportModel(Dcell[0]).Item1;
			modelsToReplace[3] = CreateGrowthModel().Item1;
			modelsToReplace[4] = CreatePressureModel(khy).Item1;
			modelsToReplace[5] = CreateCsModel().Item1;

			for (int i = 0; i < modelsToReplace.Length; i++)
			{
				if (i==0 || i == 2)
				{
					ctBuilder.IsMatrixPositiveDefinite = false;
					solversToReplace[i] = ctBuilder.BuildSolver(modelsToReplace[i]);
				}
				else
					solversToReplace[i] = builder.BuildSolver(modelsToReplace[i]);
				providersToReplace[i] = new ProblemConvectionDiffusion2((Model)modelsToReplace[i], solversToReplace[i]);
				childAnalyzersToReplace[i] = new LinearAnalyzer(modelsToReplace[i], solversToReplace[i], providersToReplace[i]);
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
		private static void UpdateNewmarkModel(Dictionary<int, IVector> accelerations, Dictionary<int, IVector> velocities, Dictionary<int, IVector> displacements, IStructuralModel[] modelsToReplace,
			ISolver[] solversToReplace, IImplicitIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
		{
			Accelerations = accelerations;
			Velocities = velocities;
			Displacements = displacements;
			if (lgElement == null) lgElement = new double[modelsToReplace[0].Elements.Count];
			foreach (Element e in structModel.Item2.elementDomains[1])
			{
				lgElement[e.ID] = 1d;
			}
			//if (accNode == null) accNode = new double[modelsToReplace[3].Nodes.Count];
			int freeDofNo = 0;
			for (int i = 0; i < structModel.Item1.Nodes.Count; i++)
			{
				aNode[i] = new double[3];
				vNode[i] = new double[3];
				uNode[i] = new double[3];
				if (structModel.Item1.GlobalDofOrdering.GlobalFreeDofs.Contains(structModel.Item1.Nodes[i], StructuralDof.TranslationX))
				{
					aNode[i][0] = Accelerations[0][freeDofNo];
					vNode[i][0] = Velocities[0][freeDofNo];
					uNode[i][0] = Displacements[0][freeDofNo];
					freeDofNo++;
				}
				if (structModel.Item1.GlobalDofOrdering.GlobalFreeDofs.Contains(structModel.Item1.Nodes[i], StructuralDof.TranslationY))
				{
					aNode[i][1] = Accelerations[0][freeDofNo];
					vNode[i][1] = Velocities[0][freeDofNo];
					uNode[i][1] = Displacements[0][freeDofNo];
					freeDofNo++;
				}
				if (structModel.Item1.GlobalDofOrdering.GlobalFreeDofs.Contains(structModel.Item1.Nodes[i], StructuralDof.TranslationZ))
				{
					aNode[i][2] = Accelerations[0][freeDofNo];
					vNode[i][2] = Velocities[0][freeDofNo];
					uNode[i][2] = Displacements[0][freeDofNo];
					freeDofNo++;
				}
			}

			foreach (Element e in structModel.Item1.Elements)
			{
				aElement[e.ID] = new double[3];
				vElement[e.ID] = new double[3];
				uElement[e.ID] = new double[3];
				foreach (Node n in e.Nodes)
				{
					for (int i = 0; i < 3; i++)
					{
						aElement[e.ID][i] += aNode[n.ID][i] / e.Nodes.Count;
						vElement[e.ID][i] += vNode[n.ID][i] / e.Nodes.Count;
						uElement[e.ID][i] += uNode[n.ID][i] / e.Nodes.Count;
					}
				}
			}
			foreach (Element e in structModel.Item2.elementDomains[1])
			{
				lgElement[e.ID] = 1d;
			}
			ReplaceLambdaGInModel(modelsToReplace[0], lgElement);
			UpdateStructuralLoads();
			solversToReplace[0] = structuralBuilder.BuildSolver(modelsToReplace[0]);
			providersToReplace[0] = new ProblemStructural(modelsToReplace[0], solversToReplace[0]);
			//solversToReplace[0].HandleMatrixWillBeSet();
			//childAnalyzersToReplace[0] = new LinearAnalyzer(modelsToReplace[0], solversToReplace[0], providersToReplace[0]);
			var increments = 2;
			var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(modelsToReplace[0], solversToReplace[0], (INonLinearProvider)providersToReplace[0], increments);
			childAnalyzerBuilder.ResidualTolerance = 1E-5;
			childAnalyzerBuilder.MaxIterationsPerIncrement = 50;
			childAnalyzerBuilder.NumIterationsForMatrixRebuild = 5;
			childAnalyzersToReplace[0] = childAnalyzerBuilder.Build();
		}
		private static Tuple<Model, IModelReader> CreateSvDModel()
		{
			ComsolMeshReader5 modelReader;
			Model model;
			if (SvDModel == null)
			{
				Console.WriteLine("Creating SvD Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "mesh446elem.mphtxt");
				int[] modelDomains = new int[] { 0, 1 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
				modelReader = new ComsolMeshReader5(filename, new double[] { 1, 1}, SvDCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating SvD Model...");
				modelReader = (ComsolMeshReader5)SvDModel.Item2;
				SvDCoefficientsCalculation(modelReader.elementDomains, SvDK, SvDU, SvDL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1, 1}, SvDK, SvDU, SvDL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}
			//double[] Grox = new double[model.Elements.Count];
			if (fsv == null)
			{
				fsv = new Dictionary<int, double>();
				for (int i = 0; i < model.Elements.Count; i++)
				{
					fsv[i] = 0d;
				}
			}
			int[] domainIDs = new int[] { 0, 1 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					var svdMaterial = new ConvectionDiffusionMaterial(1, SvDK[element.ID], SvDU[element.ID], SvDL[element.ID]);
					if (SvDElement[element.ID] >= 0.5)
						fsv[element.ID] = 24 * 3600 * (l2 * Cv0 - l4 * Svin);
					else
						fsv[element.ID] = 0;
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(svdMaterial, fsv[element.ID], ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateOxygenTransportModel(double[] k)
		{
			ComsolMeshReader4 modelReader;
			Model model;

			if (oxModel == null)
			{
				Console.WriteLine("Creating Oxygen Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "mesh446elem.mphtxt");
				int[] modelDomains = new int[] { 0, 1 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 1, 1 }, k, OxygenTransportCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Oxygen Model...");
				modelReader = (ComsolMeshReader4)oxModel.Item2;
				OxygenTransportCoefficientsCalculation(OxygenTransportU, OxygenTransportL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1, 1}, k, OxygenTransportU, OxygenTransportL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}
			if (c_oxElement == null)
			{
				c_oxElement = new double[model.Elements.Count];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					c_oxElement[i] = 0;/* 0.9673;*/
				}
			}
			if (tumcElement == null)
			{
				tumcElement = new double[model.Elements.Count];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					tumcElement[i] = 0.96;/* 0.9673;*/
				}
			}
			if (SvDElement == null)
			{
				SvDElement = new double[model.Elements.Count];
				foreach (var e in modelReader.elementDomains[0])
				{
					SvDElement[e.ID] = 0.5;
				}
				foreach (var e in modelReader.elementDomains[1])
				{
					SvDElement[e.ID] = 1;
				}
			}

			if (dd0 == null) dd0 = new double[modelReader.elementDomains[0].Count];
			foreach (var element in modelReader.elementDomains[0])
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
				dd0[element.ID] = (-2e-10) * tumcElement[element.ID] / (4 * Math.PI * Math.Pow(ri, 2)) + 1;
			}
			if (dcoxdvf == null) dcoxdvf = new double[model.Elements.Count];

			if (uXt == null) uXt = new Dictionary<int, double[]>();
			for (int i = 0; i < model.Elements.Count; i++)
			{
				uXt[i] = new double[6];
			}

			if (dpElement == null)
			{
				dpElement = new Dictionary<int, double[]>();
				for (int i = 0; i < model.Elements.Count; i++)
				{
					dpElement[i] = new double[3];
				}
			}

			if (ddpElement == null)
			{
				ddpElement = new Dictionary<int, double[]>();
				for (int i = 0; i < model.Elements.Count; i++)
				{
					ddpElement[i] = new double[3];
				}
			}
			if (dcoxElement == null)
			{
				dcoxElement = new Dictionary<int, double[]>();
				for (int i = 0; i < model.Elements.Count; i++)
				{
					dcoxElement[i] = new double[3];
				}
			}

			if (vElement == null) vElement = new Dictionary<int, double[]>();
			for (int i = 0; i < model.Elements.Count; i++)
			{
				vElement[i] = new double[6];
			}
			double[] fox = new double[model.Elements.Count];
			int[] domainIDs = new int[] { 0, 1 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					uof[element.ID] = -khy[0] * dpElement[element.ID][0] / .7 + vElement[element.ID][0];
					vof[element.ID] = -khy[0] * dpElement[element.ID][1] / .7 + vElement[element.ID][1];
					wof[element.ID] = -khy[0] * dpElement[element.ID][2] / .7 + vElement[element.ID][2];
					duof[element.ID] = -khy[0] * ddpElement[element.ID][0] / .7 + uXt[element.ID][0];
					dvof[element.ID] = -khy[0] * ddpElement[element.ID][1] / .7 + uXt[element.ID][1];
					dwof[element.ID] = -khy[0] * ddpElement[element.ID][2] / .7 + uXt[element.ID][2];
					dcoxdvf[element.ID] = dcoxElement[element.ID][0] * uof[element.ID] +
						dcoxElement[element.ID][1] * vof[element.ID] +
						dcoxElement[element.ID][2] * wof[element.ID] + c_oxElement[element.ID] * (duof[element.ID] + dvof[element.ID] + dwof[element.ID]);
					var material = new ConvectionDiffusionMaterial(1, k[domainID], OxygenTransportU[element.ID], OxygenTransportL[element.ID]);
					if (domainID == 0)
					{
                        //                 fox[element.ID] = (Dox[domainID] / Lwv * SvDElement[element.ID] * Svin * dd0[element.ID] * (1 - c_oxElement[element.ID])
                        //- ((Aox[domainID] * c_oxElement[element.ID]) / (kox[domainID] + c_oxElement[element.ID] * cvox)) * (tumcElement[element.ID] + 0.3)
                        //                     - dcoxdvf[element.ID]) * (24d * 3600d);
                        //fox[element.ID] = (Dox[domainID] / Lwv * Svin * dd0[element.ID] * SvDElement[element.ID]
                        //    - ((Aox[domainID] * c_oxElement[element.ID]) / (kox[domainID] + c_oxElement[element.ID] * cvox)) * (tumcElement[element.ID] + 0.3)
                        //    - dcoxdvf[element.ID]) * (24d * 3600d);
                        fox[element.ID] = (Dox[domainID] / Lwv * Svin * dd0[element.ID] * SvDElement[element.ID]
                            - dcoxdvf[element.ID]) * (24d * 3600d);
                    }
					else
					{
                        //fox[element.ID] = (Dox[domainID] / Lwv * SvDElement[element.ID] * Svin * (1 - c_oxElement[element.ID])
                        //    - (Aox[domainID] * c_oxElement[element.ID]) / (kox[domainID] + c_oxElement[element.ID] * cvox)
                        //    - dcoxdvf[element.ID]) * (24d * 3600d);
                        //               fox[element.ID] = (Dox[domainID] / Lwv * Svin * SvDElement[element.ID]
                        //- (Aox[domainID] * c_oxElement[element.ID]) / (kox[domainID] + c_oxElement[element.ID] * cvox)
                        //                  - dcoxdvf[element.ID]) * (24d * 3600d);
                        fox[element.ID] = (Dox[domainID] / Lwv * Svin * SvDElement[element.ID]
							- dcoxdvf[element.ID]) * (24d * 3600d);
                    }
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(material, fox[element.ID], ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateGrowthModel()
		{
			ComsolMeshReader3 modelReader;
			Model model;

			if (gModel == null)
			{
				Console.WriteLine("Creating Growth Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "mesh446elem.mphtxt");
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader3(filename, new double[] { 1 }, new double[] { 0 }, conv0, new double[] { 0 });
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Growth Model...");
				modelReader = (ComsolMeshReader3)gModel.Item2;
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { 0 }, conv0, new double[] { 0 });
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

			var materialODE = new ConvectionDiffusionMaterial(1, 0, conv0[0], 0);
			//double[] Grox = new double[model.Elements.Count];
			double[] fg = new double[model.Elements.Count];
			int[] domainIDs = new int[] { 0, };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					Grox[element.ID] = (loxc[domainID] * cvox * c_oxElement[element.ID]) / (cvox * c_oxElement[element.ID] + Koxc[domainID]);
					fg[element.ID] = 24d * 3600d * Grox[element.ID] * lgElement[element.ID] / 3d;
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, fg[element.ID], ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateCancerTransportModel(double k)
		{
			ComsolMeshReader4 modelReader;
			Model model;
			if (ctModel == null)
			{
				Console.WriteLine("Creating Cancer Transport Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "mesh446elem.mphtxt");
				int[] modelDomains = new int[] { 0 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 5 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 1 }, new double[] { k }, TumorCellsCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Cancer Transport Model...");
				modelReader = (ComsolMeshReader4)ctModel.Item2;
				TumorCellsCoefficientsCalculation(CancerTransportU, CancerTransportL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1 }, new double[] { k }, CancerTransportU, CancerTransportL);
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}
			//double[] Grox = new double[model.Elements.Count];
			if (RTumc == null)
			{
				RTumc = new Dictionary<int, double>();
				for (int i = 0; i < model.Elements.Count; i++)
				{
					RTumc[i] = 0d;
				}
			}
			int[] domainIDs = new int[] { 0, };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					var CTmaterial = new ConvectionDiffusionMaterial(1, k, CancerTransportU[element.ID], CancerTransportL[element.ID]);
					Grox[element.ID] = (loxc[domainID] * cvox * c_oxElement[element.ID]) / (cvox * c_oxElement[element.ID] + Koxc[domainID]);
					RTumc[element.ID] = 24d * 3600d * Grox[element.ID];
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(CTmaterial, RTumc[element.ID], ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
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
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "mesh446elem.mphtxt");
				int[] modelDomains = new int[] { 0, 1 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
				modelReader = new ComsolMeshReader4(filename, new double[] { 1, 1 },
					k, PressureCoefficientsCalculation);
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);

			}
			else
			{
				Console.WriteLine("Updating Pressure Model...");
				modelReader = (ComsolMeshReader4)prModel.Item2;
				PressureCoefficientsCalculation(PressureU, PressureL);
				modelReader = modelReader.UpdateModelReader(new double[] { 1 , 1}, 
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
			if (PreviousSpaceDerivatives == null)
			{
				PreviousSpaceDerivatives = new double[model.Elements.Count][];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					PreviousSpaceDerivatives[i] = new double[6];
				}
			}
			if (SpaceDerivatives == null)
			{
				SpaceDerivatives = new double[model.Elements.Count][];
				for (int i = 0; i < model.Elements.Count; i++)
				{
					SpaceDerivatives[i] = new double[6];
				}
			}
			if (PressureL == null)
			{
				for (int i = 0; i < model.Elements.Count; i++)
				{
					PressureL.Add(i, 0d);
				}
			}
			PreviousSpaceDerivatives = SpaceDerivatives.Clone() as double[][];
			SpaceDerivatives = GetStrains(model.Elements.Count);
			uXt = StructuralSpaceTimeDerivatives(SpaceDerivatives, PreviousSpaceDerivatives);

			int[] domainIDs = new int[] { 0, 1 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					var prMaterial = new ConvectionDiffusionMaterial(1, k[domainID], conv0[0], PressureL[element.ID]);
					Grox[element.ID] = (loxc[domainID] * cvox * c_oxElement[element.ID]) / (cvox * c_oxElement[element.ID] + Koxc[domainID]);
					RTumc[element.ID] = 24d * 3600d * Grox[element.ID];
					var fp = RTumc[element.ID] + pv * PressureL[element.ID] - uXt[element.ID].Sum()*24*3600;
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(prMaterial, fp, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
			}
			return new Tuple<Model, IModelReader>(model, modelReader);
		}
		private static Tuple<Model, IModelReader> CreateCsModel()
		{
			ComsolMeshReader3 modelReader;
			Model model;

			if (csModel == null)
			{
				Console.WriteLine("Creating Cs Model");
				string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "mesh446elem.mphtxt");
				int[] modelDomains = new int[] { 0, 1 };
				int[] modelBoundaries = new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
				modelReader = new ComsolMeshReader3(filename, new double[] { 1, 1}, new double[] { 0, 0 }, conv0, new double[] { l13*24*3600, 0 });
				model = modelReader.CreateModelFromFile(modelDomains, modelBoundaries);
			}
			else
			{
				Console.WriteLine("Updating Cs Model...");
				modelReader = (ComsolMeshReader3)csModel.Item2;
				modelReader = modelReader.UpdateModelReader(new double[] { 1, 1}, new double[] { 0, 0 }, conv0, new double[] { l13*24*3600, 0 });
				model = modelReader.UpdateModel(structModel.Item1, Displacements);
			}

			var materialODE = new ConvectionDiffusionMaterial(1, 0, conv0[0], l13*24*3600);
			int[] domainIDs = new int[] { 0, };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					double Ga = 0d;
					if (c_oxElement[element.ID] <= 0.5 & c_oxElement[element.ID] > 0d) Ga = 3;
					else if (c_oxElement[element.ID] > 0.5 & c_oxElement[element.ID] <= 1) Ga = 2 - c_oxElement[element.ID];
					else if (c_oxElement[element.ID] > 1d) Ga = 0.5 * c_oxElement[element.ID];
					var Cs_stD = 24*3600*((l10 * Ga * tumcElement[element.ID] * T0 * 100 / Cs0) / 2 + 
						(l10 * SvDElement[element.ID] * 100 * Cv0) / 2);
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var domainLoad = new ConvectionDiffusionDomainLoad(materialODE, Cs_stD, ThermalDof.Temperature);
					var bodyLoadElementFactory = new BodyLoadElementFactory(domainLoad, model);
					var bodyLoadElement = bodyLoadElementFactory.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElement);
				}
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
			string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "mesh446elem.mphtxt");
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
			if (dpElement == null)
			{
				dpElement = new Dictionary<int, double[]>();
				for (int i = 0; i < model.Elements.Count; i++)
				{
					dpElement[i] = new double[3];
				}
			}
			int[] domainIDs = new int[] { 0, 1 };
			foreach (int domainID in domainIDs)
			{
				foreach (Element element in modelReader.elementDomains[domainID])
				{
					var nodes = (IReadOnlyList<Node>)element.Nodes;
					var bodyLoadX = new GravityLoad(1d, -dpElement[element.ID][0], StructuralDof.TranslationX);
					var bodyLoadElementFactoryX = new BodyLoadElementFactory(bodyLoadX, model);
					var bodyLoadElementX = bodyLoadElementFactoryX.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElementX);
					var bodyLoadY = new GravityLoad(1d, -dpElement[element.ID][0], StructuralDof.TranslationY);
					var bodyLoadElementFactoryY = new BodyLoadElementFactory(bodyLoadY, model);
					var bodyLoadElementY = bodyLoadElementFactoryY.CreateElement(CellType.Tet4, nodes);
					model.BodyLoads.Add(bodyLoadElementY);
					var bodyLoadZ = new GravityLoad(1d, -dpElement[element.ID][0], StructuralDof.TranslationZ);
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
					var bodyLoadX = new GravityLoad(1d, -dpElement[element.ID][0], StructuralDof.TranslationX);
					var bodyLoadElementFactoryX = new BodyLoadElementFactory(bodyLoadX, structModel.Item1);
					var bodyLoadElementX = bodyLoadElementFactoryX.CreateElement(CellType.Tet4, nodes);
					structModel.Item1.BodyLoads.Add(bodyLoadElementX);
					var bodyLoadY = new GravityLoad(1d, -dpElement[element.ID][0], StructuralDof.TranslationY);
					var bodyLoadElementFactoryY = new BodyLoadElementFactory(bodyLoadY, structModel.Item1);
					var bodyLoadElementY = bodyLoadElementFactoryY.CreateElement(CellType.Tet4, nodes);
					structModel.Item1.BodyLoads.Add(bodyLoadElementY);
					var bodyLoadZ = new GravityLoad(1d, -dpElement[element.ID][0], StructuralDof.TranslationZ);
					var bodyLoadElementFactoryZ = new BodyLoadElementFactory(bodyLoadZ, structModel.Item1);
					var bodyLoadElementZ = bodyLoadElementFactoryZ.CreateElement(CellType.Tet4, nodes);
					structModel.Item1.BodyLoads.Add(bodyLoadElementZ);
				}
			}
		}
		private static IVectorView[] SolveModelsWithNewmark(Model[] models, IModelReader[] modelReaders)
		{
			double[] muLame = new double[] { 6e4, 2.1e4 };
			double[] poissonV = new double[] { .45, .2 };
			IDynamicMaterial[] dynamicMaterials = new DynamicMaterial[] { new DynamicMaterial(.001, 0, 0, true), new DynamicMaterial(.001, 0, 0, true) };
			structModel = CreateStructuralModel(muLame, poissonV, dynamicMaterials, 0, new double[] { 0, 0, 0 }, lgElement);//.Item1; // new Model();
			var structuralModel = structModel.Item1;
			var structuralSolver = structuralBuilder.BuildSolver(structuralModel);
			var structuralProvider = new ProblemStructural(structuralModel, structuralSolver);
			//var structuralChildAnalyzer = new LinearAnalyzer(structuralModel, structuralSolver, structuralProvider);
			var increments = 2;
			var structuralChildAnalyzerBuilder = new LoadControlAnalyzer.Builder(structuralModel, structuralSolver, structuralProvider, increments);
			structuralChildAnalyzerBuilder.ResidualTolerance = 1E-5;
			structuralChildAnalyzerBuilder.MaxIterationsPerIncrement = 50;
			structuralChildAnalyzerBuilder.NumIterationsForMatrixRebuild = 5;
			//childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
			LoadControlAnalyzer structuralChildAnalyzer = structuralChildAnalyzerBuilder.Build();
			var structuralParentAnalyzer = new NewmarkDynamicAnalyzer(UpdateNewmarkModel, structuralModel, structuralSolver,
				structuralProvider, structuralChildAnalyzer, timestep, time, 0.25, 0.5);
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
				value0[0][node.ID] = 1;
			}
			foreach (Node node in modelReaders[0].nodeDomains[0])
			{
				value0[0][node.ID] = 0.5;
			}
			foreach (Node node in models[1].Nodes)
			{
				value0[1][node.ID] = 0; /* 0.96733;*/
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

			IConvectionDiffusionIntegrationProvider[] providers = new IConvectionDiffusionIntegrationProvider[models.Length];
			IChildAnalyzer[] childAnalyzers = new IChildAnalyzer[models.Length];
			var solvers = new ISolver[models.Length];
			for (int i = 0; i < models.Length; i++)
			{
				initialValues[i] = Vector.CreateFromArray(value0[i]);
				//var builder = new DenseMatrixSolver.Builder();
				//builder.IsMatrixPositiveDefinite = false;
				if (i==0 || i == 2)
				{
					ctBuilder.IsMatrixPositiveDefinite = false;
					solvers[i] = ctBuilder.BuildSolver(models[i]);
				}
				else
					solvers[i] = builder.BuildSolver(models[i]);
				providers[i] = new ProblemConvectionDiffusion2(models[i], solvers[i]);
				childAnalyzers[i] = new LinearAnalyzer(models[i], solvers[i], providers[i]);
			}

			var parentAnalyzer = new ConvectionDiffusionImplicitDynamicAnalyzerMultiModel(UpdateModels, models, solvers,
				providers, childAnalyzers, timestep, time, structuralParentAnalyzer: structuralParentAnalyzer, initialTemperature: initialValues);
			//var parentAnalyzer = new ConvectionDiffusionImplicitDynamicAnalyzerMultiModel(UpdateModels, models, solvers,
			//	providers, childAnalyzers, timestep, time, initialTemperature: initialValues);
			parentAnalyzer.Initialize();

			for (int i = 0; i < time / timestep; i++)
			{
				parentAnalyzer.SolveTimestep(i);
				//structuralParentAnalyzer.SolveTimestep(i);
				Paraview(i);
			}

			return solvers.Select(x => x.LinearSystems[subdomainID].Solution).ToArray();
		}
	}
}