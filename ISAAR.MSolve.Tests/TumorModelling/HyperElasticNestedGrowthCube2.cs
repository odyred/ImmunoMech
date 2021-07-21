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
using System.Globalization;

namespace ISAAR.MSolve.Tests
{
    public class HyperElasticNestedGrowthCube2
    {
        private const int subdomainID = 0;
        private static readonly double[] loxc = new double[] { .07 / 24 / 3600, 1.0 / 24 / 3600 }; //1/s
        private static readonly double[] Aox = new double[] { 2200.0 / 24 / 3600, 2200.0 / 24 / 3600 }; //mol/(m^3*s)
        private static readonly double[] Dox = new double[] { 1.78e-9, 1.79e-9 }; //m^2/s
        private static readonly double[] kox = new double[] { .00464, .00464 }; //mol/m^3
        private static readonly double[] Koxc = new double[] { 0.0083, 0.0083 }; //mol/m^3
        private static double cvox = 0.2; //mol/m^3
        private static double Lwv = 5e-6; //m
        private static double[][] conv0 = new double[][] { new double[] { 0, 0, 0 }, new double[] { 0, 0, 0 } };
        //private static double fox = -((Aox * c_ox) / (kox + c_ox * cvox)) * 0.3;
        private static SkylineSolver.Builder builder = new SkylineSolver.Builder();
        //private static DenseMatrixSolver.Builder builder = new DenseMatrixSolver.Builder();
        private static SkylineSolver.Builder structuralBuilder = new SkylineSolver.Builder();
        private static double[] lgNode;
        private static double[] lgElement;
        private static double[] lg;
        //private static IModelReader structuralMR;
        private static Dictionary<int, IVector> Accelerations;
        private static Dictionary<int, IVector> Velocities;
        private static Dictionary<int, IVector> Displacements;
        private static Dictionary<int, double[]> aNode = new Dictionary<int, double[]>();
        private static Dictionary<int, double[]> aElement = new Dictionary<int, double[]>();
        private static Dictionary<int, double[]> vNode = new Dictionary<int, double[]>();
        private static Dictionary<int, double[]> vElement = new Dictionary<int, double[]>();
        private static Dictionary<int, double[]> uNode = new Dictionary<int, double[]>();
        private static Dictionary<int, double[]> uElement = new Dictionary<int, double[]>();
        private static Tuple<Model, IModelReader> oxModel, gModel, structModel;
        private static double day;
        public static double[][] Strains;
        public static double[][] Stresses;
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
            CreateGrowthModel();
            IVectorView solutions = SolveModelsWithNewmark();

            Assert.True(CompareResults(solutions));
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

                outputFile.WriteLine("          </PointData>");
                outputFile.WriteLine("          <CellData>");
                outputFile.WriteLine("              <DataArray type=\"Int32\" Name=\"element_ID\" NumberOfComponents=\"1\" format=\"ascii\">");
                for (int i = 0; i < numberOfCells; i++)
                {
                    outputFile.WriteLine($"{i + 1}");
                }
                outputFile.WriteLine("              </DataArray>");

                outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"uX\" NumberOfComponents=\"1\" format=\"ascii\">");
                for (int i = 0; i < structModel.Item1.Elements.Count; i++)
                    outputFile.WriteLine($"{Strains[i][0]}");
                outputFile.WriteLine("              </DataArray>");

                outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"vY\" NumberOfComponents=\"1\" format=\"ascii\">");
                for (int i = 0; i < structModel.Item1.Elements.Count; i++)
                    outputFile.WriteLine($"{Strains[i][1]}");
                outputFile.WriteLine("              </DataArray>");

                outputFile.WriteLine("              <DataArray type=\"Float64\" Name=\"wZ\" NumberOfComponents=\"1\" format=\"ascii\">");
                for (int i = 0; i < structModel.Item1.Elements.Count; i++)
                    outputFile.WriteLine($"{Strains[i][2]}");
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
            Strains = GetStrainsStresses(structModel.Item1.Elements.Count).Item1;
            Stresses = GetStrainsStresses(structModel.Item1.Elements.Count).Item2;
            ReplaceLambdaGInModel(modelsToReplace[0], lgElement);
            solversToReplace[0] = structuralBuilder.BuildSolver(modelsToReplace[0]);
            providersToReplace[0] = new ProblemStructural(modelsToReplace[0], solversToReplace[0]);
            //solversToReplace[0].HandleMatrixWillBeSet();
            //childAnalyzersToReplace[0] = new LinearAnalyzer(modelsToReplace[0], solversToReplace[0], providersToReplace[0]);
            var increments = 2;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(modelsToReplace[0], solversToReplace[0], (INonLinearProvider)providersToReplace[0], increments);
            childAnalyzerBuilder.ResidualTolerance = 1E-5;
            childAnalyzerBuilder.MaxIterationsPerIncrement = 50;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = 2;
            childAnalyzersToReplace[0] = childAnalyzerBuilder.Build();
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

        private static void CreateGrowthModel()
        {
            string[] lgString = System.IO.File.ReadAllLines(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "lg1.txt"));

            lg = new double[lgString.Length];
            for (int i = 0; i < lgString.Length; i++)
            {
                lg[i] = Double.Parse(lgString[i], CultureInfo.InvariantCulture);
            }
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

            ComsolMeshReader1 modelReader;
            string filename = Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "mesh446elem.mphtxt");
            if (lambdag == null)
                modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties);
            else
                modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties, lambdag);
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
            int[] domainIDs = new int[] { 0, 1 };
            foreach (int domainID in domainIDs)
            {
                foreach (Element element in modelReader.elementDomains[domainID])
                {
                    var nodes = (IReadOnlyList<Node>)element.Nodes;
                    var bodyLoadX = new GravityLoad(1d, -1d, StructuralDof.TranslationX);
                    var bodyLoadElementFactoryX = new BodyLoadElementFactory(bodyLoadX, model);
                    var bodyLoadElementX = bodyLoadElementFactoryX.CreateElement(CellType.Tet4, nodes);
                    model.BodyLoads.Add(bodyLoadElementX);
                    var bodyLoadY = new GravityLoad(1d, -1d, StructuralDof.TranslationY);
                    var bodyLoadElementFactoryY = new BodyLoadElementFactory(bodyLoadY, model);
                    var bodyLoadElementY = bodyLoadElementFactoryY.CreateElement(CellType.Tet4, nodes);
                    model.BodyLoads.Add(bodyLoadElementY);
                    var bodyLoadZ = new GravityLoad(1d, -1d, StructuralDof.TranslationZ);
                    var bodyLoadElementFactoryZ = new BodyLoadElementFactory(bodyLoadZ, model);
                    var bodyLoadElementZ = bodyLoadElementFactoryZ.CreateElement(CellType.Tet4, nodes);
                    model.BodyLoads.Add(bodyLoadElementZ);
                }
            }
            return new Tuple<Model, IModelReader>(model, modelReader);
        }

        private static IVectorView SolveModelsWithNewmark()
        {
            const double timestep = 1;
            const double time = 30;
            double[] muLame = new double[] { 6e4, 2.1e4 };
            double[] poissonV = new double[] { .45, .2 };
            IDynamicMaterial[] dynamicMaterials = new DynamicMaterial[] { new DynamicMaterial(1, 0, 0, true), new DynamicMaterial(1, 0, 0, true) };
            structModel = CreateStructuralModel(muLame, poissonV, dynamicMaterials, 0, new double[] { 0, 0, 0 }, lgElement);//.Item1; // new Model();
            var structuralModel = structModel.Item1;
            var structuralSolver = structuralBuilder.BuildSolver(structuralModel);
            var structuralProvider = new ProblemStructural(structuralModel, structuralSolver);
            //var structuralChildAnalyzer = new LinearAnalyzer(structuralModel, structuralSolver, structuralProvider);
            var increments = 2;
            var structuralChildAnalyzerBuilder = new LoadControlAnalyzer.Builder(structuralModel, structuralSolver, structuralProvider, increments);
            structuralChildAnalyzerBuilder.ResidualTolerance = 1E-5;
            structuralChildAnalyzerBuilder.MaxIterationsPerIncrement = 50;
            structuralChildAnalyzerBuilder.NumIterationsForMatrixRebuild = 2;
            //childAnalyzerBuilder.SubdomainUpdaters = new[] { new NonLinearSubdomainUpdater(model.SubdomainsDictionary[subdomainID]) }; // This is the default
            LoadControlAnalyzer structuralChildAnalyzer = structuralChildAnalyzerBuilder.Build();
            var structuralParentAnalyzer = new NewmarkDynamicAnalyzer(UpdateNewmarkModel, structuralModel, structuralSolver,
                structuralProvider, structuralChildAnalyzer, timestep, time, 0.25, 0.5);
            structuralParentAnalyzer.Initialize();

            for (int i = 0; i < time / timestep; i++)
            {
                if (lgElement == null)
                {
                    lgElement = new double[structModel.Item1.Elements.Count];
                    for (int j = 0; j < lgElement.Length; j++)
                    {
                        lgElement[j] = 1d;
                    }
                }
                for (int eID = 0; eID < structModel.Item2.elementDomains[0].Count; eID++)
                {
                    lgElement[eID] = lg[i];
                }
                structuralParentAnalyzer.SolveTimestep(i);
                Paraview(i);
            }

            return structuralSolver.LinearSystems[subdomainID].Solution;
        }
    }
}