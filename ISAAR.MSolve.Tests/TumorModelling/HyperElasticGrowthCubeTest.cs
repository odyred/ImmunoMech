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
    public class HyperElasticGrowthCubeTest
    {
        private const int subdomainID = 0;
        private static double[] lambdag;
        private static double[] currentLambdag;
        private const double a = 1;
        private const double b = 1;
        private static DenseMatrixSolver.Builder builder = new DenseMatrixSolver.Builder();
        private static SkylineSolver.Builder structuralBuilder = new SkylineSolver.Builder();
        private static Dictionary<int, IVector>[] Solutions;
        private static Dictionary<int, IVector> Accelerations;
        private static Dictionary<int, IVector> Velocities;
        private static Dictionary<int, IVector> Displacements;
        private static Tuple<Model, IModelReader> gModel, structModel;
        private static string inputFile = "mesh446elem.mphtxt";
        private static int NewtonRaphsonIncrements = 2;
        private static int NewtonRaphosnIterations = 50;
        private static double NewtonRaphsonTolerarance = 1e-5;
        private static int NewtonRaphsonIterForMatrixRebuild = 10;
        private static double MultiModelAnalyzerTolerance = 5e-3;
        private const double timestep = 1;
        private const double time = 30;
        private static double[][] nodalCoordinates;

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
            IVectorView solution = SolveModelsWithNewmark();

            Assert.True(CompareResults(solution));
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
                    outputFile.WriteLine($"{nodalCoordinates[i][0]} {nodalCoordinates[i][1]} {nodalCoordinates[i][2]} ");
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
                    double dist = Math.Sqrt(Math.Pow(nodalCoordinates[i][0] - structModel.Item1.Nodes[i].X, 2) +
                                            Math.Pow(nodalCoordinates[i][1] - structModel.Item1.Nodes[i].Y, 2) +
                                            Math.Pow(nodalCoordinates[i][2] - structModel.Item1.Nodes[i].Z, 2));
                    outputFile.WriteLine($"{dist} ");
                }
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
        private static void UpdateNodes(Model structuralModel, Dictionary<int, IVector> displacements = null)
        {//only works for tet4 now
            //this.displacements = displacements;
            if (displacements != null)
            {
                for (int id = 0; id < displacements.Count; id++)
                {
                    int count = 0;
                    foreach (Node n in structuralModel.Nodes)
                    {
                        nodalCoordinates[n.ID] = new double[3];
                        double x = n.X;
                        double y = n.Y;
                        double z = n.Z;

                        if (structuralModel.GlobalDofOrdering.GlobalFreeDofs.Contains(n, StructuralDof.TranslationX))
                        {
                            x += displacements[id][count];
                            count += 1;
                        }
                        if (structuralModel.GlobalDofOrdering.GlobalFreeDofs.Contains(n, StructuralDof.TranslationY))
                        {
                            y += displacements[id][count];
                            count += 1;
                        }
                        if (structuralModel.GlobalDofOrdering.GlobalFreeDofs.Contains(n, StructuralDof.TranslationZ))
                        {
                            z += displacements[id][count];
                            count += 1;
                        }
                        nodalCoordinates[n.ID][0] = x;
                        nodalCoordinates[n.ID][1] = y;
                        nodalCoordinates[n.ID][2] = z;
                    }
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

        private static void UpdateNewmarkModel(Dictionary<int, IVector> accelerations, Dictionary<int, IVector> velocities, Dictionary<int, IVector> displacements, IStructuralModel[] modelsToReplace,
            ISolver[] solversToReplace, IImplicitIntegrationProvider[] providersToReplace, IChildAnalyzer[] childAnalyzersToReplace)
        {
            Displacements = displacements;
            ReplaceLambdaGInModel(modelsToReplace[0], currentLambdag);
            solversToReplace[0] = structuralBuilder.BuildSolver(modelsToReplace[0]);
            providersToReplace[0] = new ProblemStructural(modelsToReplace[0], solversToReplace[0]);
            solversToReplace[0].HandleMatrixWillBeSet();
            var increments = NewtonRaphsonIncrements;
            var childAnalyzerBuilder = new LoadControlAnalyzer.Builder(modelsToReplace[0], solversToReplace[0], (INonLinearProvider)providersToReplace[0], increments);
            childAnalyzerBuilder.ResidualTolerance = NewtonRaphsonTolerarance;
            childAnalyzerBuilder.MaxIterationsPerIncrement = NewtonRaphosnIterations;
            childAnalyzerBuilder.NumIterationsForMatrixRebuild = NewtonRaphsonIterForMatrixRebuild;
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
            string[] lgString = System.IO.File.ReadAllLines(Path.Combine(Directory.GetCurrentDirectory(), "InputFiles", "TumorGrowthModel", "lg.txt"));

            lambdag = new double[lgString.Length];
            for (int i = 0; i < lgString.Length; i++)
            {
                lambdag[i] = Double.Parse(lgString[i], CultureInfo.InvariantCulture);
            }
        }
        private static Tuple<Model, IModelReader> CreateStructuralModel(double[] MuLame, double[] PoissonV, IDynamicMaterial[] commonDynamicMaterialProperties,
            double b, double[] l)
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
            if (currentLambdag == null)
            {
                modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties);
            }
            else
            {
                modelReader = new ComsolMeshReader1(filename, C1, C2, bulkModulus, commonDynamicMaterialProperties, currentLambdag);
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
            return new Tuple<Model, IModelReader>(model, modelReader);
        }
        private static IVectorView SolveModelsWithNewmark()
        {
            double[] muLame = new double[] { 6e4, 2.1e4 };
            double[] poissonV = new double[] { .45, .2 };
            IDynamicMaterial[] dynamicMaterials = new DynamicMaterial[] { new DynamicMaterial(.001, 0, 0, true), new DynamicMaterial(.001, 0, 0, true) };
            structModel = CreateStructuralModel(muLame, poissonV, dynamicMaterials, 0, new double[] { 0, 0, 0 });//.Item1; // new Model();
            currentLambdag = new double[structModel.Item1.Elements.Count];
            foreach (Element e in structModel.Item2.elementDomains[1])
            {
                currentLambdag[e.ID] = 1;
            }
            CreateGrowthModel();
            nodalCoordinates = new double[structModel.Item1.Nodes.Count][];
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
            var structuralParentAnalyzer = new NewmarkDynamicAnalyzer(UpdateNewmarkModel, structuralModel, structuralSolver,
                structuralProvider, structuralChildAnalyzer, timestep, time, 0.25, 0.5);
            structuralParentAnalyzer.Initialize();

            for (int i = 0; i < time / timestep; i++)
            {
                foreach (Element e in structModel.Item2.elementDomains[0])
                {
                    currentLambdag[e.ID] = lambdag[i];
                }
                structuralParentAnalyzer.SolveTimestep(i);
                UpdateNodes(structModel.Item1, Displacements);
                Paraview(i);
            }

            return structuralSolver.LinearSystems[0].Solution;
        }
    }
}