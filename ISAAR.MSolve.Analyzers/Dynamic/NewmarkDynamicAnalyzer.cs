﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

//TODO: Optimization: I could avoid initialization and GC of some vectors by reusing existing ones.
//TODO: Use a base class for implicit time integration methods (perhaps to together with explicit)
namespace ISAAR.MSolve.Analyzers.Dynamic
{
    public class NewmarkDynamicAnalyzer : INonLinearParentAnalyzer //TODO: why is this non linear
    {
        private readonly double beta, gamma, timeStep, totalTime;
        private readonly double a0, a1, a2, a3, a4, a5, a6, a7;
        private IStructuralModel model;
        private IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private ISolver solver;
        private IImplicitIntegrationProvider provider;
        private Dictionary<int, IVector> rhs = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> uu = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> uum = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> uc = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> ucc = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> u = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> v = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> v1 = new Dictionary<int, IVector>();
        private Dictionary<int, IVector> v2 = new Dictionary<int, IVector>();
        private readonly Action<Dictionary<int, IVector>, Dictionary<int, IVector>, Dictionary<int, IVector>, IStructuralModel[], ISolver[], IImplicitIntegrationProvider[], IChildAnalyzer[]> CreateNewModel;
        IStructuralModel[] modelsForReplacement = new IStructuralModel[1];
        ISolver[] solversForReplacement = new ISolver[1];
        IImplicitIntegrationProvider[] providersForReplacement = new IImplicitIntegrationProvider[1];
        IChildAnalyzer[] childAnalyzersForReplacement = new IChildAnalyzer[1];

        public NewmarkDynamicAnalyzer(Action<Dictionary<int, IVector>, Dictionary<int, IVector>, Dictionary<int, IVector>, IStructuralModel[], ISolver[], IImplicitIntegrationProvider[], IChildAnalyzer[]> modelCreator, 
            IStructuralModel model, ISolver solver, IImplicitIntegrationProvider provider,
            IChildAnalyzer childAnalyzer, double timeStep, double totalTime, double alpha, double delta)
        {
            this.CreateNewModel = modelCreator;
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.provider = provider;
            this.ChildAnalyzer = childAnalyzer;
            this.beta = alpha;
            this.gamma = delta;
            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.ChildAnalyzer.ParentAnalyzer = this;

            modelsForReplacement[0] = model;
            solversForReplacement[0] = solver;
            providersForReplacement[0] = provider;
            childAnalyzersForReplacement[0] = childAnalyzer;

            // Initialize coefficients. It would make sense for them to be initialized in a different function, if they could 
            // change during the analysis
            a0 = 1 / (alpha * timeStep * timeStep);
            a1 = delta / (alpha * timeStep);
            a2 = 1 / (alpha * timeStep);
            a3 = 1 / (2 * alpha) - 1;
            a4 = delta / alpha - 1;
            a5 = timeStep * 0.5 * (delta / alpha - 2);
            a6 = timeStep * (1 - delta);
            a7 = delta * timeStep;
        }

        private NewmarkDynamicAnalyzer(IStructuralModel model, ISolver solver, IImplicitIntegrationProvider provider,
            IChildAnalyzer childAnalyzer, double timeStep, double totalTime, double alpha, double delta)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.provider = provider;
            this.ChildAnalyzer = childAnalyzer;
            this.beta = alpha;
            this.gamma = delta;
            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.ChildAnalyzer.ParentAnalyzer = this;

            // Initialize coefficients. It would make sense for them to be initialized in a different function, if they could 
            // change during the analysis
            a0 = 1 / (alpha * timeStep * timeStep);
            a1 = delta / (alpha * timeStep);
            a2 = 1 / (alpha * timeStep);
            a3 = 1 / (2 * alpha) - 1;
            a4 = delta / alpha - 1;
            a5 = timeStep * 0.5 * (delta / alpha - 2);
            a6 = timeStep * (1 - delta);
            a7 = delta * timeStep;
        }

        public Dictionary<int, IAnalyzerLog[]> Logs => null; //TODO: this can't be right
        public Dictionary<int, ImplicitIntegrationAnalyzerLog> ResultStorages { get; }
            = new Dictionary<int, ImplicitIntegrationAnalyzerLog>();

        public IChildAnalyzer ChildAnalyzer { get; set; }

        public void BuildMatrices()
        {
            var coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = a0,
                Damping = a1,
                Stiffness = 1
            };
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.Matrix = provider.LinearCombinationOfMatricesIntoStiffness(coeffs, linearSystem.Subdomain);
            }
        }

        public IVector GetOtherRhsComponents(ILinearSystem linearSystem, IVector currentSolution)
        {
            #region old code
            //// u[id]: old solution
            //// v[id]: current solution
            //// vv: old acceleration
            //// v2: current acceleration
            //// v1: current velocity
            ////double vv = v2[id].Data[j];
            ////v2[id].Data[j] = a0 * (v[id].Data[j] - u[id].Data[j]) - a2 * v1[id].Data[j] - a3 * vv;
            ////v1[id].Data[j] += a6 * vv + a7 * v2[id].Data[j];

            //int id = subdomain.ID;
            //Vector<double> currentAcceleration = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> currentVelocity = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uu = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uc = new Vector<double>(subdomain.Solution.Length);
            //for (int j = 0; j < subdomain.Rhs.Length; j++)
            //{
            //    currentAcceleration.Data[j] = a0 * (currentSolution[j] - v[id].Data[j]) - a2 * v1[id].Data[j] - a3 * v2[id].Data[j];
            //    currentVelocity.Data[j] = v1[id].Data[j] + a6 * v2[id].Data[j] + a7 * currentAcceleration.Data[j];
            //    uu.Data[j] = a0 * currentSolution[j] + a2 * currentVelocity.Data[j] + a3 * currentAcceleration.Data[j];
            //    uc.Data[j] = a1 * currentSolution[j] + a4 * currentVelocity.Data[j] + a5 * currentAcceleration.Data[j];
            //}

            //Vector<double> tempResult = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> result = new Vector<double>(subdomain.Solution.Length);
            //provider.MassMatrixVectorProduct(subdomain, uu, tempResult.Data);
            //result.Add(tempResult);

            //provider.DampingMatrixVectorProduct(subdomain, uc, tempResult.Data);
            //result.Add(tempResult);

            //return result.Data;

            //Vector<double> uu = new Vector<double>(subdomain.Solution.Length);
            //Vector<double> uc = new Vector<double>(subdomain.Solution.Length);
            //int id = subdomain.ID;
            //for (int j = 0; j < subdomain.Rhs.Length; j++)
            //{
            //    uu.Data[j] = -a0 * (v[id].Data[j] - currentSolution[j]) - a2 * v1[id].Data[j] - a3 * v2[id].Data[j];
            //    uc.Data[j] = -a1 * (v[id].Data[j] - currentSolution[j]) - a4 * v1[id].Data[j] - a5 * v2[id].Data[j];
            //}
            //provider.MassMatrixVectorProduct(subdomain, uu, tempResult.Data);
            //result.Add(tempResult);
            //provider.DampingMatrixVectorProduct(subdomain, uc, tempResult.Data);
            //result.Add(tempResult);

            ////CalculateRhsImplicit(subdomain, result.Data, false);
            ////result.Scale(-1d);
            #endregion

            // result = a0 * (M * u) + a1 * (C * u) 
            IVector result = provider.MassMatrixVectorProduct(linearSystem.Subdomain, currentSolution);
            IVector temp = provider.DampingMatrixVectorProduct(linearSystem.Subdomain, currentSolution);
            result.LinearCombinationIntoThis(a0, temp, a1);
            return result;
        }

        public void Initialize(bool isFirstAnalysis = true)
        {
            if (isFirstAnalysis)
            {
                // The order in which the next initializations happen is very important.
                model.ConnectDataStructures();
                solver.OrderDofs(false);
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.Reset(); // Necessary to define the linear system's size 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
            }
            else
            {
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    //TODO: Perhaps these shouldn't be done if an analysis has already been executed. The model will not be 
                    //      modified. Why should the linear system be?
                    linearSystem.Reset();
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
            }

            //TODO: Perhaps this should be called by the child analyzer
            BuildMatrices();

            // Loads must be created after building the matrices.
            //TODO: Some loads may not have to be recalculated each time the stiffness changes.
            model.AssignLoads(solver.DistributeNodalLoads);
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }

            //InitializeCoefficients();
            InitializeInternalVectors();
            //InitializeMatrices();
            InitializeRhs();

            if (ChildAnalyzer == null) throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
            ChildAnalyzer.Initialize(isFirstAnalysis);
        }

        public void InitializeInternal(bool isFirstAnalysis = true)
        {
            if (isFirstAnalysis)
            {
                // The order in which the next initializations happen is very important.
                model.ConnectDataStructures();
                solver.OrderDofs(false);
                //provider.ClearMatrices();
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    linearSystem.Reset(); // Necessary to define the linear system's size 
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
            }
            else
            {
                foreach (ILinearSystem linearSystem in linearSystems.Values)
                {
                    //TODO: Perhaps these shouldn't be done if an analysis has already been executed. The model will not be 
                    //      modified. Why should the linear system be?
                    linearSystem.Reset();
                    linearSystem.Subdomain.Forces = Vector.CreateZero(linearSystem.Size);
                }
            }

            //TODO: Perhaps this should be called by the child analyzer
            BuildMatrices();

            // Loads must be created after building the matrices.
            //TODO: Some loads may not have to be recalculated each time the stiffness changes.
            model.AssignLoads(solver.DistributeNodalLoads);
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector = linearSystem.Subdomain.Forces;
            }

            //InitializeCoefficients();
            //InitializeMatrices();
            InitializeRhs();

            if (ChildAnalyzer == null) throw new InvalidOperationException("Newmark analyzer must contain an embedded analyzer.");
            ChildAnalyzer.Initialize(isFirstAnalysis);
        }
        public void SolveTimestep(int i)
        {
            Debug.WriteLine("Newmark step: {0}", i);

            if (CreateNewModel != null)
            {
                CreateNewModel(uu, uc, v, modelsForReplacement, solversForReplacement, providersForReplacement, childAnalyzersForReplacement);
                model = modelsForReplacement[0];
                solver = solversForReplacement[0];
                linearSystems = solver.LinearSystems;
                provider = providersForReplacement[0];
                ChildAnalyzer = childAnalyzersForReplacement[0];
                ChildAnalyzer.ParentAnalyzer = this;

                InitializeInternal(true);
            }
            IDictionary<int, IVector> rhsVectors = provider.GetRhsFromHistoryLoad(i);
            foreach (var l in linearSystems.Values) l.RhsVector = rhsVectors[l.Subdomain.ID];
            InitializeRhs();
            CalculateRhsImplicit();

            DateTime start = DateTime.Now;
            ChildAnalyzer.Solve();
            DateTime end = DateTime.Now;

            UpdateVelocityAndAcceleration(i);
            UpdateResultStorages(start, end);
            // Print output results in *.txt file
            using (var fileName = new StreamWriter(@"C:\Users\Ody\Documents\Marie Curie\comsolModels\MSolveHyperelasticDynamicsImplicitResults.txt", true))
            {
                double currentTime = ((i + 1) * timeStep);
                string strTimeStep = currentTime.ToString();
                //var totalSolution = ChildAnalyzer.Responses[0][0];
                var totalSolution = ChildAnalyzer.Responses[0][16];
                string strTotalSolution = totalSolution.ToString();
                fileName.WriteLine(strTimeStep + ", " + strTotalSolution);
            }
            //if (CreateNewModel != null)
            //{
            //    CreateNewModel(uu, uc, v, modelsForReplacement, solversForReplacement, providersForReplacement, childAnalyzersForReplacement);
            //    model = modelsForReplacement[0];
            //    solver = solversForReplacement[0];
            //    linearSystems = solver.LinearSystems;
            //    provider = providersForReplacement[0];
            //    ChildAnalyzer = childAnalyzersForReplacement[0];
            //    ChildAnalyzer.ParentAnalyzer = this;

            //    InitializeInternal(true);
            //}
        }

        public void Solve()
        {
            int numTimeSteps = (int)(totalTime / timeStep);
            for (int i = 0; i < numTimeSteps; ++i)
            {
                SolveTimestep(i);
            }
        }

        private void CalculateRhsImplicit()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.RhsVector = CalculateRhsImplicit(linearSystem, true);

                #region old code
                //int id = subdomain.ID;
                //for (int i = 0; i < subdomain.Rhs.Length; i++)
                //{
                //    //uu[id].Data[i] = v[id].Data[i] + a2a0 * v1[id].Data[i] + a3a0 * v2[id].Data[i];
                //    uu[id].Data[i] = a0 * v[id].Data[i] + a2 * v1[id].Data[i] + a3 * v2[id].Data[i];
                //    uc[id].Data[i] = a1 * v[id].Data[i] + a4 * v1[id].Data[i] + a5 * v2[id].Data[i];
                //}
                ////provider.Ms[id].Multiply(uu[id], uum[id].Data);
                //provider.MassMatrixVectorProduct(subdomain, uu[id], uum[id].Data);
                //provider.DampingMatrixVectorProduct(subdomain, uc[id], ucc[id].Data);
                //for (int i = 0; i < subdomain.Rhs.Length; i++)
                //    subdomain.Rhs[i] = rhs[id].Data[i] + uum[id].Data[i] + ucc[id].Data[i];
                #endregion
            }
            #region Old Fortran code
            //If (ptAnal.tDynamic.bHasDamping) Then
            //    atDynUU(iMesh).afValues = atDynV(iMesh).afValues + fA2A0 * atDynV1(iMesh).afValues + fA3A0 * atDynV2(iMesh).afValues
            //    atDynUC(iMesh).afValues = atDynV(iMesh).afValues + fA4A1 * atDynV1(iMesh).afValues + fA5A1 * atDynV2(iMesh).afValues
            //Else
            //    atDynUU(iMesh).afValues = atDynV(iMesh).afValues + fA2A0 * atDynV1(iMesh).afValues + fA3A0 * atDynV2(iMesh).afValues
            //End If
            //If (Present(ProcessVel)) Then
            //    Call ProcessVel(iMesh, Size(atDynUU(iMesh).afValues), atDynUU(iMesh).afValues, atDynUUM(iMesh).afValues)
            //Else
            //    Call MatrixVectorMultiplication(iMesh, Size(atRhs(iMesh).afValues), 2, atDynUU(iMesh).afValues, atDynUUM(iMesh).afValues)
            //End If
            //If (ptAnal.tDynamic.bHasDamping) Then
            //    If (Present(ProcessAcc)) Then
            //        Call ProcessAcc(iMesh, Size(atDynUU(iMesh).afValues), atDynUC(iMesh).afValues, atDynUCC(iMesh).afValues)
            //    Else
            //        Call MatrixVectorMultiplication(iMesh, Size(atRhs(iMesh).afValues), 3, atDynUC(iMesh).afValues, atDynUCC(iMesh).afValues)
            //    End If
            //    atRhs(iMesh).afValues = atRhs(iMesh).afValues + atDynUUM(iMesh).afValues + atDynUCC(iMesh).afValues
            //Else
            //    atRhs(iMesh).afValues = atRhs(iMesh).afValues + atDynUUM(iMesh).afValues
            //End If
            #endregion
        }

        private IVector CalculateRhsImplicit(ILinearSystem linearSystem, bool addRhs)
        {
            //TODO: instead of creating a new Vector and then trying to set ILinearSystem.RhsVector, clear it and operate on it.
            int id = linearSystem.Subdomain.ID;

            // uu = a0 * v + a2 * v1 + a3 * v2
            uu[id] = v[id].LinearCombination(a0, v1[id], a2);
            uu[id].AxpyIntoThis(v2[id], a3);

            // uc = a1 * v + a4 * v1 + a5 * v2
            uc[id] = v[id].LinearCombination(a1, v1[id], a4);
            uc[id].AxpyIntoThis(v2[id], a5);

            uum[id] = provider.MassMatrixVectorProduct(linearSystem.Subdomain, uu[id]);
            ucc[id] = provider.DampingMatrixVectorProduct(linearSystem.Subdomain, uc[id]);

            IVector rhsResult = uum[id].Add(ucc[id]);
            if (addRhs) rhsResult.AddIntoThis(rhs[id]);
            return rhsResult;
        }

        //private void InitializeCoefficients()
        //{
        //    if (delta < 0.5) throw new ArgumentException("Newmark delta has to be bigger than 0.5.");
        //    double aLimit = 0.25 * Math.Pow(0.5 + delta, 2);
        //    if (alpha < aLimit) throw new ArgumentException($"Newmark alpha has to be bigger than {aLimit}.");

        //    a0 = 1 / (alpha * timeStep * timeStep);
        //    a1 = delta / (alpha * timeStep);
        //    a2 = 1 / (alpha * timeStep);
        //    a3 = 1 / (2 * alpha) - 1;
        //    a4 = delta / alpha - 1;
        //    a5 = timeStep * 0.5 * (delta / alpha - 2);
        //    a6 = timeStep * (1 - delta);
        //    a7 = delta * timeStep;
        //}

        private void InitializeInternalVectors()
        {
            uu.Clear();
            uum.Clear();
            uc.Clear();
            ucc.Clear();
            u.Clear();
            v.Clear();
            v1.Clear();
            v2.Clear();
            rhs.Clear();

            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                uu.Add(id, linearSystem.CreateZeroVector());
                uum.Add(id, linearSystem.CreateZeroVector());
                uc.Add(id, linearSystem.CreateZeroVector());
                ucc.Add(id, linearSystem.CreateZeroVector());
                u.Add(id, linearSystem.CreateZeroVector());
                //v.Add(id, linearSystem.CreateZeroVector());
                v1.Add(id, linearSystem.CreateZeroVector());
                v2.Add(id, linearSystem.CreateZeroVector());
                rhs.Add(id, linearSystem.CreateZeroVector());

                // Account for initial conditions coming from a previous solution. 
                //TODO: This doesn't work as intended. The solver (previously the LinearSystem) initializes the solution to zero.
                if (linearSystem.Solution != null) v[id] = linearSystem.Solution.Copy();
                else v.Add(id, linearSystem.CreateZeroVector());
            }
        }

        private void InitializeRhs()
        {
            ImplicitIntegrationCoefficients coeffs = new ImplicitIntegrationCoefficients
            {
                Mass = a0, Damping = a1, Stiffness = 1
            };
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                provider.ProcessRhs(coeffs, linearSystem.Subdomain, linearSystem.RhsVector);
                rhs[linearSystem.Subdomain.ID] = linearSystem.RhsVector.Copy(); //TODO: copying the vectors is wasteful.
            }
        }

        private void UpdateResultStorages(DateTime start, DateTime end)
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                if (ResultStorages.ContainsKey(id))
                    if (ResultStorages[id] != null)
                        foreach (var l in ChildAnalyzer.Logs[id])
                            ResultStorages[id].StoreResults(start, end, l);
            }
        }

        private void UpdateVelocityAndAcceleration(int timeStep)
        {
            var externalVelocities = provider.GetVelocitiesOfTimeStep(timeStep);
            var externalAccelerations = provider.GetAccelerationsOfTimeStep(timeStep);

            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                int id = linearSystem.Subdomain.ID;
                u[id].CopyFrom(v[id]); //TODO: this copy can be avoided by pointing to v[id] and then v[id] = null;
                //v[id].CopyFrom(linearSystem.Solution);
                v[id].CopyFrom(ChildAnalyzer.Responses[id]);
                IVector vv = v2[id].Add(externalAccelerations[id]);

                // v2 = a0 * (v - u) - a2 * v1 - a3 * vv
                v2[id] = v[id].Subtract(u[id]);
                v2[id].LinearCombinationIntoThis(a0, v1[id], -a2);
                v2[id].AxpyIntoThis(vv, -a3);

                // v1 = v1 + externalVelocities + a6 * vv + a7 * v2
                v1[id].AddIntoThis(externalVelocities[id]);
                v1[id].AxpyIntoThis(vv, a6);
                v1[id].AxpyIntoThis(v2[id], a7);

                if ((timeStep + 1) % 1 == 0)
                {
                    string path0 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                    //string path1 = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput\temperature0.txt";
                    //string path = @"C:\Users\Ody\Documents\Marie Curie\comsolModels\MsolveOutput";
                    var path2 = Path.Combine(path0, $"displacement{timeStep}.txt");
                    var writer = new LinearAlgebra.Output.FullVectorWriter() { ArrayFormat = Array1DFormat.PlainVertical };
                    writer.WriteToFile(v[id], path2);
                    //writer.WriteToFile(temperature[id][0], path1);

                    //File.AppendAllLines(path1, new string[] { temperature[id][0].ToString() }, Encoding.UTF8);
                    //File.AppendAllLines(path2, new string[] { temperature[id][340].ToString() }, Encoding.UTF8);
                }

            }
        }

        public class Builder
        {
            private readonly double timeStep, totalTime; //TODO: perhaps totalTime should be numTimeSteps
            private readonly IChildAnalyzer childAnalyzer;
            private readonly IStructuralModel model;
            private readonly ISolver solver;
            private readonly IImplicitIntegrationProvider provider;
            private double beta = 0.25, gamma = 0.5; // constant acceleration is the default

            public Builder(IStructuralModel model, ISolver solver, IImplicitIntegrationProvider provider,
                IChildAnalyzer childAnalyzer, double timeStep, double totalTime)
            {
                this.model = model;
                this.solver = solver;
                this.provider = provider;
                this.childAnalyzer = childAnalyzer;

                this.timeStep = timeStep;
                this.totalTime = totalTime;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="beta">
            /// Used in the intepolation between the accelerations of the previous and current time step, in order to obtain the 
            /// current displacements. Also called alpha by Bathe.
            /// </param>
            /// <param name="gamma">
            /// Used in the intepolation between the accelerations of the previous and current time step, in order to obtain the 
            /// current velocities. Also called delta by Bathe.
            /// </param>
            /// <param name="allowConditionallyStable">
            /// If set to true, the user must make sure that the time step chosen is lower than the critical step size 
            /// corresponding to these particular <paramref name="beta"/>, <paramref name="gamma"/> parameters.
            /// </param>
            public void SetNewmarkParameters(double beta, double gamma, bool allowConditionallyStable = false)
            {
                if (!allowConditionallyStable)
                {
                    if (gamma < 0.5) throw new ArgumentException(
                        "Newmark delta has to be bigger than 0.5 to ensure unconditional stability.");
                    if (beta < 0.25) throw new ArgumentException(
                        "Newmark alpha has to be bigger than 0.25 to ensure unconditional stability.");

                    // No idea where Bathe got this from.
                    //double bLimit = 0.25 * Math.Pow(0.5 + delta, 2);
                    //if (beta < bLimit) throw new ArgumentException(
                    //$"Newmark beta has to be bigger than {bLimit} to ensure unconditional stability.");
                }
                if (gamma < 0.5) throw new ArgumentException("Newmark delta has to be bigger than 0.5.");
                double aLimit = 0.25 * Math.Pow(0.5 + gamma, 2);
                if (beta < aLimit) throw new ArgumentException($"Newmark alpha has to be bigger than {aLimit}.");

                this.gamma = gamma;
                this.beta = beta;
            }

            /// <summary>
            /// Central diffences: gamma = 1/2, beta = 0. Newmark results in central diffences, a conditionally stable explicit 
            /// method. To ensure stability, the time step must be &lt;= the critical step size = 2 / w,  where w is the maximum 
            /// natural radian frequency. It would be more efficient to use an explicit dynamic analyzer. 
            /// </summary>
            public void SetNewmarkParametersForCentralDifferences()
            {
                gamma = 0.5;
                beta = 0.0;
            }

            /// <summary>
            /// Constant acceleration (also called average acceleration or trapezoid rule): gamma = 1/2, beta = 1/4. 
            /// This is the most common scheme and is unconditionally stable. In this analyzer, it is used as the default.
            /// </summary>
            public void SetNewmarkParametersForConstantAcceleration()
            {
                gamma = 0.5;
                beta = 0.25;
            }

            /// <summary>
            /// Linear acceleration: gamma = 1/2, beta = 1/6. This is more accurate than the default constant acceleration, 
            /// but it conditionally stable. To ensure stability, the time step must be &lt;= the critical step size = 3.464 / w 
            /// = 0.551 * T, where w is the maximum natural radian frequency and T is the minimum natural period.
            /// </summary>
            public void SetNewmarkParametersForLinearAcceleration()
            {
                gamma = 0.5;
                beta = 1.0 / 6.0;
            }

            public NewmarkDynamicAnalyzer Build()
                => new NewmarkDynamicAnalyzer(model, solver, provider, childAnalyzer, timeStep, totalTime, beta, gamma);
        }
    }
}
