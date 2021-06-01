﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging;
using ISAAR.MSolve.Solvers;

namespace ISAAR.MSolve.Analyzers.NonLinear
{
    public class LoadControlAnalyzer : NonLinearAnalyzerBase
    {
        private LoadControlAnalyzer(IStructuralModel model, ISolver solver, INonLinearProvider provider,
            IReadOnlyDictionary<int, INonLinearSubdomainUpdater> subdomainUpdaters,
            int numIncrements, int maxIterationsPerIncrement, int numIterationsForMatrixRebuild, double residualTolerance) : 
            base(model, solver, provider, subdomainUpdaters, numIncrements, maxIterationsPerIncrement, 
                numIterationsForMatrixRebuild, residualTolerance)
        { }
        
        public override void Solve()
        {
            InitializeLogs();

            DateTime start = DateTime.Now;
            UpdateInternalVectors();//TODOMaria this divides the externally applied load by the number of increments and scatters it to all subdomains and stores it in the class subdomain dictionary and total external load vector
            for (int increment = 0; increment < numIncrements; increment++)
            {
                double errorNorm = 0;
                ClearIncrementalSolutionVector();//TODOMaria this sets du to 0
                UpdateRhs(increment);//TODOMaria this copies the residuals stored in the class dictionary to the subdomains

                double firstError = 0;
                int iteration = 0;
                for (iteration = 0; iteration < maxIterationsPerIncrement; iteration++)
                {
                    if (iteration == maxIterationsPerIncrement - 1) return;
                    if (Double.IsNaN(errorNorm)) return;
                    solver.Solve();
                    //double rhsNormIt = solver.LinearSystems.First().Value.RhsVector.Norm2();
                    //double xNormIt = solver.LinearSystems.First().Value.Solution.Norm2();
                    Dictionary<int, IVector> internalRhsVectors = CalculateInternalRhs(increment, iteration);
                    double residualNormCurrent = UpdateResidualForcesAndNorm(increment, internalRhsVectors); // This also sets the rhs vectors in linear systems.
                    errorNorm = globalRhsNormInitial != 0 ? residualNormCurrent / globalRhsNormInitial : 0;// (rhsNorm*increment/increments) : 0;//TODOMaria this calculates the internal force vector and subtracts it from the external one (calculates the residual)
                    //Console.WriteLine($"Increment {increment}, iteration {iteration}: norm2(error) = {errorNorm}");

                    if (iteration == 0) firstError = errorNorm;

                    if (TotalDisplacementsPerIterationLog != null) TotalDisplacementsPerIterationLog.StoreDisplacements(uPlusdu);

                    if (errorNorm < residualTolerance)
                    {
                        foreach (var subdomainLogPair in IncrementalLogs)
                        {
                            int subdomainID = subdomainLogPair.Key;
                            TotalLoadsDisplacementsPerIncrementLog log = subdomainLogPair.Value;
                            log.LogTotalDataForIncrement(increment, iteration, errorNorm,
                                uPlusdu[subdomainID], internalRhsVectors[subdomainID]);
                        }
                        break;
                    }

                    SplitResidualForcesToSubdomains();//TODOMaria scatter residuals to subdomains
                    if ((iteration + 1) % numIterationsForMatrixRebuild == 0) // Matrix rebuilding should be handled in another way. E.g. in modified NR, it must be done at each increment.
                    {
                        provider.Reset();
                        BuildMatrices();
                    }
                }
                //double rhsNormInc = solver.LinearSystems.First().Value.RhsVector.Norm2();
                //double xNormInc = solver.LinearSystems.First().Value.Solution.Norm2();
                Debug.WriteLine("NR {0}, first error: {1}, exit error: {2}", iteration, firstError, errorNorm);
                SaveMaterialStateAndUpdateSolution();
            }
            //            ClearMaterialStresses();

            // TODO: Logging should be done at each iteration. And it should be done using pull observers
            DateTime end = DateTime.Now;
            StoreLogResults(start, end);
        }

        public class Builder: NonLinearAnalyzerBuilderBase
        {
            public Builder(IStructuralModel model, ISolver solver, INonLinearProvider provider, int numIncrements) :
                base(model, solver, provider, numIncrements)
            {
                MaxIterationsPerIncrement = 1000;
                NumIterationsForMatrixRebuild = 1;
                ResidualTolerance = 1E-3;
            }

            public LoadControlAnalyzer Build()
            {
                return new LoadControlAnalyzer(model, solver, provider, SubdomainUpdaters,
                    numIncrements, maxIterationsPerIncrement, numIterationsForMatrixRebuild, residualTolerance);
            }
        }
    }
}
