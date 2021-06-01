﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Solvers.Assemblers;
using ISAAR.MSolve.Solvers.Ordering;
using ISAAR.MSolve.Solvers.Ordering.Reordering;

namespace ISAAR.MSolve.Solvers.Direct
{
    /// <summary>
    /// Direct solver for models with only 1 subdomain. Uses Cholesky factorization on sparse symmetric positive definite 
    /// matrices stored in Skyline format.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineSolver : SingleSubdomainSolverBase<SkylineMatrix>
    {
        private readonly double factorizationPivotTolerance;

        private bool factorizeInPlace = true;
        private bool mustFactorize = true;
        private LdlSkyline factorizedMatrix;

        private SkylineSolver(IStructuralModel model, double factorizationPivotTolerance, IDofOrderer dofOrderer):
            base(model, dofOrderer, new SkylineAssembler(), "SkylineSolver")
        {
            this.factorizationPivotTolerance = factorizationPivotTolerance;
        }

        public override void HandleMatrixWillBeSet()
        {
            mustFactorize = true;
            factorizedMatrix = null;
        }

        public override void Initialize() { }

        public override void PreventFromOverwrittingSystemMatrices() => factorizeInPlace = false;

        /// <summary>
        /// Solves the linear system with back-forward substitution. If the matrix has been modified, it will be refactorized.
        /// </summary>
        public override void Solve()
        {
            var watch = new Stopwatch();
            if (linearSystem.SolutionConcrete == null) linearSystem.SolutionConcrete = linearSystem.CreateZeroVectorConcrete();
            //else linearSystem.Solution.Clear(); // no need to waste computational time on this in a direct solver

            // Factorization
            if (mustFactorize)
            {
                watch.Start();
                factorizedMatrix = linearSystem.Matrix.FactorLdl(factorizeInPlace, factorizationPivotTolerance);
                watch.Stop();
                Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);
                watch.Reset();
                mustFactorize = false;
            }

            // Substitutions
            watch.Start();
            factorizedMatrix.SolveLinearSystem(linearSystem.RhsConcrete, linearSystem.SolutionConcrete);
            watch.Stop();
            Logger.LogTaskDuration("Back/forward substitutions", watch.ElapsedMilliseconds);
            Logger.IncrementAnalysisStep();
        }

        protected override Matrix InverseSystemMatrixTimesOtherMatrix(IMatrixView otherMatrix)
        {
            var watch = new Stopwatch();

            // Factorization
            if (mustFactorize)
            {
                watch.Start();
                factorizedMatrix = linearSystem.Matrix.FactorLdl(factorizeInPlace, factorizationPivotTolerance);
                watch.Stop();
                Logger.LogTaskDuration("Matrix factorization", watch.ElapsedMilliseconds);
                watch.Reset();
                mustFactorize = false;
            }

            // Substitutions
            watch.Start();
            int systemOrder = linearSystem.Matrix.NumColumns;
            int numRhs = otherMatrix.NumColumns;
            var solutionVectors = Matrix.CreateZero(systemOrder, numRhs);
            if (otherMatrix is Matrix otherDense)
            {
                factorizedMatrix.SolveLinearSystems(otherDense, solutionVectors);
            }
            else
            {
                try
                {
                    // If there is enough memory, copy the RHS matrix to a dense one, to speed up computations. 
                    //TODO: must be benchmarked, if it is actually more efficient than solving column by column.
                    Matrix rhsVectors = otherMatrix.CopyToFullMatrix();
                    factorizedMatrix.SolveLinearSystems(rhsVectors, solutionVectors);
                }
                catch (InsufficientMemoryException) //TODO: what about OutOfMemoryException?
                {
                    // Solve each linear system separately, to avoid copying the RHS matrix to a dense one.
                    Vector solutionVector = linearSystem.CreateZeroVectorConcrete();
                    for (int j = 0; j < numRhs; ++j)
                    {
                        if (j != 0) solutionVector.Clear();
                        Vector rhsVector = otherMatrix.GetColumn(j);
                        factorizedMatrix.SolveLinearSystem(rhsVector, solutionVector);
                        solutionVectors.SetSubcolumn(j, solutionVector);
                    }
                }
            }
            watch.Stop();
            Logger.LogTaskDuration("Back/forward substitutions", watch.ElapsedMilliseconds);
            Logger.IncrementAnalysisStep();
            return solutionVectors;
        }

        public class Builder : ISolverBuilder
        {
            public Builder() { }

            public IDofOrderer DofOrderer { get; set; } 
                = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());

            public double FactorizationPivotTolerance { get; set; } = 1E-15;

            ISolver ISolverBuilder.BuildSolver(IStructuralModel model) => BuildSolver(model);

            public SkylineSolver BuildSolver(IStructuralModel model)
            {
                return new SkylineSolver(model, FactorizationPivotTolerance, DofOrderer);
            }
        }
    }
}
