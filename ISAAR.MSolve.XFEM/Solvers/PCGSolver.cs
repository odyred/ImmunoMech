﻿using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.LinearSystems;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Algorithms;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Preconditioning;
using ISAAR.MSolve.LinearAlgebra.LinearSystems.Statistics;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Assemblers;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Solvers
{
    class PCGSolver: SolverBase
    {
        private readonly double maxIterationsOverOrder; //roughly
        private readonly double tolerance;

        public PCGSolver(Model2D model, double maxIterationsOverOrder, double tolerance) : base(model)
        {
            this.maxIterationsOverOrder = maxIterationsOverOrder;
            this.tolerance = tolerance;
        }

        public override void Solve()
        {
            var watch = new Stopwatch();
            watch.Start();

            // Interleaced and separate dof enumerators seem to have similar performance.
            DofOrderer = InterleavedDofOrderer.Create(model);
            //DofOrderer = DofOrdererSeparate.Create(model);

            var assembler = new GlobalCSRAssembler();
            (DOKRowMajor Kuu, DOKRowMajor Kuc) = assembler.BuildGlobalMatrix(model, DofOrderer);
            Vector rhs = CalcEffectiveRhs(Kuc);


            int maxIterations = (int)Math.Ceiling(Kuu.NumColumns * maxIterationsOverOrder);

            // PCG
            //var pcg = new PreconditionedConjugateGradient(maxIterations, tolerance);
            //(double[] diagonal, int firstZeroIdx) = Kuu.GetDiagonalAsArray(); // Preconditioner could be abstracted, but I think it depends on the solver.
            //var preconditioner = new JacobiPreconditioner(diagonal);
            //var preconditioner = new IdentityPreconditioner(true);
            //(Vector x, IterativeStatistics statistics) = pcg.Solve(Kuu.BuildCSRMatrix(true), rhs, preconditioner);
            
            // CG
            var cg = new ConjugateGradient(maxIterations, tolerance);
            (Vector x, IterativeStatistics statistics) = cg.Solve(Kuu.BuildCSRMatrix(true), rhs);

            Console.WriteLine(statistics);
            Solution = x;

            watch.Stop();
            Logger.SolutionTimes.Add(watch.ElapsedMilliseconds);
        }
    }
}
