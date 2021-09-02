using System;
using System.Collections.Generic;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.Analyzers.NonLinear;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Logging.Interfaces;
using ISAAR.MSolve.Solvers;
using ISAAR.MSolve.Solvers.LinearSystems;

namespace ISAAR.MSolve.Analyzers
{
    public class StaticAnalyzer : INonLinearParentAnalyzer
    {
        public  IReadOnlyDictionary<int, ILinearSystem> linearSystems;
        private IStructuralModel model;
        private IStaticProvider provider;
        private ISolver solver;
        private readonly Action<IStructuralModel[], ISolver[], IStaticProvider[], IChildAnalyzer[]> CreateNewModel;
        private readonly Action<IChildAnalyzer[]> UpdateSolution;
        IStructuralModel[] modelsForReplacement = new IStructuralModel[1];
        ISolver[] solversForReplacement = new ISolver[1];
        IStaticProvider[] providersForReplacement = new IStaticProvider[1];
        IChildAnalyzer[] childAnalyzersForReplacement = new IChildAnalyzer[1];

        public StaticAnalyzer(Action<IStructuralModel[], ISolver[], IStaticProvider[], IChildAnalyzer[]> modelCreator,
            Action<IChildAnalyzer[]> solutionUpdater, IStructuralModel model, ISolver solver, IStaticProvider provider,
            IChildAnalyzer childAnalyzer)
        {
            this.CreateNewModel = modelCreator;
            this.UpdateSolution = solutionUpdater;
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.provider = provider;
            this.ChildAnalyzer = childAnalyzer;
            this.ChildAnalyzer.ParentAnalyzer = this;

            modelsForReplacement[0] = model;
            solversForReplacement[0] = solver;
            providersForReplacement[0] = provider;
            childAnalyzersForReplacement[0] = childAnalyzer;
        }
        public StaticAnalyzer(IStructuralModel model, ISolver solver, IStaticProvider provider,
            IChildAnalyzer childAnalyzer)
        {
            this.model = model;
            this.linearSystems = solver.LinearSystems;
            this.solver = solver;
            this.provider = provider;
            this.ChildAnalyzer = childAnalyzer;
            this.ChildAnalyzer.ParentAnalyzer = this;
        }

        public Dictionary<int, IAnalyzerLog[]> Logs { get; } = new Dictionary<int, IAnalyzerLog[]>();

        public IChildAnalyzer ChildAnalyzer { get; set; }

        public void BuildMatrices()
        {
            foreach (ILinearSystem linearSystem in linearSystems.Values)
            {
                linearSystem.Matrix = provider.CalculateMatrix(linearSystem.Subdomain);
            }
        }

        public IVector GetOtherRhsComponents(ILinearSystem linearSystem, IVector currentSolution)
        {
            //TODO: use a ZeroVector class that avoid doing useless operations or refactor this method. E.g. let this method 
            // alter the child analyzer's rhs vector, instead of the opposite (which is currently done).
            return linearSystem.CreateZeroVector();
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

            if (ChildAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            ChildAnalyzer.Initialize(isFirstAnalysis);
        }

        public void Solve()
        {
            if (CreateNewModel != null)
            {
                CreateNewModel(modelsForReplacement, solversForReplacement, providersForReplacement, childAnalyzersForReplacement);
                model = modelsForReplacement[0];
                solver = solversForReplacement[0];
                linearSystems = solver.LinearSystems;
                provider = providersForReplacement[0];
                ChildAnalyzer = childAnalyzersForReplacement[0];
                ChildAnalyzer.ParentAnalyzer = this;

                Initialize(true);
            }
            if (ChildAnalyzer == null) throw new InvalidOperationException("Static analyzer must contain an embedded analyzer.");
            ChildAnalyzer.Solve();
            if (UpdateSolution != null)
            {
                UpdateSolution(childAnalyzersForReplacement);
            }
        }
    }
}
