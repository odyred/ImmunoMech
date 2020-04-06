using ISAAR.MSolve.Analyzers;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System.Linq;

//TODO: time logging must be refactored
//TODO: perhaps this belongs to Solvers.Assemblers, since the vector type depends on the solver. In that case, the 
//      elementMatrixProvider should be injected by the problem/provider.
namespace ISAAR.MSolve.Problems
{
    /// <summary>
    /// Calculates the equivalent nodal forces (at the subdomain level) due to Dirichlet boundary conditions.
    /// Authors: Maria Tavlaki
    /// </summary>
    public class DirichletEquivalentLoadsConvectionDiffusion : IDirichletEquivalentLoadsAssembler
    {
        private IBoundaryElementMatrixProvider elementProvider; //TODO: not sure if df = K * du is the best way to calcuate df.

        public DirichletEquivalentLoadsConvectionDiffusion(IBoundaryElementMatrixProvider elementProvider)
        {
            this.elementProvider = elementProvider;
        }

        public IVector GetEquivalentNodalLoads(ISubdomain subdomain, IVectorView solution, double constraintScalingFactor) 
        {
            //var times = new Dictionary<string, TimeSpan>();
            //var totalStart = DateTime.Now;
            //times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            //times.Add("element", TimeSpan.Zero);
            //times.Add("addition", TimeSpan.Zero);

            var subdomainEquivalentForces = Vector.CreateZero(subdomain.FreeDofOrdering.NumFreeDofs);
            foreach (IConvectionDiffusionBoundaryElement element in subdomain.Elements.Where(x=>x is IConvectionDiffusionBoundaryElement)) //TODO: why go through all the elements? Most of them will not have Dirichlet bc.
            {
                //var elStart = DateTime.Now;
                IMatrix elementK = elementProvider.Matrix(element);

                //double[] localSolution = subdomain.CalculateElementDisplacements(element, solution);
                double[] localdSolution = subdomain.CalculateElementIncrementalConstraintDisplacements(element, constraintScalingFactor);

                var elementEquivalentForces = elementK.Multiply(localdSolution);

                subdomain.FreeDofOrdering.AddVectorElementToSubdomain(element, elementEquivalentForces, subdomainEquivalentForces);

                //times["addition"] += DateTime.Now - elStart;
            }

            //var totalTime = DateTime.Now - totalStart;

            return subdomainEquivalentForces;
        }        
    }
}
