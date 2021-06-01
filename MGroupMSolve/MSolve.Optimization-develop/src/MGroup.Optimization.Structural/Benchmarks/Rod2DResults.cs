using MGroup.FEM.Elements;
using MGroup.FEM.Entities;
using MGroup.Solvers.LinearSystems;

namespace MGroup.Optimization.Structural.Benchmarks
{
    public class Rod2DResults
    {
        private readonly Subdomain subdomain;
        private readonly ILinearSystem linearSystem;

        public Rod2DResults(Subdomain subdomain, ILinearSystem linearSystem)
        {
            this.subdomain = subdomain;
            this.linearSystem = linearSystem;
        }

        public double AxialRod2DStress(Element element)
        {
            double[] localDisplacements = 
                subdomain.FreeDofOrdering.ExtractVectorElementFromSubdomain(element, linearSystem.Solution);
            Rod2D rod = (Rod2D)element.ElementType;
            return  rod.CalculateAxialStress(element, localDisplacements, null);
        }
    }
}
