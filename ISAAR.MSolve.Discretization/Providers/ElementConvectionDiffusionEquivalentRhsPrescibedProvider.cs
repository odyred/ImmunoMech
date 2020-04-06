using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

//TODO: Should this be in the Problems project?
namespace ISAAR.MSolve.Discretization.Providers
{
    public class ElementConvectionDiffusionEquivalentRhsPrescibedProvider : IBoundaryElementMatrixProvider
    {
        public IMatrix Matrix(IConvectionDiffusionBoundaryElement element) => element.RHSPrescribedMatrix(element);

        //IMatrix IElementMatrixProvider.Matrix(IElement element)
        //{
        //    throw new System.NotImplementedException();
        //}
    }
}
