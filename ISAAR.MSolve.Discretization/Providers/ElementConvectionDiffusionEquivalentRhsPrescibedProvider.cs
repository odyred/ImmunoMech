using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.FEM.Interfaces;

//TODO: Should this be in the Problems project?
namespace ISAAR.MSolve.Discretization.Providers
{
    public class ElementConvectionDiffusionEquivalentRhsPrescibedProvider : IBoundaryElementMatrixProvider
    {
        public IMatrix Matrix(IElement element) => element.ElementType.StiffnessMatrix(element);

        //IMatrix IElementMatrixProvider.Matrix(IElement element)
        //{
        //    throw new System.NotImplementedException();
        //}
    }
}
