using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IConvectionDiffusionBoundaryElement : IElement
    {
        IMatrix RHSPrescribedMatrix(IConvectionDiffusionBoundaryElement element);
        IMatrix RHSFLuxMatrix(IElement element);
    }
}
