using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using System.Collections.Generic;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IConvectionDiffusionElement: IFiniteElement
    {
        IMatrix MassTransportConductivityMatrix(IElement element);
        IMatrix DiffusionConductivityMatrix(IElement element);
        IMatrix SecondSpaceDerivativeXMatrix(IElement element);
        IMatrix SecondSpaceDerivativeYMatrix(IElement element);
        IMatrix SecondSpaceDerivativeZMatrix(IElement element);
    }
}
