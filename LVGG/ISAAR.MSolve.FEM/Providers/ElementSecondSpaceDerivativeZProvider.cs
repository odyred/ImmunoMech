﻿using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.FEM.Interfaces;
using System.Collections.Generic;

namespace ISAAR.MSolve.FEM.Providers
{
    public class ElementSecondSpaceDerivativeZProvider : IElementMatrixProvider
    {
        public IMatrix Matrix(IElement element)
        {
            IConvectionDiffusionElement elementType = (IConvectionDiffusionElement)element.ElementType;
            return elementType.SecondSpaceDerivativeZMatrix(element);
        }
    }
}
