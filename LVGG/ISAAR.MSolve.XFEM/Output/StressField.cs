﻿using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Output
{
    class StressField: IOutputField
    {
        public Tensor2D EvaluateAt(XContinuumElement2D element, NaturalPoint point,
            Vector standardDisplacements, Vector enrichedDisplacements)
        {
            EvalInterpolation2D evaluatedInterpolation = element.Interpolation.EvaluateAllAt(element.Nodes, point);
            Matrix2by2 displacementGradient = element.CalculateDisplacementFieldGradient(
                point, evaluatedInterpolation, standardDisplacements, enrichedDisplacements);
            Matrix constitutive =
                element.Material.CalculateConstitutiveMatrixAt(point, evaluatedInterpolation);

            return element.CalculateStressTensor(displacementGradient, constitutive);
        }
    }
}
