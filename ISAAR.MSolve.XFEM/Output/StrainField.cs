﻿using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;

namespace ISAAR.MSolve.XFEM.Output
{
    //TODO: Would be better to do it with deformation matrices, but requires unification of vectors and matrices
    class StrainField : IOutputField
    {
        public Tensor2D EvaluateAt(XContinuumElement2D element, NaturalPoint point,
            Vector standardDisplacements, Vector enrichedDisplacements)
        {
            EvalInterpolation2D evaluatedInterpolation = element.Interpolation.EvaluateAllAt(element.Nodes, point);
            Matrix2by2 displacementGradient = element.CalculateDisplacementFieldGradient(
                point, evaluatedInterpolation, standardDisplacements, enrichedDisplacements);

            double strainXX = displacementGradient[0, 0];
            double strainYY = displacementGradient[1, 1];
            double strainXY = 0.5 * (displacementGradient[0, 1] + displacementGradient[1, 0]);

            return new Tensor2D(strainXX, strainYY, strainXY);
        }
    }
}
