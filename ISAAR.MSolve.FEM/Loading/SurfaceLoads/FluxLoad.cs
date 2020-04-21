using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Loading.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Loading.SurfaceLoads
{
    public class FluxLoad : ISurfaceLoad
    {
        private readonly double _flux;

        public FluxLoad(double flux)
        {
            _flux = flux;
        }

        public Table<INode, IDofType, double> CalculateSurfaceLoad(IIsoparametricInterpolation2D interpolation, IQuadrature2D integration, IReadOnlyList<Node> nodes)
        {
            int numDofs = nodes.Count;
            var stiffness = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                interpolation.EvaluateFunctionsAtGaussPoints(integration);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                interpolation.EvaluateNaturalGradientsAtGaussPoints(integration);

            for (int gp = 0; gp < integration.IntegrationPoints.Count; ++gp)
            {
                //var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                Vector shapeFunctionMatrix =Vector.CreateFromArray(shapeFunctions[gp]);
                Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
                Vector IV = Vector.CreateWithValue(nodes.Count, 1);
                Matrix IM = Matrix.CreateIdentity(nodes.Count);
                for (int k = 0; k < nodes.Count; k++)
                {
                    //xGaussPoint += shapeFunctionMatrix[gp, k] * this.Nodes[k].X;
                    //yGaussPoint += shapeFunctionMatrix[gp, k] * this.Nodes[k].Y;
                    //zGaussPoint += shapeFunctionMatrix[gp, k] * this.Nodes[k].Z;
                    jacobianMatrix[0, 0] += shapeGradientsNatural[gp][k, 0] * nodes[k].X;
                    jacobianMatrix[0, 1] += shapeGradientsNatural[gp][k, 0] * nodes[k].Y;
                    jacobianMatrix[0, 2] += shapeGradientsNatural[gp][k, 0] * nodes[k].Z;
                    jacobianMatrix[1, 0] += shapeGradientsNatural[gp][k, 1] * nodes[k].X;
                    jacobianMatrix[1, 1] += shapeGradientsNatural[gp][k, 1] * nodes[k].Y;
                    jacobianMatrix[1, 2] += shapeGradientsNatural[gp][k, 1] * nodes[k].Z;
                }
                Vector tangentVector1 = jacobianMatrix.GetRow(0);
                Vector tangentVector2 = jacobianMatrix.GetRow(1);
                Vector normalVector = tangentVector1.CrossProduct(tangentVector2);

                Matrix partial = IV.TensorProduct(shapeFunctionMatrix) * IM;

                //Vector surfaceBasisVector1 = Vector.CreateZero(3);
                //surfaceBasisVector1[0] = jacobianMatrix[0, 0];
                //surfaceBasisVector1[1] = jacobianMatrix[0, 1];
                //surfaceBasisVector1[2] = jacobianMatrix[0, 2];

                //Vector surfaceBasisVector2 = Vector.CreateZero(3);
                //surfaceBasisVector2[0] = jacobianMatrix[1, 0];
                //surfaceBasisVector2[1] = jacobianMatrix[1, 1];
                //surfaceBasisVector2[2] = jacobianMatrix[1, 2];

                //Vector surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
                //var jacdet = surfaceBasisVector3.Norm2();
                var jacdet = normalVector.Norm2();

                double dA = jacdet * integration.IntegrationPoints[gp].Weight;
                stiffness.AxpyIntoThis(partial, dA);
            }
            var appliedDisplacements = Vector.CreateWithValue(nodes.Count, _flux);
            var fluxLoad = stiffness * appliedDisplacements;

            var table = new Table<INode, IDofType, double>();
            for (int i = 0; i < nodes.Count; i++)
            {
                table.TryAdd(nodes[i], ThermalDof.Temperature, fluxLoad[i]);
            }
            return table;
        }
    }
}
