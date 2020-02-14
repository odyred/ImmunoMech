using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.FEM.Interpolation.Jacobians;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Elements.BoundaryConditionElements
{
    public class SurfaceBoundary3D : IFiniteElement,  ICell<Node>
    {
        private readonly IDofType[][] dofTypes; //TODO: this should not be stored for each element. Instead store it once for each Quad4, Tri3, etc. Otherwise create it on the fly.
        private readonly ThermalMaterial material;
        //private readonly Dictionary<GaussPoint, ThermalMaterial> materialsAtGaussPoints;


        public SurfaceBoundary3D(double thickness, IReadOnlyList<Node> nodes, IIsoparametricInterpolation2D interpolation,
            IQuadrature2D quadratureForStiffness, IQuadrature2D quadratureForConsistentMass,
            IGaussPointExtrapolation2D gaussPointExtrapolation,
            ThermalMaterial material)
        {
            this.material = material;
            this.GaussPointExtrapolation = gaussPointExtrapolation;
            this.Nodes = nodes;
            this.Interpolation = interpolation;
            this.QuadratureForConsistentMass = quadratureForConsistentMass;
            this.QuadratureForStiffness = quadratureForStiffness;
            this.Thickness = thickness;

            dofTypes = new IDofType[nodes.Count][];
            for (int i = 0; i < interpolation.NumFunctions; ++i) dofTypes[i] = new IDofType[] { ThermalDof.Temperature };
        }

        public CellType CellType => Interpolation.CellType;
        public ElementDimensions ElementDimensions => ElementDimensions.TwoD;

        public int ID => throw new NotImplementedException(
            "Element type codes should be in a settings class. Even then it's a bad design choice");

        public IGaussPointExtrapolation2D GaussPointExtrapolation { get; }
        public IIsoparametricInterpolation2D Interpolation { get; }
        public IReadOnlyList<Node> Nodes { get; }
        public IQuadrature2D QuadratureForConsistentMass { get; }
        public IQuadrature2D QuadratureForStiffness { get; }
        public double Thickness { get; }

        public bool MaterialModified => throw new NotImplementedException();

        public IElementDofEnumerator DofEnumerator { get; set; } = new GenericDofEnumerator();

        public IMatrix MassMatrix(IElement element)
        {
            return BuildCapacityMatrix();
        }

        public Matrix BuildCapacityMatrix()
        {
            int numDofs = Nodes.Count;
            var capacity = Matrix.CreateZero(numDofs, numDofs);
            //IReadOnlyList<double[]> shapeFunctions =
            //    Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            //IReadOnlyList<Matrix> shapeGradientsNatural =
            //    Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            //for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            //{
            //    Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
            //    Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
            //    var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
            //    double dA = jacobian.DirectDeterminant * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
            //    capacity.AxpyIntoThis(partial, dA);
            //}

            ////WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            //capacity.ScaleIntoThis(Thickness * material.Density * material.SpecialHeatCoeff);
            return capacity;
        }
        public Matrix BuildConvectionMatrix()
        {
            int numDofs = Nodes.Count;
            var convection = Matrix.CreateZero(numDofs, numDofs);
            IReadOnlyList<double[]> shapeFunctions =
                Interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForConsistentMass);
            IReadOnlyList<Matrix> shapeGradientsNatural =
                Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);

            for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
            {
                Matrix shapeFunctionMatrix = BuildShapeFunctionMatrix(shapeFunctions[gp]);
                Matrix partial = shapeFunctionMatrix.Transpose() * shapeFunctionMatrix;
                Matrix jacobianMatrix = Matrix.CreateZero(2, 3);
                double xGaussPoint = 0;
                double yGaussPoint = 0;
                double zGaussPoint = 0;
                for (int k = 0; k < this.Nodes.Count; k++)
                {
                    //xGaussPoint += shapeFunctionMatrix[gp, k] * this.Nodes[k].X;
                    //yGaussPoint += shapeFunctionMatrix[gp, k] * this.Nodes[k].Y;
                    //zGaussPoint += shapeFunctionMatrix[gp, k] * this.Nodes[k].Z;
                    jacobianMatrix[0, 0] += shapeGradientsNatural[gp][k, 0] * this.Nodes[k].X;
                    jacobianMatrix[0, 1] += shapeGradientsNatural[gp][k, 0] * this.Nodes[k].Y;
                    jacobianMatrix[0, 2] += shapeGradientsNatural[gp][k, 0] * this.Nodes[k].Z;
                    jacobianMatrix[1, 0] += shapeGradientsNatural[gp][k, 1] * this.Nodes[k].X;
                    jacobianMatrix[1, 1] += shapeGradientsNatural[gp][k, 1] * this.Nodes[k].Y;
                    jacobianMatrix[1, 2] += shapeGradientsNatural[gp][k, 1] * this.Nodes[k].Z;
                }

                Vector surfaceBasisVector1 = Vector.CreateZero(3);
                surfaceBasisVector1[0] = jacobianMatrix[0, 0];
                surfaceBasisVector1[1] = jacobianMatrix[0, 1];
                surfaceBasisVector1[2] = jacobianMatrix[0, 2];

                Vector surfaceBasisVector2 = Vector.CreateZero(3);
                surfaceBasisVector2[0] = jacobianMatrix[1, 0];
                surfaceBasisVector2[1] = jacobianMatrix[1, 1];
                surfaceBasisVector2[2] = jacobianMatrix[1, 2];

                Vector surfaceBasisVector3 = surfaceBasisVector1.CrossProduct(surfaceBasisVector2);
                var jacdet = surfaceBasisVector3.Norm2();
                
//                var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
                double dA = jacdet * QuadratureForConsistentMass.IntegrationPoints[gp].Weight;
                convection.AxpyIntoThis(partial, dA);
            }

            //WARNING: the following needs to change for non uniform density. Perhaps the integration order too.
            convection.ScaleIntoThis(material.ThermalConvection);
            return convection;
        }

        //public Matrix BuildConductivityMatrix()
        //{
        //    int numDofs = Nodes.Count;
        //    var conductivity = Matrix.CreateZero(numDofs, numDofs);
        //    IReadOnlyList<Matrix> shapeGradientsNatural =
        //        Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

        //    for (int gp = 0; gp < QuadratureForStiffness.IntegrationPoints.Count; ++gp)
        //    {
        //        // Calculate the necessary quantities for the integration
        //        //Matrix constitutive = (Matrix)(materialsAtGaussPoints[gp].ConstitutiveMatrix); // ugly cast will be removed along with the legacy Matrix classes
        //        var jacobian = new IsoparametricJacobian3D(Nodes, shapeGradientsNatural[gp]);
        //        Matrix shapeGradientsCartesian =
        //            jacobian.TransformNaturalDerivativesToCartesian(shapeGradientsNatural[gp]);
        //        Matrix deformation = BuildDeformationMatrix(shapeGradientsCartesian);

        //        // Contribution of this gauss point to the element stiffness matrix
        //        Matrix partialK = deformation.Transpose() * deformation;
        //        //Matrix partialΚ = deformation.Transpose() * (constitutive * deformation);
        //        //partialK.Scale(materialsAtGaussPoints[gaussPoint].ThermalConductivity);

        //        double dA = jacobian.DirectDeterminant * QuadratureForStiffness.IntegrationPoints[gp].Weight; //TODO: this is used by all methods that integrate. I should cache it.
        //        conductivity.AxpyIntoThis(partialK, dA * material.ThermalConductivity);
        //    }

        //    conductivity.ScaleIntoThis(Thickness);
        //    return conductivity;
        //}

        // Provatidis uses two distinct vectors K = N,x^T * k * N,x + N,y^T * k * N,y
        //private (Matrix dNdX, Matrix dNdY) CalcdNdx(EvalShapeGradients2D shapeGrad)
        //{
        //    int n = Nodes.Count;
        //    var dNdX = new double[n, 1];
        //    var dNdY = new double[n, 1];
        //    for (int i = 0; i < n; ++i)
        //    {
        //        dNdX[i, 0] = shapeGrad[i][0];
        //        dNdY[i, 0] = shapeGrad[i][1];
        //    }
        //    return (new Matrix(dNdX), new Matrix(dNdY));
        //}

        private Matrix BuildDeformationMatrix(Matrix shapeGradientsCartesian)
        {
            //TODO: isn't this just the transpose of [dNi/dxj]?
            var deformation = Matrix.CreateZero(2, Nodes.Count);
            for (int nodeIdx = 0; nodeIdx < Nodes.Count; ++nodeIdx)
            {
                deformation[0, nodeIdx] = shapeGradientsCartesian[nodeIdx, 0];
                deformation[1, nodeIdx] = shapeGradientsCartesian[nodeIdx, 1];
            }
            return deformation;
        }

        /// <summary>
        /// The shape function matrix is 1-by-n, where n = is the number of shape functions.
        /// </summary>
        public Matrix BuildShapeFunctionMatrix(double[] shapeFunctions) //TODO: reconsider this. As it is, it just returns the shape functions in a Matrix
        {
            //var array2D = new double[1, shapeFunctions.Length];
            //for (int i = 0; i < shapeFunctions.Length; ++i)
            //{
            //    array2D[0, i] = shapeFunctions[i];
            //}
            //return Matrix.CreateFromArray(array2D);
            return Matrix.CreateFromArray(shapeFunctions, 1, shapeFunctions.Length);
        }

        public IReadOnlyList<IReadOnlyList<IDofType>> GetElementDofTypes(IElement element) => dofTypes;

        public void ResetMaterialModified()
        {
            throw new NotImplementedException();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            throw new NotImplementedException();
        }

        public void SaveMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialState()
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialStresses()
        {
            throw new NotImplementedException();
        }

        public IMatrix StiffnessMatrix(IElement element)
        {
            return BuildConvectionMatrix();
        }

        public IMatrix DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
        }


    }
}
