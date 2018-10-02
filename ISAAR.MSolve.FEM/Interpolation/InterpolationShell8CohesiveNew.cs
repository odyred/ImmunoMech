﻿using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interpolation
{
    public class InterpolationShell8CohesiveNew: IsoparametricInterpolation2DBase
    {
        private static readonly InterpolationShell8CohesiveNew uniqueInstance = new InterpolationShell8CohesiveNew();

        private readonly Dictionary<IQuadrature2D, IReadOnlyList<double[,]>> cachedN3AtGPs;

        private InterpolationShell8CohesiveNew() : base(CellType2D.Quad8,8)
        {
            cachedN3AtGPs = new Dictionary<IQuadrature2D, IReadOnlyList<double[,]>>();
            NodalNaturalCoordinates = new NaturalPoint2D[]
            {
                //TODO: validate this
                new NaturalPoint2D(1, 1),
                new NaturalPoint2D(-1, 1),
                new NaturalPoint2D(-1, -1),
                new NaturalPoint2D(1, -1),
                new NaturalPoint2D(0, 1),
                new NaturalPoint2D(-1, 0),
                new NaturalPoint2D(0, -1),
                new NaturalPoint2D(1, 0)
            };
        }
        public override IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        public static InterpolationShell8CohesiveNew UniqueInstance => uniqueInstance;

        public override IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes) => throw new NotImplementedException();

        protected override double[] EvaluateAt(double ksi, double heta)
        {
            double[] N1gp = new double[8]; //8=nShapeFunctions;

            N1gp[4] = 0.5 * (1 - Math.Pow(ksi, 2)) * (1 + heta);
            N1gp[5] = 0.5 * (1 - Math.Pow(heta, 2)) * (1 - ksi);
            N1gp[6] = 0.5 * (1 - Math.Pow(ksi, 2)) * (1 - heta);
            N1gp[7] = 0.5 * (1 - Math.Pow(heta, 2)) * (1 + ksi);

            N1gp[0] = 0.25 * (1 + ksi) * (1 + heta) - 0.5 * N1gp[4] - 0.5 * N1gp[7];
            N1gp[1] = 0.25 * (1 - ksi) * (1 + heta) - 0.5 * N1gp[4] - 0.5 * N1gp[5];
            N1gp[2] = 0.25 * (1 - ksi) * (1 - heta) - 0.5 * N1gp[5] - 0.5 * N1gp[6];
            N1gp[3] = 0.25 * (1 + ksi) * (1 - heta) - 0.5 * N1gp[6] - 0.5 * N1gp[7];

            return N1gp;
        }

        protected override double[,] EvaluateGradientsAt(double ksi, double heta)
        {
            var shapeFunctionDerivativesGp = new double[2, 8]; // 2(0:ws pros ksi 1:ws pros heta) 8(arithmos shape functions ara kai komvwn)

            shapeFunctionDerivativesGp[0, 4] = (-ksi) * (1 + heta);
            shapeFunctionDerivativesGp[0, 5] = -0.5 * (1 - Math.Pow(heta, 2));
            shapeFunctionDerivativesGp[0, 6] = 0.5 * (-2 * ksi) * (1 - heta);
            shapeFunctionDerivativesGp[0, 7] = 0.5 * (1 - Math.Pow(heta, 2));
            shapeFunctionDerivativesGp[0, 0] = +0.25 * (1 + heta) - 0.5 * shapeFunctionDerivativesGp[0, 4] - 0.5 * shapeFunctionDerivativesGp[0, 7];
            shapeFunctionDerivativesGp[0, 1] = -0.25 * (1 + heta) - 0.5 * shapeFunctionDerivativesGp[0, 4] - 0.5 * shapeFunctionDerivativesGp[0, 5];
            shapeFunctionDerivativesGp[0, 2] = -0.25 * (1 - heta) - 0.5 * shapeFunctionDerivativesGp[0, 5] - 0.5 * shapeFunctionDerivativesGp[0, 6];
            shapeFunctionDerivativesGp[0, 3] = +0.25 * (1 - heta) - 0.5 * shapeFunctionDerivativesGp[0, 6] - 0.5 * shapeFunctionDerivativesGp[0, 7];

            shapeFunctionDerivativesGp[1, 4] = 0.5 * (1 - Math.Pow(ksi, 2));
            shapeFunctionDerivativesGp[1, 5] = 0.5 * (-2 * heta) * (1 - ksi);
            shapeFunctionDerivativesGp[1, 6] = 0.5 * (1 - Math.Pow(ksi, 2)) * (-1);
            shapeFunctionDerivativesGp[1, 7] = 0.5 * (-2 * heta) * (1 + ksi);
            shapeFunctionDerivativesGp[1, 0] = +0.25 * (1 + ksi) - 0.5 * shapeFunctionDerivativesGp[1, 4] - 0.5 * shapeFunctionDerivativesGp[1, 7];
            shapeFunctionDerivativesGp[1, 1] = +0.25 * (1 - ksi) - 0.5 * shapeFunctionDerivativesGp[1, 4] - 0.5 * shapeFunctionDerivativesGp[1, 5];
            shapeFunctionDerivativesGp[1, 2] = -0.25 * (1 - ksi) - 0.5 * shapeFunctionDerivativesGp[1, 5] - 0.5 * shapeFunctionDerivativesGp[1, 6];
            shapeFunctionDerivativesGp[1, 3] = -0.25 * (1 + ksi) - 0.5 * shapeFunctionDerivativesGp[1, 6] - 0.5 * shapeFunctionDerivativesGp[1, 7];

           return shapeFunctionDerivativesGp;
        }

        public IReadOnlyList<double[,]> EvaluateN3ShapeFunctionsReorganized(IQuadrature2D quadrature)
        {
            bool isCached = cachedN3AtGPs.TryGetValue(quadrature,
                out IReadOnlyList<double[,]> N3AtGPs);
            if (isCached) return N3AtGPs;
            else
            {
                IReadOnlyList<Vector> N1 = EvaluateFunctionsAtGaussPoints(quadrature);
                N3AtGPs = GetN3ShapeFunctionsReorganized(quadrature, N1);
                cachedN3AtGPs.Add(quadrature, N3AtGPs);
                return N3AtGPs;
            }
        }

        private IReadOnlyList<double[,]> GetN3ShapeFunctionsReorganized(IQuadrature2D quadrature, IReadOnlyList<Vector> N1)
        {
            //TODO reorganize cohesive shell  to use only N1 (not reorganised)

            int nGaussPoints = quadrature.IntegrationPoints.Count;
            var N3 = new double[nGaussPoints][,]; // shapeFunctionsgpData
            for (int npoint = 0; npoint < nGaussPoints; npoint++)
            {
                double ksi = quadrature.IntegrationPoints[npoint].Xi;
                double heta = quadrature.IntegrationPoints[npoint].Eta;
                double[,] N3gp = new double[3, 24]; ; //8=nShapeFunctions;

                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 8; m++)
                    { N3gp[l, l + 3 * m] = N1[npoint][m]; }
                }
                N3[npoint] = N3gp;
            }
            return N3;

        }

    }
}
