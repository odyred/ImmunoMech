﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation.Inverse;
using ISAAR.MSolve.Geometry.Coordinates;

//TODO: cache the natural shape functions at gauss point sets
namespace ISAAR.MSolve.FEM.Interpolation
{
    /// <summary>
    /// Isoparametric interpolation of a triangular finite element with 3 nodes. Linear shape functions. 
    /// Implements Singleton pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class InterpolationTri3 : IsoparametricInterpolation2DBase
    {
        private static readonly InterpolationTri3 uniqueInstance = new InterpolationTri3();

        private InterpolationTri3() : base(3)
        {
            NodalNaturalCoordinates = new NaturalPoint2D[]
            {
                new NaturalPoint2D(0.0, 0.0),
                new NaturalPoint2D(1.0, 0.0),
                new NaturalPoint2D(0.0, 1.0)
            };
        }

        /// <summary>
        /// The coordinates of the finite element's nodes in the natural (element local) coordinate system. The order of these
        /// nodes matches the order of the shape functions and is always the same for each element.
        /// </summary>
        public override IReadOnlyList<NaturalPoint2D> NodalNaturalCoordinates { get; }

        /// <summary>
        /// Get the unique <see cref="InterpolationQuad4"/> object for the whole program. Thread safe.
        /// </summary>
        public static InterpolationTri3 UniqueInstance { get { return uniqueInstance; } }

        /// <summary>
        /// The inverse mapping of this interpolation, namely from global cartesian to natural (element local) coordinate system.
        /// </summary>
        /// <param name="nodes">The nodes of the finite element in the global cartesian coordinate system.</param>
        /// <returns></returns>
        public override IInverseInterpolation2D CreateInverseMappingFor(IReadOnlyList<Node2D> nodes)
        {
            return new InverseInterpolationTri3(nodes);
        }

        protected override sealed double[] EvaluateAt(double xi, double eta)
        {
            var values = new double[3];
            values[0] = 1 - xi - eta;
            values[1] = xi;
            values[2] = eta;
            return values;
        }

        protected override sealed double[,] EvaluateGradientsAt(double xi, double eta)
        {
            var derivatives = new double[3, 2];
            derivatives[0, 0] = -1;
            derivatives[0, 1] = -1;
            derivatives[1, 0] = 1;
            derivatives[1, 1] = 0;
            derivatives[2, 0] = 0;
            derivatives[2, 1] = 1;
            return derivatives;
        }
    }
}
