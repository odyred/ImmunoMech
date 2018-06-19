﻿using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;
using Xunit;

namespace ISAAR.MSolve.Tests.FEM
{
    /// <summary>
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ContinuumElementTests
    {
        private static readonly ContinuumElement2DFactory factory = new ContinuumElement2DFactory();

        private static readonly IReadOnlyList<Node2D> nodeSet1 = new Node2D[]
        {
            new Node2D(0, -1.0, -1.0),
            new Node2D(1, +1.0, -1.0),
            new Node2D(2, +1.0, +1.0),
            new Node2D(3, -1.0, +1.0)
        };

        private static readonly IReadOnlyList<Node2D> nodeSet2 = new Node2D[]
        {
            new Node2D(0,  0.0,  0.0),
            new Node2D(1, 20.0,  0.0),
            new Node2D(2, 20.0, 10.0),
            new Node2D(3,  0.0, 10.0)
        };

        private static readonly IReadOnlyList<Node2D> nodeSet3 = new Node2D[]
        {
            new Node2D(0, 0.2, 0.3),
            new Node2D(1, 2.2, 1.5),
            new Node2D(2, 3.0, 2.7),
            new Node2D(3, 0.7, 2.0)
        };

        private static readonly ElasticMaterial2D material1 = new ElasticMaterial2D
        {
            YoungModulus = 2e6,
            PoissonRatio = 0.3,
            Thickness = 1.0,
            StressState = "plstress"
        };

        private static readonly ElasticMaterial2D material2 = new ElasticMaterial2D
        {
            YoungModulus = 1.0,
            PoissonRatio = 0.25,
            Thickness = 1.0,
            StressState = "plstress"
        };

        [Fact]
        private static void TestQuad4Stiffness1()
        {
            ContinuumElement2D quad4 = factory.CreateQuad4(nodeSet1, material1);
            IMatrix2D K = quad4.StiffnessMatrix(null);
            double[,] expectedK = new double[,]
            {
                {  989010.98901099,  357142.85714286, -604395.60439561, -27472.527472528, -494505.49450549, -357142.85714286,  109890.10989011,  27472.527472528 },
                {  357142.85714286,  989010.98901099,  27472.527472528,  109890.10989011, -357142.85714286, -494505.49450549, -27472.527472528, -604395.60439561 },
                { -604395.60439561,  27472.527472528,  989010.98901099, -357142.85714286,  109890.10989011, -27472.527472528, -494505.49450549,  357142.85714286 },
                { -27472.527472528,  109890.10989011, -357142.85714286,  989010.98901099,  27472.527472528, -604395.60439561,  357142.85714286, -494505.49450549 },
                { -494505.49450549, -357142.85714286,  109890.10989011,  27472.527472528,  989010.98901099,  357142.85714286, -604395.60439561, -27472.527472528 },
                { -357142.85714286, -494505.49450549, -27472.527472528, -604395.60439561,  357142.85714286,  989010.98901099,  27472.527472528,  109890.10989011 },
                {  109890.10989011, -27472.527472528, -494505.49450549,  357142.85714286, -604395.60439561,  27472.527472528,  989010.98901099, -357142.85714286 },
                {  27472.527472528, -604395.60439561,  357142.85714286, -494505.49450549, -27472.527472528,  109890.10989011, -357142.85714286,  989010.98901099 }
            };
            Assert.True(AreMatricesEqual(K, new Matrix2D(expectedK), 1e-10));
        }

        [Fact]
        private static void TestQuad4Stiffness2()
        {
            ContinuumElement2D quad4 = factory.CreateQuad4(nodeSet2, material1);
            IMatrix2D K = quad4.StiffnessMatrix(null);
            double[,] expectedK = new double[,]
            {
                {  879120.87912088,  357142.85714286, -109890.10989011, -27472.527472527, -439560.43956044, -357142.85714286, -329670.32967033,  27472.527472527 },
                {  357142.85714286,  1593406.5934066,  27472.527472528,  604395.60439560, -357142.85714286, -796703.29670329, -27472.527472528, -1401098.9010989 },
                { -109890.10989011,  27472.527472528,  879120.87912088, -357142.85714286, -329670.32967033, -27472.527472527, -439560.43956044,  357142.85714286 },
                { -27472.527472527,  604395.60439560, -357142.85714286,  1593406.5934066,  27472.527472528, -1401098.9010989,  357142.85714286, -796703.29670329 },
                { -439560.43956044, -357142.85714286, -329670.32967033,  27472.527472528,  879120.87912088,  357142.85714286, -109890.10989011, -27472.527472527 },
                { -357142.85714286, -796703.29670329, -27472.527472527, -1401098.9010989,  357142.85714286,  1593406.5934066,  27472.527472528,  604395.60439560 },
                { -329670.32967033, -27472.527472528, -439560.43956044,  357142.85714286, -109890.10989011,  27472.527472528,  879120.87912088, -357142.85714286 },
                {  27472.527472527, -1401098.9010989,  357142.85714286, -796703.29670329, -27472.527472527,  604395.60439560, -357142.85714286,  1593406.5934066 }
            };
            Assert.True(AreMatricesEqual(K, new Matrix2D(expectedK), 1e-10));
        }

        [Fact]
        private static void TestQuad4Stiffness3()
        {
            ContinuumElement2D quad4 = factory.CreateQuad4(nodeSet3, material1);
            IMatrix2D K = quad4.StiffnessMatrix(null);
            double[,] expectedK = new double[,]
            {
                {  514198.06499808, -10764.170693892, -403744.28140248,  6179.6240659003,  136202.13267487, -257206.34711690, -246655.91627047,  261790.89374489 },
                { -10764.170693892,  877762.34956189,  61124.679010955,  241708.83734231, -257206.34711690, -50430.336321830,  206845.83879984, -1069040.8505824 },
                { -403744.28140248,  61124.679010955,  2618366.6203953, -1268741.0140479, -648721.49301375,  372500.34071924, -1565900.8459791,  835115.99431770 },
                {  6179.6240659003,  241708.83734231, -1268741.0140479,  3119360.6646083,  427445.39566429, -1580482.4587671,  835115.99431770, -1780587.0431835 },
                {  136202.13267487, -257206.34711690, -648721.49301375,  427445.39566429,  691579.93709070, -83847.039187745, -179060.57675182, -86392.009359646 },
                { -257206.34711690, -50430.336321830,  372500.34071924, -1580482.4587671, -83847.039187745,  1103398.3531728, -31446.954414592,  527514.44191615 },
                { -246655.91627047,  206845.83879984, -1565900.8459791,  835115.99431770, -179060.57675182, -31446.954414592,  1991617.3390014, -1010514.8787030 },
                {  261790.89374489, -1069040.8505824,  835115.99431770, -1780587.0431835, -86392.009359646,  527514.44191615, -1010514.8787030,  2322113.4518497 }
            };
            Assert.True(AreMatricesEqual(K, new Matrix2D(expectedK), 1e-10));
        }

        [Fact]
        private static void TestQuad4Stiffness4()
        {
            ContinuumElement2D quad4 = factory.CreateQuad4(nodeSet3, material2);
            IMatrix2D K = quad4.StiffnessMatrix(null);
            double[,] expectedK = new double[,]
            {
                { 0.2592059570522744,       -0.005023279657148771,  -0.1906544880785286,    -0.017628995948735127,  0.06474697564228854,    -0.12002962865455298,   -0.13329844461603435,   0.1426819042604369      },
                { -0.005023279657148771,    0.42886928984872097,    0.049037670717931546,   0.11055696733570626,    -0.12002962865455298,   -0.022348176556173396,  0.0760152375937702,     -0.5170780806282538     },
                { -0.1906544880785286,      0.049037670717931546,   1.3012408988907096,     -0.5920791398890163,    -0.33356025755449087,   0.15332067182282197,    -0.77702615325769,      0.3897207973482628      },
                { -0.017628995948735127,    0.11055696733570626,    -0.5920791398890163,    1.5350381195234322,     0.21998733848948865,    -0.7683820415727372,    0.3897207973482628,     -0.8772130452864011     },
                { 0.06474697564228854,      -0.12002962865455298,   -0.33356025755449087,   0.21998733848948865,    0.3475567568780642,     -0.03912861828761283,   -0.07874347496586184,   -0.060829091547322856   },
                { -0.12002962865455298,     -0.022348176556173396,  0.15332067182282197,    -0.7683820415727372,    -0.03912861828761283,   0.5397386843830521,     0.005837575119343828,   0.25099153374585864     },
                { -0.13329844461603435,     0.0760152375937702,     -0.77702615325769,      0.3897207973482628,     -0.07874347496586184,   0.005837575119343828,   0.9890680728395862,     -0.4715736100613768     },
                { 0.1426819042604369,       -0.5170780806282538,    0.3897207973482628,     -0.8772130452864011,    -0.060829091547322856,  0.25099153374585864,    -0.4715736100613768,    1.1432995921687963      }
            };
            Assert.True(AreMatricesEqual(K, new Matrix2D(expectedK), 1e-10));
        }

        private static bool AreMatricesEqual(IMatrix2D matrix1, IMatrix2D matrix2, double tolerance)
        {
            if ((matrix1.Rows != matrix2.Rows) || (matrix1.Columns != matrix2.Columns)) return false;
            for (int i = 0; i < matrix1.Rows; ++i)
            {
                for (int j = 0; j < matrix1.Columns; ++j)
                {
                    if (!AreValuesEqual(matrix1[i, j], matrix2[i, j], tolerance)) return false;
                }
            }
            return true;
        }

        private static bool AreValuesEqual(double value1, double value2, double tolerance)
        {
            if (Math.Abs(value2) <= tolerance) // Can't divide with expected ~= 0. 
            {
                if (Math.Abs(value1) <= tolerance) return true;
                else return false;
            }
            else return (Math.Abs(1.0 - value1 / value2) < tolerance) ? true : false;
        }
    }
}
