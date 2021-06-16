using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;


namespace ISAAR.MSolve.Materials
{
    /// <summary>
    /// Deformation Gradient based implementation of a Mooney Rivlin hyperelastic material and a 
    /// NeoHookian hyperelastic material for 
    /// Authors Gerasimos Sotiropoulos
    /// </summary>
    public class HyperElasticMaterial3DDefGrad_new : IContinuumMaterial3DDefGrad
    {
        private readonly double[] strains = new double[6];
        private double[] stresses = new double[6];
        private Matrix constitutiveMatrix = null;
        public double C1 { get; set; }
        public double C2 { get; set; }

        public double k_cons { get; set; }
        public double[] Coordinates { get; set; }

        #region IFiniteElementMaterial Members

        public int ID
        {
            get { return 1; }
        }

        public bool Modified
        {
            get { return false; }
        }

        public void ResetModified()
        {
        }

        #endregion

        #region IFiniteElementMaterial3D Members

        public double[] Stresses { get { return (stresses); } }

        public IMatrixView ConstitutiveMatrix
        {
            get
            {
                if (constitutiveMatrix == null) UpdateMaterial(new double[9] { 1, 1, 1, 0, 0, 0, 0, 0, 0 });
                return constitutiveMatrix;
            }
        }

        public double YoungModulus => throw new NotImplementedException();

        public double PoissonRatio => throw new NotImplementedException();
        private double CalcDet2D(double[,] a)
        {
            return a[0, 0] * a[1, 1] - a[1, 0] * a[0, 1];
        }
        private double[,] ArrayTransposeMultiplyArray(double[,] a)
        {
            return new double[3, 3]
            {
                {a[0,0]*a[0,0] + a[1,0]*a[1,0] + a[2,0]*a[2,0], a[0,0]*a[0,1] + a[1,0]*a[1,1] + a[2,0]*a[2,1], a[0,0]*a[0,2] + a[1,0]*a[1,2] + a[2,0]*a[2,2]},
                {a[0,1]*a[0,0] + a[1,1]*a[1,0] + a[2,1]*a[2,0], a[0,1]*a[0,1] + a[1,1]*a[1,1] + a[2,1]*a[2,1], a[0,1]*a[0,2] + a[1,1]*a[1,2] + a[2,1]*a[2,2]},
                {a[0,2]*a[0,0] + a[1,2]*a[1,0] + a[2,2]*a[2,0], a[0,2]*a[0,1] + a[1,2]*a[1,1] + a[2,2]*a[2,1], a[0,2]*a[0,2] + a[1,2]*a[1,2] + a[2,2]*a[2,2]}
            };
        }
        private double[,] ArrayMultiplyArrayTranspose(double[] v)
        {
            var a00 = v[0];
            var a01 = v[3];
            var a02 = v[6];

            var a10 = v[7];
            var a11 = v[1];
            var a12 = v[4];

            var a20 = v[5];
            var a21 = v[8];
            var a22 = v[2];


            return new double[3, 3]
            {
                {a00*a00 + a01*a01 + a02*a02,
                 a00*a10+a01*a11+a02*a12,
                a00*a20+a01*a21+a02*a22},
                { a10*a00+a11*a01+a12*a02,
                  a10*a10+a11*a11+a12*a12,
                a10*a20+a11*a21+a12*a21
                },
                { a20*a00+a21*a01+a22*a02,
                a20*a10+a21*a11+a22*a12,
                a20*a20+a21*a21+a22*a22}
                };
        }
        //private double CalcDet3D(double[,] a)
        //{
        //    return a[0, 0] * (a[1, 1] * a[2, 2] - a[1, 2] * a[2, 1]) - a[0, 1] * (a[1, 0] * a[2, 2] - a[2, 0] * a[1, 2])
        //        + a[0, 2] * (a[0, 1] * a[2, 1] - a[2, 0] * a[1, 1]);
        //}
        private double CalcDet3D(double[,] a)
        {
            return a[0, 0] * a[1, 1] * a[2, 2] + a[0, 1] * a[1, 2] * a[2, 0] + a[0, 2] * a[1, 0] * a[2, 1]
                 - a[0, 0] * a[1, 2] * a[2, 1] - a[0, 1] * a[1, 0] * a[2, 2] - a[0, 2] * a[1, 1] * a[2, 0];
        }
        private double[,] Array6by6Inverse(double[,] a)
        {
            double[,] inv;
            double det = CalcDet3D(a);
            double[,] a00 = new double[,] { { a[1, 1], a[1, 2] }, { a[2, 1], a[2, 2] } };
            double[,] a01 = new double[,] { { a[1, 0], a[1, 2] }, { a[2, 0], a[2, 2] } };
            double[,] a02 = new double[,] { { a[1, 0], a[1, 1] }, { a[2, 0], a[2, 1] } };
            double[,] a10 = new double[,] { { a[0, 1], a[0, 2] }, { a[2, 1], a[2, 2] } };
            double[,] a11 = new double[,] { { a[0, 0], a[0, 2] }, { a[2, 0], a[2, 2] } };
            double[,] a12 = new double[,] { { a[0, 0], a[0, 1] }, { a[2, 0], a[2, 1] } };
            double[,] a20 = new double[,] { { a[0, 1], a[0, 2] }, { a[1, 1], a[1, 2] } };
            double[,] a21 = new double[,] { { a[0, 0], a[0, 2] }, { a[1, 0], a[1, 2] } };
            double[,] a22 = new double[,] { { a[0, 0], a[0, 1] }, { a[1, 0], a[1, 1] } };
            inv = new double[,]
            {
                {CalcDet2D(a00)/det, -CalcDet2D(a10)/det, CalcDet2D(a20)/det},
                {-CalcDet2D(a01)/det, CalcDet2D(a11)/det, -CalcDet2D(a21)/det},
                {CalcDet2D(a02)/det, -CalcDet2D(a12)/det, CalcDet2D(a22)/det},
            };
            return inv;
        }
        private double[,] Array3by3TimesScalar(double s, double[,] a)
        {
            return new double[,] { { s*a[0,0], s*a[0,1], s*a[0,2]},
            { s*a[1,0], s*a[1,1], s*a[1,2]},
            { s*a[2,0], s*a[2,1], s*a[2,2]}};
        }

        private double[,] Array6by6TimesScalar(double s, double[,] a)
        {
            return new double[,] { { s*a[0,0], s*a[0,1], s*a[0,2],
                    s*a[0,3], s*a[0,4],s*a[0,5]},
            { s*a[1,0], s*a[1,1], s*a[1,2],
                    s*a[1,3], s*a[1,4], s*a[1,5]},
            { s*a[2,0], s*a[2,1], s*a[2,2],
                    s*a[2,3], s*a[2,4], s*a[2,5]},
            { s*a[3,0], s*a[3,1], s*a[3,2],
                    s*a[3,3], s*a[3,4], s*a[3,5]},
            {s* a[4,0], s*a[4,1], s*a[4,2],
                    s*a[4,3], s*a[4,4], s*a[4,5]},
            { s*a[5,0], s*a[5,1], s*a[5,2],
                    s*a[5,3], s*a[5,4], s*a[5,5]}};
        }
        private double[,] Array3by3Sum(double[,] a, double[,] b)
        {
            return new double[,] { { a[0,0]+b[0,0], a[0,1]+b[0,1],a[0,2]+b[0,2]},
            { a[1,0]+b[1,0], a[1,1]+b[1,1],a[1,2]+b[1,2]},
            { a[2,0]+b[2,0], a[2,1]+b[2,1],a[2,2]+b[2,2]}};
        }
        private double[,] Array6by6Sum(double[,] a, double[,] b)
        {
            return new double[,] { { a[0,0]+b[0,0], a[0,1]+b[0,1],a[0,2]+b[0,2],a[0,3]+b[0,3],
                a[0,4]+b[0,4],a[0,5]+b[0,5]},
            { a[1,0]+b[1,0], a[1,1]+b[1,1],a[1,2]+b[1,2],a[1,3]+b[1,3],
                a[1,4]+b[1,4],a[1,5]+b[1,5]},
            { a[2,0]+b[2,0], a[2,1]+b[2,1],a[2,2]+b[2,2],a[2,3]+b[2,3],
                a[2,4]+b[2,4],a[2,5]+b[2,5]},
            { a[3,0]+b[3,0], a[3,1]+b[3,1],a[3,2]+b[3,2],a[3,3]+b[3,3],
                a[3,4]+b[3,4],a[3,5]+b[3,5]},
            { a[4,0]+b[4,0], a[4,1]+b[4,1],a[4,2]+b[4,2],a[4,3]+b[4,3],
                a[4,4]+b[4,4],a[4,5]+b[4,5]},
            { a[5,0]+b[5,0], a[5,1]+b[5,1],a[5,2]+b[5,2],a[5,3]+b[5,3],
                a[5,4]+b[5,4],a[5,5]+b[5,5]},};
        }
        private double[,] TensorProduct3by3(double[] a, double[] b)
        {
            return new double[,] { { a[0]*b[0], a[0] * b[1], a[0] * b[2]},
            { a[1]*b[0], a[1] * b[1], a[1] * b[2]},
            { a[2]*b[0], a[2] * b[1], a[2] * b[2]}};
        }
        private double[,] TensorProduct6by6(double[] a, double[] b)
        {
            return new double[,] { { a[0]*b[0], a[0] * b[1], a[0] * b[2],
                    a[0] * b[3], a[0] * b[4],a[0] * b[5]},
            { a[1]*b[0], a[1] * b[1], a[1] * b[2],
                    a[1] * b[3], a[1] * b[4], a[1] * b[5]},
            { a[2]*b[0], a[2] * b[1], a[2] * b[2],
                    a[2] * b[3], a[2] * b[4], a[2] * b[5]},
            { a[3]*b[0], a[3] * b[1], a[3] * b[2],
                    a[3] * b[3], a[3] * b[4], a[3] * b[5]},
            { a[4]*b[0], a[4] * b[1], a[4] * b[2],
                    a[4] * b[3], a[4] * b[4], a[4] * b[5]},
            { a[5]*b[0], a[5] * b[1], a[5] * b[2],
                    a[5] * b[3], a[5] * b[4], a[5] * b[5]}};
        }
        public void UpdateMaterial(double[] DefGradVec)
        {
            var deformationGradient = new double[3, 3] { { DefGradVec [0], DefGradVec[3], DefGradVec[6] },
                                                { DefGradVec [7], DefGradVec[1], DefGradVec[4] },
                                                { DefGradVec [5], DefGradVec[8], DefGradVec[2] }};

            //var rCG = deformationGradient.MultiplyRight(deformationGradient, true, false);
            var rCG = ArrayTransposeMultiplyArray(deformationGradient);
            var rbG = ArrayMultiplyArrayTranspose(DefGradVec);

            //var rCGvec = Vector.CreateFromArray(new double[] { rCG[0, 0], rCG[1, 1], rCG[2, 2], rCG[0, 1], rCG[1, 2], rCG[2, 0] });
            var rCG_inv = Array6by6Inverse(rCG);
            //var rCG_inv_Vec = Vector.CreateFromArray(new double[] { rCG_inv[0, 0], rCG_inv[1, 1], rCG_inv[2, 2], rCG_inv[0, 1], rCG_inv[1, 2], rCG_inv[2, 0] });

            var I_1 = rCG[0, 0] + rCG[1, 1] + rCG[2, 2];
            var vec_12 = new double[3]; var vec_22 = new double[3];
            var I_2 = 0.5 * (Math.Pow(I_1, 2) - (rCG[0, 0] * rCG[0, 0] + rCG[1, 1] * rCG[1, 1] + rCG[2, 2] * rCG[2, 2])
                - 2 * (rCG[0, 1] * rCG[0, 1] + rCG[1, 2] * rCG[1, 2] + rCG[2, 0] * rCG[2, 0]));
            var I_3 = CalcDet3D(rCG);
            var J_3 = CalcDet3D(deformationGradient);

            var I_1_st = new double[,] { { 2, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } };
            var I_2_st = new double[,] { {2 * (I_1-rCG[0,0]), -2 * rCG[0, 1], -2 * rCG[0, 2]},
                                              {-2 * rCG[1,0], 2*(I_1-rCG[1, 1]), -2 * rCG[1, 2]},
                                              {-2 * rCG[2,0], -2 * rCG[2,1], 2 * (I_1-rCG[2, 2])}};
            //var I_3_st = Array3by3TimesScalar(2 * I_3, rCG_inv);
            var I_3_st = new double[,] { {2 * I_3 * rCG_inv[0,0], 2 * I_3 *rCG_inv[0, 1], 2 * I_3 *rCG_inv[0, 2]},
                                              {2 * I_3 * rCG_inv[1,0], 2 * I_3 * rCG_inv[1,1], 2 * I_3 *rCG_inv[1, 2]},
                                              {2 * I_3 * rCG_inv[2,0], 2 * I_3 * rCG_inv[2,1], 2 * I_3 *rCG_inv[2, 2]}};
            var I_1_st_vec = new double[] { 2, 2, 2, 0, 0, 0 };
            var I_2_st_vec = new double[] { I_2_st[0, 0], I_2_st[1, 1], I_2_st[2, 2], I_2_st[0, 1], I_2_st[1, 2], I_2_st[2, 0] };
            var I_3_st_vec = new double[] { I_3_st[0, 0], I_3_st[1, 1], I_3_st[2, 2], I_3_st[0, 1], I_3_st[1, 2], I_3_st[2, 0] };

            double a = Math.Pow(I_3, (-1d / 3d));
            double b = (1d / 3d) * I_1 * Math.Pow(I_3, -4d / 3d);
            double c = Math.Pow(I_3, -2d / 3d);
            double d = (2d / 3d) * I_2 * Math.Pow(I_3, -5d / 3d);
            double e = 0.5 * Math.Pow(I_3, -1d / 2d);
            var aTimesI_1_st = Array3by3TimesScalar(a, I_1_st);
            var bTimesI_3_st = Array3by3TimesScalar(-b, I_3_st);
            var J_1_st = Array3by3Sum(aTimesI_1_st, bTimesI_3_st);
            //var J_1_st = new double[,] { { a * I_1_st[0, 0] - b * I_3_st[0, 0], a * I_1_st[0, 1] - b * I_3_st[0, 1], a * I_1_st[0, 2] - b * I_3_st[0, 2] },
            //                            { a * I_1_st[1, 0] - b * I_3_st[1, 0], a * I_1_st[1, 1] - b * I_3_st[1, 1], a * I_1_st[1, 2] - b * I_3_st[1, 2] },
            //                            { a * I_1_st[2, 0] - b * I_3_st[2, 0], a * I_1_st[2, 1] - b * I_3_st[2, 1], a * I_1_st[2, 2] - b * I_3_st[2, 2] }};

            var cTimesI_2_st = Array3by3TimesScalar(c, I_2_st);
            var dTimesI_3_st = Array3by3TimesScalar(-d, I_3_st);
            var J_2_st = Array3by3Sum(cTimesI_2_st, dTimesI_3_st);

            //var J_2_st = new double[,] { { c * I_2_st[0, 0] - d * I_3_st[0, 0], c * I_2_st[0, 1] - d * I_3_st[0, 1], c * I_2_st[0, 2] - d * I_3_st[0, 2] },
            //                            { c * I_2_st[1, 0] - d * I_3_st[1, 0], c * I_2_st[1, 1] - d * I_3_st[1, 1], c * I_2_st[1, 2] - d * I_3_st[1, 2] },
            //                            { c * I_2_st[2, 0] - d * I_3_st[2, 0], c * I_2_st[2, 1] - d * I_3_st[2, 1], c * I_2_st[2, 2] - d * I_3_st[2, 2] }};
            var J_3_st = Array3by3TimesScalar(e, I_3_st);

            //var J_3_st = new double[,] { { e * I_3_st[0, 0], e * I_3_st[0, 1], e * I_3_st[0, 2]},
            //                            { e * I_3_st[1, 0], e * I_3_st[1, 1], e * I_3_st[1, 2]},
            //                            { e * I_3_st[2, 0], e * I_3_st[2, 1], e * I_3_st[2, 2]}};

            var J_3_st_vec = new double[] { J_3_st[0, 0], J_3_st[1, 1], J_3_st[2, 2], J_3_st[0, 1], J_3_st[1, 2], J_3_st[2, 0] };

            var I_2_stst = new double[6, 6] { {0, 4, 4, 0, 0, 0},
                                            {  4, 0, 4, 0, 0, 0 },
                                            {  4, 4, 0, 0, 0, 0 },
                                            {  0, 0, 0,-2, 0, 0 },
                                            {  0, 0, 0, 0,-2, 0 },
                                            {  0, 0, 0, 0, 0,-2 } };

            var I_3_stst = new double[6, 6]{ {             0, 4 * rCG[2, 2], 4 * rCG[1, 1],             0, -4 * rCG[1, 2],            0 },
                                             { 4 * rCG[2, 2],             0, 4 * rCG[0, 0],             0,              0,-4 * rCG[2, 0] },
                                             { 4 * rCG[1, 1], 4 * rCG[0, 0],             0,-4 * rCG[0, 1],              0,            0 },
                                                         { 0,             0,-4 * rCG[0, 1],-2 * rCG[2, 2],  2 * rCG[2, 0], 2 * rCG[1, 2] },
                                            { -4 * rCG[1, 2],             0,             0, 2 * rCG[2, 0], -2 * rCG[0, 0], 2 * rCG[0, 1] },
                                                         { 0,-4 * rCG[2, 0],             0, 2 * rCG[1, 2],  2 * rCG[0, 1],-2 * rCG[1, 1] }};

            var f = (-1d / 3d) * Math.Pow(I_3, -4d / 3d);
            var g = (4d / 9d) * I_1 * Math.Pow(I_3, -7d / 3d);
            var I_13 = TensorProduct6by6(I_1_st_vec, I_3_st_vec);
            var I_31 = TensorProduct6by6(I_3_st_vec, I_1_st_vec);
            var I_33 = TensorProduct6by6(I_3_st_vec, I_3_st_vec);
            var I_23 = TensorProduct6by6(I_2_st_vec, I_3_st_vec);
            var gTimesI_33 = Array6by6TimesScalar(g, I_33);
            var I_1TimesI_3_stst = Array6by6TimesScalar(I_1, I_3_stst);
            var I_13PlusI_31 = Array6by6Sum(I_13, I_31);
            var I_13PLusI_31PlusI_1I_3_stst = Array6by6Sum(I_13PlusI_31, I_1TimesI_3_stst);
            var fI_13PLusI_31PlusI_1I_3_stst = Array6by6TimesScalar(f, I_13PLusI_31PlusI_1I_3_stst);
            var J_1_stst = Array6by6Sum(fI_13PLusI_31PlusI_1I_3_stst, gTimesI_33);
            var h = -(2d / 3d) * (Math.Pow(I_3, -5d / 3d));
            var i = (10d / 9d) * I_2 * Math.Pow(I_3, -8d / 3d);
            var cTimesI_2_stst = Array6by6TimesScalar(c, I_2_stst);
            var I_32 = TensorProduct6by6(I_3_st_vec, I_2_st_vec);
            var I_23PlusI_32 = Array6by6Sum(I_23, I_32);
            var I_2TimesI_3_stst = Array6by6TimesScalar(I_2, I_3_stst);
            var I_23PlusI_32PlusI_2TimesI_3_stst = Array6by6Sum(I_23PlusI_32, I_2TimesI_3_stst);
            var hTimesI_23PlusI_32PlusI_2I_3_stst = Array6by6TimesScalar(h, I_23PlusI_32PlusI_2TimesI_3_stst);
            var iTimesI_33 = Array6by6TimesScalar(i, I_33);
            var cTimesI_2_ststPlusiI_33 = Array6by6Sum(cTimesI_2_stst, iTimesI_33);
            var J_2_stst = Array6by6Sum(hTimesI_23PlusI_32PlusI_2I_3_stst, cTimesI_2_ststPlusiI_33);
            var j = (-1d / 4d) * Math.Pow(I_3, -3d / 2d);
            var k = (1d / 2d) * Math.Pow(I_3, -1d / 2d);
            var jTimesI_33 = Array6by6TimesScalar(j, I_33);
            var kTimesI_3_stst = Array6by6TimesScalar(k, I_3_stst);
            var J_3_stst = Array6by6Sum(jTimesI_33, kTimesI_3_stst);

            var s1 = Array3by3TimesScalar(C1, J_1_st);
            var s2 = Array3by3TimesScalar(C2, J_2_st);
            var s3 = Array3by3TimesScalar(k_cons * (J_3 - 1d), J_3_st);
            var Spk = Array3by3Sum(s1, s2);
            Spk = Array3by3Sum(Spk, s3);

            var Spk_vec = new double[] { Spk[0, 0], Spk[1, 1], Spk[2, 2], Spk[0, 1], Spk[1, 2], Spk[2, 0] };

            //einai to Cklrs_pavla_vec
            var ss1 = Array6by6TimesScalar(C1, J_1_stst);
            var ss2 = Array6by6TimesScalar(C2, J_2_stst);
            var J_33 = TensorProduct6by6(J_3_st_vec, J_3_st_vec);
            var ss3 = Array6by6TimesScalar(k_cons, J_33);
            var ss4 = Array6by6TimesScalar(J_3 - 1, J_3_stst);
            var Cons = Array6by6Sum(ss1, ss2);
            Cons = Array6by6Sum(Cons, ss3);
            Cons = Array6by6Sum(Cons, ss4);

            stresses = Spk_vec.Copy();
            constitutiveMatrix = Matrix.CreateFromArray(Cons);
        }

        public void ClearState()
        {
            //throw new NotImplementedException();
        }

        public void SaveState()
        {
            //throw new NotImplementedException();
        }

        public void ClearStresses()
        {
            //throw new NotImplementedException();
        }

        #endregion

        #region ICloneable Members

        public object Clone()
        {
            return new HyperElasticMaterial3DDefGrad() { C1 = this.C1, C2 = this.C2, k_cons = this.k_cons };
        }

        object ICloneable.Clone()
        {
            return Clone();
        }

        #endregion

    }
}