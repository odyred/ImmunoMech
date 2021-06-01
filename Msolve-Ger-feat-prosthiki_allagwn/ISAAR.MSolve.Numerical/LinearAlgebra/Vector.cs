﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System.Threading;
using System.IO;
using System.Globalization;
using System.Collections.Specialized;
using System.Diagnostics;
using ISAAR.MSolve.Numerical;
using ISAAR.MSolve.Numerical.Interfaces;

namespace ISAAR.MSolve.Numerical.LinearAlgebra
{
    public class Vector : IVector, IFileWriteable
    {
        private readonly double[] data;

        public Vector(int length)
        {
            data = new double[length];
        }

        public Vector(double[] data)
        {
            this.data = data;
        }

        public double[] Data
        {
            get { return data; }
        }

        public double Norm
        {
            get
            {
                double norm = 0;
                double[] v1Data = data as double[];
                for (int i = 0; i < this.Length; i++)
                    norm += v1Data[i] * v1Data[i];
                return Math.Sqrt(norm);
            }
        }

        public static Vector operator *(double s1, Vector v2)
        {
            double[] v2Data = v2.data as double[];
            double[] result = new double[v2.Length];
            for (int i = 0; i < v2.Length; i++) result[i] = s1 * v2Data[i];
            return (new Vector(result)) as Vector;
        }

        private static double DoDot(int segment, double[] d1, double[] d2)
        {
            int l = (int)(d1.Length / VectorExtensions.AffinityCount);
            int start = segment * l;
            int finish = (segment + 1) * l;
            if (segment == VectorExtensions.AffinityCount - 1) finish = d1.Length;
            double result = 0;
            for (int i = start; i < finish; i++) result += d1[i] * d2[i];
            return result;
        }

        public static double operator *(Vector v1, Vector v2)
        {
            double result = 0;
            double[] v1Data = v1.data;
            double[] v2Data = v2.data;
            for (int i = 0; i < v1.Length; i++) result += v1Data[i] * v2Data[i];

            //int iProcs = VectorExtensions.AffinityCount;
            //double[] results = new double[iProcs];
            //Parallel.For(0, iProcs, i => { results[i] = DoDot(i, v1Data, v2Data); });
            //for (int i = 0; i < iProcs; i++) result += results[i];
            return result;
        }

        public static double[] operator +(Vector v1, Vector v2)
        {
            double[] result = new double[v1.Length];
            double[] v1Data = v1.data;
            double[] v2Data = v2.data;
            for (int i = 0; i < v1.Length; i++) result[i] = v1Data[i] + v2Data[i];
            return result;
        }

        public static double[] operator -(Vector v1, Vector v2)
        {
            double[] result = new double[v1.Length];
            double[] v1Data = v1.data;
            double[] v2Data = v2.data;
            for (int i = 0; i < v1.Length; i++) result[i] = v1Data[i] - v2Data[i];
            return result;
        }

        public void Scale(double scale)
        {
            double[] vData = data;
            for (int i = 0; i < data.Length; i++) vData[i] *= scale;
        }

        public void Add(Vector v)
        {
            double[] v1Data = data;
            double[] v2Data = v.Data;
            for (int i = 0; i < data.Length; i++) v1Data[i] += v2Data[i];
        }

        public void Subtract(Vector v)
        {
            double[] v1Data = data;
            double[] v2Data = v.Data;
            for (int i = 0; i < data.Length; i++) v1Data[i] -= v2Data[i];
        }

        #region IVector<T> Members

        public int Length
        {
            get { return data.Length; }
        }

        public double this[int x]
        {
            get { return data[x]; }
            set { data[x] = value; }
        }

        public double DotProduct(IVector y)
        {
            double result = 0;
            for (int i = 0; i < data.Length; i++)
                result += data[i] * y[i];
            return result;
        }

        public void Multiply(double coefficient)
        {
            double[] vData = data as double[];
            for (int i = 0; i < this.Length; i++) vData[i] *= coefficient;
        }

        public void CopyTo(Array array, int index)
        {
            data.CopyTo(array, index);
        }

        public void CopyFrom(int startIndex, int length, IVector fromVector, int fromStartIndex)
        {
            for (int i = 0; i < length; i++)
                data[i + startIndex] = fromVector[i + fromStartIndex];
        }

        public void Clear()
        {
            Array.Clear(data, 0, data.Length);
        }

        public void WriteToFile(string name)
        {
            double[] mData = data as double[];

            StreamWriter sw = new StreamWriter(name);
            foreach (double d in mData)
                sw.WriteLine(d.ToString("g17", new CultureInfo("en-US", false).NumberFormat));
            sw.Close();
        }

        #endregion
    }
}
