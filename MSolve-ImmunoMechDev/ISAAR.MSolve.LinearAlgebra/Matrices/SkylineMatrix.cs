﻿using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using ISAAR.MSolve.LinearAlgebra.Commons;
using ISAAR.MSolve.LinearAlgebra.Exceptions;
using ISAAR.MSolve.LinearAlgebra.Triangulation;
using ISAAR.MSolve.LinearAlgebra.Output.Formatting;
using ISAAR.MSolve.LinearAlgebra.Providers.Managed;
using ISAAR.MSolve.LinearAlgebra.Reduction;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using static ISAAR.MSolve.LinearAlgebra.LibrarySettings;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using System.Linq;
using ISAAR.MSolve.LinearAlgebra.Reordering;

//TODO: Also linear combinations with other matrix types may be useful, e.g. Skyline (K) with diagonal (M), but I think 
//      that for global matrices, this should be done through concrete class to use DoEntrywiseIntoThis methods. 
//TODO: Checks like: col - row <= colHeight can be written more efficiently without calculating the height:
//      entryOffset = diagOffsets[col] + col - row, entryOffset <= diagOffsets[col+1]
//TODO: Throw MatrixDataOverwrittenException whenever necessary. Possibly make the values and diagOffsets readonly again.
//TODO: In most algorithms, I can cache in local variables for each column the height, diagOffset and diagOffset + colIdx
//TODO: Most algorithms implemented here should be moved to a class the holds the implementations and called from there.
namespace ISAAR.MSolve.LinearAlgebra.Matrices
{
    /// <summary>
    /// Symmetric sparse matrix stored in Skyline format (3-array version). Only the non-zero entries of the upper triangle are 
    /// stored. The Skyline format is optimized for Cholesky factorizations. To build a <see cref="SkylineMatrix"/> conveniently, 
    /// use <see cref="Builders.SkylineBuilder"/>.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class SkylineMatrix: IMatrix, ISparseMatrix, ISymmetricMatrix
    {
        /// <summary>
        /// Contains the non-zero entries of the matrix's upper triangle in column major order, starting from the diagonal and 
        /// going upwards. Its length is nnz.
        /// </summary>
        private double[] values;

        /// <summary>
        /// Contains the indices into values of the diagonal entries of the matrix. Its length = order + 1, with the last entry
        /// being equal to nnz.
        /// </summary>
        private int[] diagOffsets;

        private bool isOverwritten = false;

        private SkylineMatrix(int order, double[] values, int[] diagOffsets)
        {
            this.values = values;
            this.diagOffsets = diagOffsets;
            this.NumColumns = order;
        }

        /// <summary>
        /// The number of columns of the matrix. 
        /// </summary>
        public int NumColumns { get; }

        /// <summary>
        /// The number of rows of the matrix.
        /// </summary>
        public int NumRows { get { return NumColumns; } }

        /// <summary>
        /// The internal array that stores the non-zero entries of the matrix's upper triangle in column major order, 
        /// starting from the diagonal and going upwards. Its length is equal to the number of non-zero entries. 
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal double[] RawValues => values;

        /// <summary>
        /// The internal array that stores the indices into <see cref="RawValues"/> of the diagonal entries of the matrix. 
        /// Its length = order + 1, with the last entry being equal to nnz.
        /// It should only be used for passing the raw array to linear algebra libraries.
        /// </summary>
        internal int[] RawDiagOffsets => diagOffsets;

        /// <summary>
        /// See <see cref="IIndexable2D.this[int, int]"/>.
        /// </summary>
        public double this[int rowIdx, int colIdx] // TODO: should I add index bound checking?
        {
            get
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                int diagOffset = diagOffsets[colIdx];
                int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffset - 1; // excluding diagonal
                int entryHeight = colIdx - rowIdx; // excluding diagonal
                if (entryHeight > maxColumnHeight) return 0.0; // outside stored non zero pattern
                else return values[diagOffset + entryHeight];
            }
            set
            {
                ProcessIndices(ref rowIdx, ref colIdx);
                int diagOffset = diagOffsets[colIdx];
                int maxColumnHeight = diagOffsets[colIdx + 1] - diagOffset - 1; // excluding diagonal
                int entryHeight = colIdx - rowIdx; // excluding diagonal
                if (entryHeight > maxColumnHeight)
                {
                    throw new SparsityPatternModifiedException($"In column {colIdx} only rows [{maxColumnHeight}, {colIdx}]"
                        + $" can be changed, but you are trying to set entry ({rowIdx}, {colIdx})");
                }
                else values[diagOffset + entryHeight] = value;
            }
        }

        /// <summary>
        /// Initializes a new <see cref="SkylineMatrix"/> that contains the non zero entries of the upper triangle of 
        /// <paramref name="original"/>. In skyline format, some zero entries will be explicitly stored.
        /// </summary>
        /// <param name="original">The matrix that will be copied. It must be symmetric.</param>
        /// <param name="tolerance">
        /// The tolerance used to determine if an entry is zero. It will also be used to check, if <paramref name="original"/>
        /// is symmetric (<paramref name="original"/>[i, j] == <paramref name="original"/>[j, i]).
        /// </param>
        public static SkylineMatrix CreateFromArray(double[,] original, double tolerance = 1E-10)
        {
            if (!original.IsSymmetric(tolerance)) throw new ArgumentException("The original matrix must be symmetric.");

            // Indexing array
            var comparer = new ValueComparer(tolerance);
            int order = original.GetLength(0);
            var diagOffsets = new int[order + 1];
            //diagOffsets[0] = 0; // by default;
            for (int j = 0; j < order; ++j)
            {
                int colHeight = 0;
                for (int i = j - 1; i >= 0; --i)
                {
                    if (!comparer.AreEqual(0.0, original[i, j])) colHeight = j - i;
                }
                diagOffsets[j + 1] = diagOffsets[j] + colHeight + 1;
            }

            // Values array
            int nnz = diagOffsets[order];
            var values = new double[nnz];
            for (int j = 0; j < order; ++j)
            {
                int colHeight = diagOffsets[j + 1] - diagOffsets[j] - 1;
                for (int t = 0; t <= colHeight; ++t)
                {
                    int i = j - t; // row index of Aij
                    int offsetAij = diagOffsets[j] + t;
                    values[offsetAij] = original[i, j];
                }
            }

            return new SkylineMatrix(order, values, diagOffsets);
        }

        /// <summary>
        /// Initializes a new <see cref="SkylineMatrix"/> with the specified dimensions and the provided arrays 
        /// (<paramref name="values"/> and <paramref name="diagOffsets"/>) as its internal data.
        /// </summary>
        /// <param name="order">The number of rows/columns of the new matrix.</param>
        /// <param name="values">Contains the non zero superdiagonal entries of the matrix in column major order, starting from 
        ///     the diagonal and going upwards.</param>
        /// <param name="diagOffsets">Contains the indices into <paramref name="values"/> of the diagonal entries of the matrix. 
        ///     Its length is <paramref name="order"/> + 1, with the last entry being equal to nnz.</param>
        /// <param name="checkInput">If true, the provided arrays will be checked to make sure they are valid Skyline arrays, 
        ///     which is safer. If false, no such check will take place, which is faster.</param>
        /// <param name="copyArrays">If true, the provided arrays will be copied and the new <see cref="SkylineMatrix"/> instance 
        ///     will have references to the copies, which is safer. If false, the new matrix will have references to the 
        ///     provided arrays themselves, which is faster.</param>
        public static SkylineMatrix CreateFromArrays(int order, double[] values, int[] diagOffsets, 
            bool checkInput, bool copyArrays = false)
        {
            if (checkInput)
            {
                if (diagOffsets.Length != order + 1)
                {
                    throw new ArgumentException("The length of the Skyline diagonal offsets array must be equal to the number of"
                        + " rows/columns + 1, but was " + diagOffsets.Length);
                }
                if (diagOffsets[diagOffsets.Length - 1] != values.Length)
                {
                    throw new ArgumentException("The last entry of the Skyline diagonal offsets array must be equal to the number"
                        + " of non zero entries, but was " + diagOffsets[diagOffsets.Length - 1]);
                }
            }
            if (copyArrays)
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                int[] diagOffsetsCopy = new int[diagOffsets.Length];
                Array.Copy(diagOffsets, diagOffsetsCopy, diagOffsets.Length);
                return new SkylineMatrix(order, valuesCopy, diagOffsetsCopy);
            }
            else return new SkylineMatrix(order, values, diagOffsets);
        }

        /// <summary>
        /// Initializes a new <see cref="SkylineMatrix"/> that contains the non zero entries of the upper triangle of 
        /// <paramref name="original"/>. In skyline format, some zero entries will be explicitly stored.
        /// </summary>
        /// <param name="original">The matrix that will be copied. It must be symmetric.</param>
        /// <param name="tolerance">
        /// The tolerance used to determine if an entry is zero. It will also be used to check, if <paramref name="original"/>
        /// is symmetric (<paramref name="original"/>[i, j] == <paramref name="original"/>[j, i]).
        /// </param>
        public static SkylineMatrix CreateFromMatrix(IIndexable2D original, double tolerance = 1E-10)
        {
            if (!original.IsSymmetric(tolerance)) throw new ArgumentException("The original matrix must be symmetric.");
            
            // Indexing array
            var comparer = new ValueComparer(tolerance);
            int order = original.NumColumns;
            var diagOffsets = new int[order + 1];
            //diagOffsets[0] = 0; // by default;
            for (int j = 0; j < order; ++j)
            {
                int colHeight = 0;
                for (int i = j - 1; i >= 0; --i)
                {
                    if (!comparer.AreEqual(0.0, original[i, j])) colHeight = j - i;
                }
                diagOffsets[j+1] = diagOffsets[j] + colHeight + 1;
            }

            // Values array
            int nnz = diagOffsets[order];
            var values = new double[nnz];
            for (int j = 0; j < order; ++j)
            {
                int colHeight = diagOffsets[j + 1] - diagOffsets[j] - 1;
                for (int t = 0; t <= colHeight; ++t)
                {
                    int i = j - t; // row index of Aij
                    int offsetAij = diagOffsets[j] + t;
                    values[offsetAij] = original[i, j];
                }
            }

            return new SkylineMatrix(order, values, diagOffsets);
        }

        /// <summary>
        /// Initializes a new <see cref="SkylineMatrix"/> with the specified dimensions and the sparsity pattern defined by 
        /// <paramref name="diagOffsets"/>. The stored entries will initially be 0.
        /// </summary>
        /// <param name="order">The number of rows/columns of the new matrix.</param>
        /// <param name="diagOffsets">Contains the indices into <paramref name="values"/> of the diagonal entries of the matrix. 
        ///     Its length is <paramref name="order"/> + 1, with the last entry being equal to nnz.</param>
        /// <param name="checkInput">If true, <paramref name="diagOffsets"/> will be checked to make sure it is a valid Skyline  
        ///     array, which is safer. If false, no such check will take place, which is faster.</param>
        public static SkylineMatrix CreateZeroWithPattern(int order, int[] diagOffsets, bool checkInput)
        {
            if (checkInput)
            {
                if (diagOffsets.Length != order + 1)
                {
                    throw new ArgumentException("The length of the Skyline diagonal offsets array must be equal to the number of"
                        + " rows/columns + 1, but was " + diagOffsets.Length);
                }
            }
            int nnz = diagOffsets[diagOffsets.Length] - 1;
            return new SkylineMatrix(order, new double[nnz], diagOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Axpy(IMatrixView, double)"/>.
        /// </summary>
        public IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is SkylineMatrix otherSKY) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherSKY))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    Array.Copy(this.values, resultValues, values.Length);
                    Blas.Daxpy(values.Length, otherCoefficient, otherSKY.values, 0, 1, resultValues, 0, 1);
                    return new SkylineMatrix(NumColumns, resultValues, this.diagOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, 1.0, otherMatrix, otherCoefficient);
        }

        /// <summary>
        /// Performs the following operation for 0 &lt;= j &lt; <see cref="NumColumns"/>, 0 &lt;= i &lt;= j:
        /// result[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix is written to a new <see cref="SkylineMatrix"/> and then returned.
        /// </summary>
        /// <param name="otherMatrix">A matrix with the same indexing array as this <see cref="SkylineMatrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="SparsityPatternModifiedException">Thrown if <paramref name="otherMatrix"/> has different 
        ///     indexing array than this instance.</exception>
        public SkylineMatrix Axpy(SkylineMatrix otherMatrix, double otherCoefficient)
        {
            // Conceptually it is not wrong to so this, even if the indexers are different, but how would I implement it.
            if (!HasSameIndexer(otherMatrix))
            {
                throw new SparsityPatternModifiedException("Only allowed if the indexing arrays are the same");
            }
            //TODO: Perhaps this should be done using mkl_malloc and BLAS copy. 
            double[] resultValues = new double[values.Length];
            Array.Copy(this.values, resultValues, values.Length);
            Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, resultValues, 0, 1);
            // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
            return new SkylineMatrix(NumColumns, resultValues, this.diagOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrix.AxpyIntoThis(IMatrixView, double)"/>.
        /// </summary>
        public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is SkylineMatrix casted) AxpyIntoThis(casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                 "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// Performs the following operation for th non-zero entries (i, j), such that 0 &lt;= j &lt; <see cref="NumColumns"/> 
        /// 0 &lt;= i &lt;= j:
        /// this[i, j] = <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j] + this[i, j]. 
        /// The resulting matrix overwrites the entries of this <see cref="CscMatrix"/> instance.
        /// </summary>
        /// <param name="otherMatrix">A matrix for whose active columns (non-zero part of the column) are shorter or equal to the
        ///     active columns of this <see cref="SkylineMatrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="SparsityPatternModifiedException">Thrown if <paramref name="otherMatrix"/> has taller active columns 
        ///     than this <see cref="SkylineMatrix"/> instance.</exception>
        public void AxpyIntoThis(SkylineMatrix otherMatrix, double otherCoefficient)
        {
            if (HasSameIndexer(otherMatrix)) // no need to check dimensions if the indexing arrays are the same
            {
                Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, this.values, 0, 1);
            }
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                for (int j = 0; j < NumColumns; ++j)
                {
                    // Check if the current column of this matrix is tall enough
                    int thisDiagOffset = this.diagOffsets[j];
                    int otherDiagOffset = otherMatrix.diagOffsets[j];
                    int thisColTop = j - this.diagOffsets[j + 1] + thisDiagOffset + 1;
                    int otherColTop = j - this.diagOffsets[j + 1] + otherDiagOffset + 1;
                    if (thisColTop < otherColTop) throw new SparsityPatternModifiedException(
                        $"Column {j} of this matrix is shorter, which would result in overflow");

                    // Do the operation between the two columns. The column of this matrix is taller or equal to the other one.
                    for (int i = j; i >= otherColTop; --i) // non zero entries of shortest=other column, including diagonal
                    {
                        this.values[thisDiagOffset + j - i] += otherCoefficient * otherMatrix.values[otherDiagOffset + j - i];
                    }
                    // Don't do anything to the non zero entries of the above the shortest=other column (this[i,j] += a*0)
                }
            }
        }

        /// <summary>
        /// See <see cref="IMatrix.Clear"/>.
        /// </summary>
        public void Clear() => Array.Clear(values, 0, values.Length);

        /// <summary>
        /// See <see cref="IMatrixView.Copy(bool)"/>.
        /// </summary>
        IMatrix IMatrixView.Copy(bool copyIndexingData) => Copy(copyIndexingData);

        /// <summary>
        /// Copies the entries of this matrix.
        /// </summary>
        /// <param name="copyIndexingData">
        /// If true, all data of this object will be copied. If false, only the array containing the values of the stored 
        /// matrix entries will be copied. The new matrix will reference the same indexing arrays as this one.
        /// </param>
        public SkylineMatrix Copy(bool copyIndexingData)
        {
            var valuesCopy = new double[this.values.Length];
            Array.Copy(this.values, valuesCopy, this.values.Length);

            if (!copyIndexingData) return new SkylineMatrix(NumColumns, valuesCopy, this.diagOffsets);
            else
            {
                var diagOffsetsCopy = new int[this.diagOffsets.Length];
                Array.Copy(this.diagOffsets, diagOffsetsCopy, this.diagOffsets.Length);
                return new SkylineMatrix(NumColumns, valuesCopy, diagOffsetsCopy);
            }
        }

        /// <summary>
        /// Copies the entries of the matrix into a 2-dimensional array. The returned array has length(0) = <see cref="NumRows"/> 
        /// and length(1) = <see cref="NumColumns"/>. 
        /// </summary>
        public double[,] CopyToArray2D()
        {
            double[,] array2D = new double[NumColumns, NumColumns];
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
                array2D[j, j] = values[colOffset]; // diagonal entry
                for (int i = columnTop; i < j; ++i) // non zero entries stored above diagonal
                {
                    double value = values[colOffset + j - i];
                    array2D[j, i] = value;
                    array2D[i, j] = value;
                }
            }
            return array2D;
        }

        /// <summary>
        /// See <see cref="IMatrixView.CopyToFullMatrix()"/>
        /// </summary>
        public Matrix CopyToFullMatrix()
        {
            Matrix fullMatrix = Matrix.CreateZero(this.NumColumns, this.NumColumns);
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
                fullMatrix[j, j] = values[colOffset]; // diagonal entry
                for (int i = columnTop; i < j; ++i) // non zero entries stored above diagonal
                {
                    double value = values[colOffset + j - i];
                    fullMatrix[j, i] = value;
                    fullMatrix[i, j] = value; 
                }
            }
            return fullMatrix;
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.CountNonZeros"/>
        /// </summary>
        public int CountNonZeros() => values.Length;

        /// <summary>
        /// See <see cref="IEntrywiseOperableView2D{TMatrixIn, TMatrixOut}.DoEntrywise(TMatrixIn, Func{double, double, double})"/>.
        /// </summary>
        public IMatrix DoEntrywise(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is SkylineMatrix otherSKY) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherSKY))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    for (int i = 0; i < values.Length; ++i)
                    {
                        resultValues[i] = binaryOperation(this.values[i], otherSKY.values[i]);
                    }
                    return new SkylineMatrix(NumColumns, resultValues, diagOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.DoEntrywise(this, other, binaryOperation);
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperable2D{TMatrixIn}.DoEntrywiseIntoThis(TMatrixIn, Func{double, double, double})"/>.
        /// </summary>
        public void DoEntrywiseIntoThis(IMatrixView other, Func<double, double, double> binaryOperation)
        {
            if (other is SkylineMatrix sky)
            {
                if (HasSameIndexer(sky)) // no need to check dimensions if the indexing arrays are the same
                {
                    for (int i = 0; i < values.Length; ++i) this.values[i] = binaryOperation(this.values[i], sky.values[i]);
                }
                else
                {
                    Preconditions.CheckSameMatrixDimensions(this, other);
                    for (int j = 0; j < NumColumns; ++j)
                    {
                        // Check if the current column of this matrix is tall enough
                        int thisDiagOffset = this.diagOffsets[j];
                        int otherDiagOffset = sky.diagOffsets[j];
                        int thisColTop = j - this.diagOffsets[j + 1] + thisDiagOffset + 1;
                        int otherColTop = j - this.diagOffsets[j + 1] + otherDiagOffset + 1;
                        if (thisColTop < otherColTop) throw new SparsityPatternModifiedException(
                            $"Column {j} of this matrix is shorter, which would result in overflow");

                        // Do the operation between the two columns. The column of this matrix is taller or equal to the other one.
                        for (int i = j; i >= otherColTop; --i) // non zero entries of shortest column, including diagonal
                        {
                            int thisIndex = thisDiagOffset + j - i;
                            this.values[thisIndex] = binaryOperation(this.values[thisIndex], sky.values[otherDiagOffset + j - i]);
                        }
                        for (int i = otherColTop - 1; i >= thisColTop; --i) // non zero entries of the above the shortest column
                        {
                            int thisIndex = thisDiagOffset + j - i;
                            this.values[thisIndex] = binaryOperation(this.values[thisIndex], sky.values[otherDiagOffset + j - i]);
                        }
                    }
                }
            }
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperableView2D{TMatrixIn, TMatrixOut}.DoToAllEntries(Func{double, double})"/>.
        /// </summary>
        public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
        {
            // Only apply the operation on non zero entries
            double[] newValues = new double[values.Length];
            for (int i = 0; i < values.Length; ++i) newValues[i] = unaryOperation(values[i]);

            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0)) // The same sparsity pattern can be used.
            {
                // Copy the index arrays. TODO: See if we can use the same index arrays (e.g. if this class does not change them (it shouldn't))
                int[] diagOffsetsCopy = new int[diagOffsets.Length];
                Array.Copy(diagOffsets, diagOffsetsCopy, diagOffsets.Length);
                return new SkylineMatrix(NumColumns, newValues, diagOffsetsCopy);
            }
            else // The sparsity is destroyed. Revert to a full matrix.
            {
                return new SkylineMatrix(NumColumns, newValues, diagOffsets).CopyToFullMatrix();
            }
        }

        /// <summary>
        /// See <see cref="IEntrywiseOperable2D{TMatrixIn}.DoToAllEntriesIntoThis(Func{double, double})"/>.
        /// </summary>
        public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
        {
            if (new ValueComparer(1e-10).AreEqual(unaryOperation(0.0), 0.0))
            {
                for (int i = 0; i < values.Length; ++i) values[i] = unaryOperation(values[i]);
            }
            else
            {
                throw new SparsityPatternModifiedException("This operation will change the sparsity pattern");
            }
        }

        /// <summary>
        /// See <see cref="ISparseMatrix.EnumerateNonZeros"/>.
        /// </summary>
        public IEnumerable<(int row, int col, double value)> EnumerateNonZeros()
        {
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j + 1] + colOffset + 1;
                yield return (j, j, values[colOffset]); // diagonal entry
                for (int i = columnTop; i < j; ++i)
                {
                    double value = values[colOffset + j - i];
                    yield return (i, j, value);
                    yield return (j, i, value);
                }
            }
        }

        /// <summary>
        /// See <see cref="IIndexable2D.Equals(IIndexable2D, double)"/>.
        /// </summary>
        public bool Equals(IIndexable2D other, double tolerance = 1E-13)
        {
            if ((this.NumRows != other.NumRows) || (this.NumColumns != other.NumColumns)) return false;
            var comparer = new ValueComparer(1e-13);
            for (int j = 0; j < NumColumns; ++j)
            {
                int colOffset = diagOffsets[j];
                int columnTop = j - diagOffsets[j+1] + colOffset + 1;
                for (int i = 0; i < columnTop; ++i) // zero entries above stored column
                {
                    if (!( comparer.AreEqual(0.0, other[i, j]) && comparer.AreEqual(0.0, other[j, i]) )) return false;
                }
                for (int i = columnTop; i < j; ++i) // non zero entries of column, excluding diafonal
                {
                    double value = values[colOffset + j - i];
                    if (!(comparer.AreEqual(value, other[i, j]) && comparer.AreEqual(value, other[j, i]))) return false;
                }
                if (!comparer.AreEqual(values[colOffset], other[j, j])) return false; // non zero diagonal entry
            }
            return true; // At this point all entries have been checked and are equal
        }

        /// <summary>
        /// Calculate the Cholesky factorization. The matrix must be positive definite, otherwise an
        /// <see cref="IndefiniteMatrixException"/> will be thrown. If <paramref name="inPlace"/> is set to true, this object 
        /// must not be used again, otherwise a <see cref="NullReferenceException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">
        /// False, to copy the internal non zero entries before factorization. True, to overwrite them with the factorized data, 
        /// thus saving memory and time. However, that will make this object unusable, so you MUST NOT call any other members 
        /// afterwards.
        /// </param>
        /// <param name="pivotTolerance">
        /// If a diagonal entry is closer to zero than this tolerance, an <see cref="IndefiniteMatrixException"/> exception will
        /// be thrown.
        /// </param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not positive definite.</exception>
        /// <exception cref="NullReferenceException">
        /// Thrown if a member of his instance is accessed after this method is called.
        /// </exception>"
        public CholeskySkyline FactorCholesky(bool inPlace, double tolerance = CholeskySkyline.PivotTolerance)
        {
            if (inPlace)
            {
                var factor = CholeskySkyline.Factorize(NumColumns, values, diagOffsets, tolerance);
                // Set the skyline arrays to null to force NullReferenceException if they are accessed again.
                // TODO: perhaps there is a better way to handle this.
                values = null;
                diagOffsets = null;
                isOverwritten = true;
                return factor;
            }
            else
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                return CholeskySkyline.Factorize(NumColumns, valuesCopy, diagOffsets, tolerance);
            }
        }

        /// <summary>
        /// Calculate the LDL factorization. The matrix must be invertible, otherwise a <see cref="SingularMatrixException"/> 
        /// will be thrown. This method succeeds for all symmetric positive definite matrices, but not for all symmetric ones. 
        /// If <paramref name="inPlace"/> is set to true, this object must not be used again, otherwise a 
        /// <see cref="NullReferenceException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">
        /// False, to copy the internal non zero entries before factorization. True, to overwrite them with the factorized data, 
        /// thus saving memory and time. However, that will make this object unusable, so you MUST NOT call any other members 
        /// afterwards.
        /// </param>
        /// <param name="tolerance">
        /// If a diagonal entry is closer to zero than this tolerance, an <see cref="SingularMatrixException"/> exception will 
        /// be thrown.
        /// </param>
        /// <exception cref="SingularMatrixException">Thrown if the matrix is not invertible.</exception>
        /// <exception cref="NullReferenceException">
        /// Thrown if a member of his instance is accessed after this method is called.
        /// </exception>"
        public LdlSkyline FactorLdl(bool inPlace, double tolerance = LdlSkyline.PivotTolerance)
        {
            if (inPlace)
            {
                var factor = LdlSkyline.Factorize(NumColumns, values, diagOffsets, tolerance);
                // Set the skyline arrays to null to force NullReferenceException if they are accessed again.
                // TODO: perhaps there is a better way to handle this.
                values = null;
                diagOffsets = null;
                isOverwritten = true;
                return factor;
            }
            else
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                return LdlSkyline.Factorize(NumColumns, valuesCopy, diagOffsets, tolerance);
            }
        }

        /// <summary>
        /// Applies the Cholesky factorization to the independent columns of a symmetric positive semi-definite matrix,
        /// sets the dependent ones equal to columns of the identity matrix and return the nullspace of the matrix. Requires 
        /// extra memory for the basis vectors of the nullspace. If <paramref name="inPlace"/> is set to true, this object 
        /// must not be used again, otherwise a <see cref="NullReferenceException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">
        /// False, to copy the internal non zero entries before factorization. True, to overwrite them with the factorized data, 
        /// thus saving memory and time. However, that will make this object unusable, so you MUST NOT call any other members 
        /// afterwards.
        /// </param>
        /// <param name="pivotTolerance">
        /// If a diagonal entry is &lt;= <paramref name="pivotTolerance"/> it means that the corresponding column is dependent 
        /// on the rest. The Cholesky factorization only applies to independent column, while dependent ones are used to compute
        /// the nullspace. Therefore it is important to select a tolerance that will identify small pivots that result from 
        /// singularity, but not from ill-conditioning.
        /// </param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not positive definite.</exception>
        /// <exception cref="NullReferenceException">
        /// Thrown if a member of his instance is accessed after this method is called.
        /// </exception>"
        public SemidefiniteCholeskySkyline FactorSemidefiniteCholesky(bool inPlace, 
            double pivotTolerance = SemidefiniteCholeskySkyline.PivotTolerance)
        {
            if (inPlace)
            {
                var factor = SemidefiniteCholeskySkyline.Factorize(NumColumns, values, diagOffsets, pivotTolerance);
                // Set the skyline arrays to null to force NullReferenceException if they are accessed again.
                // TODO: perhaps there is a better way to handle this.
                values = null;
                diagOffsets = null;
                isOverwritten = true;
                return factor;
            }
            else
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                return SemidefiniteCholeskySkyline.Factorize(NumColumns, valuesCopy, diagOffsets, pivotTolerance);
            }
        }

        /// <summary>
        /// Applies the LDL factorization to the independent columns of a symmetric positive semi-definite matrix,
        /// sets the dependent ones equal to columns of the identity matrix and return the nullspace of the matrix. Requires 
        /// extra memory for the basis vectors of the nullspace. If <paramref name="inPlace"/> is set to true, this object 
        /// must not be used again, otherwise a <see cref="NullReferenceException"/> will be thrown.
        /// </summary>
        /// <param name="inPlace">
        /// False, to copy the internal non zero entries before factorization. True, to overwrite them with the factorized data, 
        /// thus saving memory and time. However, that will make this object unusable, so you MUST NOT call any other members 
        /// afterwards.
        /// </param>
        /// <param name="pivotTolerance">
        /// If a diagonal entry is &lt;= <paramref name="pivotTolerance"/> it means that the corresponding column is dependent 
        /// on the rest. The Cholesky factorization only applies to independent column, while dependent ones are used to compute
        /// the nullspace. Therefore it is important to select a tolerance that will identify small pivots that result from 
        /// singularity, but not from ill-conditioning.
        /// </param>
        /// <exception cref="IndefiniteMatrixException">Thrown if the matrix is not positive definite.</exception>
        /// <exception cref="NullReferenceException">
        /// Thrown if a member of his instance is accessed after this method is called.
        /// </exception>"
        public SemidefiniteLdlSkyline FactorSemidefiniteLdl(bool inPlace,
            double pivotTolerance = SemidefiniteLdlSkyline.PivotTolerance)
        {
            if (inPlace)
            {
                var factor = SemidefiniteLdlSkyline.Factorize(NumColumns, values, diagOffsets, pivotTolerance);
                // Set the skyline arrays to null to force NullReferenceException if they are accessed again.
                // TODO: perhaps there is a better way to handle this.
                values = null;
                diagOffsets = null;
                isOverwritten = true;
                return factor;
            }
            else
            {
                double[] valuesCopy = new double[values.Length];
                Array.Copy(values, valuesCopy, values.Length);
                return SemidefiniteLdlSkyline.Factorize(NumColumns, valuesCopy, diagOffsets, pivotTolerance);
            }
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetColumn(int)"/>.
        /// </summary>
        public Vector GetColumn(int colIndex)
        {
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            Preconditions.CheckIndexCol(this, colIndex);
            return Vector.CreateFromArray(SkylineSlicing.GetColumn(values, diagOffsets, colIndex));
        }

        /// <summary>
        /// Returns a <see cref="Vector"/> with the entries of the matrix's main diagonal.
        /// </summary>
        public Vector GetDiagonal() => Vector.CreateFromArray(GetDiagonalAsArray(), false);

        /// <summary>
        /// Returns an array with the entries of the matrix's main diagonal.
        /// </summary>
        public double[] GetDiagonalAsArray()
        {
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            return SkylineSlicing.GetDiagonal(values, diagOffsets);
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetRow(int)"/>.
        /// </summary>
        public Vector GetRow(int rowIndex) => GetColumn(rowIndex);

        /// <summary>
        /// See <see cref="ISparseMatrix.GetSparseFormat"/>.
        /// </summary>
        public SparseFormat GetSparseFormat()
        {
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            var format = new SparseFormat();
            format.RawValuesTitle = "Values";
            format.RawValuesArray = values;
            format.RawIndexArrays.Add("Diagonal offsets", diagOffsets);
            return format;
        }

        /// <summary>
        /// See <see cref="ISliceable2D.GetSubmatrix(int[], int[])"/>.
        /// </summary>
        public IMatrix GetSubmatrix(int[] rowIndices, int[] colIndices) => GetSubmatrixFull(rowIndices, colIndices);

        /// <summary>
        /// See <see cref="ISliceable2D.GetSubmatrix(int, int, int, int)"/>.
        /// </summary>
        public IMatrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
        {
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            int[] rowIndices = Enumerable.Range(rowStartInclusive, rowEndExclusive - rowStartInclusive).ToArray();
            int[] colIndices = Enumerable.Range(colStartInclusive, colEndExclusive - colStartInclusive).ToArray();
            return GetSubmatrix(rowIndices, colIndices);
        }

        public CscMatrix GetSubmatrixCsc(int[] rowIndices, int[] colIndices)
        {
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            return SkylineSlicing.GetSubmatrixCsc(values, diagOffsets, rowIndices, colIndices);
        }

        public Matrix GetSubmatrixFull(int[] rowIndices, int[] colIndices)
        {
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            return DenseStrategies.GetSubmatrix(this, rowIndices, colIndices);
        }

        public Matrix GetSubmatrixSymmetricFull(int[] indices)
        {
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            return SkylineSlicing.GetSubmatrixSymmetricFull(values, diagOffsets, indices);
            //// I am not sure that the above is faster than: 
            //return DenseStrategies.GetSubmatrix(this, indices, indices);
        }

        public SymmetricMatrix GetSubmatrixSymmetricPacked(int[] indices)
        {
            //TODO: perhaps this can be combined with the CSC and full version to get all 2 submatrices needed for 
            //      Schur complements more efficiently.
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            return SkylineSlicing.GetSubmatrixSymmetricPacked(values, diagOffsets, indices);
        }

        public SparsityPatternSymmetric GetSubmatrixSymmetricPattern(int[] indices)
        {
            //TODO: perhaps this can be combined with the CSC and full version to get all 2 submatrices needed for 
            //      Schur complements more efficiently.
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            return SkylineSlicing.GetSubmatrixSymmetricPattern(values, diagOffsets, indices);
        }

        public SkylineMatrix GetSubmatrixSymmetricSkyline(int[] indices) 
        {
            //TODO: perhaps this can be combined with the CSC and full version to get all 2 submatrices needed for 
            //      Schur complements more efficiently.
            if (isOverwritten) throw new MatrixDataOverwrittenException();
            return SkylineSlicing.GetSubmatrixSymmetricSkyline(values, diagOffsets, indices);
        }
        
        /// <summary>
        /// See <see cref="IMatrixView.LinearCombination(double, IMatrixView, double)"/>.
        /// </summary>
        public IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is SkylineMatrix otherSKY) // In case both matrices have the exact same index arrays
            {
                if (HasSameIndexer(otherSKY))
                {
                    // Do not copy the index arrays, since they are already spread around. TODO: is this a good idea?
                    double[] resultValues = new double[values.Length];
                    if (thisCoefficient == 1.0)
                    {
                        Array.Copy(this.values, resultValues, values.Length);
                        Blas.Daxpy(values.Length, otherCoefficient, otherSKY.values, 0, 1, this.values, 0, 1);
                    }
                    else if (otherCoefficient == 1.0)
                    {
                        Array.Copy(otherSKY.values, resultValues, values.Length);
                        Blas.Daxpy(values.Length, thisCoefficient, this.values, 0, 1, resultValues, 0, 1);
                    }
                    else
                    {
                        Array.Copy(this.values, resultValues, values.Length);
                        BlasExtensions.Daxpby(values.Length, otherCoefficient, otherSKY.values, 0, 1,
                            thisCoefficient, resultValues, 0, 1);
                    }
                    return new SkylineMatrix(NumColumns, resultValues, this.diagOffsets);
                }
            }

            // All entries must be processed. TODO: optimizations may be possible (e.g. only access the nnz in this matrix)
            return DenseStrategies.LinearCombination(this, thisCoefficient, otherMatrix, otherCoefficient);
        }

        /// <summary>
        /// See <see cref="IMatrix.LinearCombinationIntoThis(double, IMatrixView, double)"/>.
        /// </summary>
        public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
        {
            if (otherMatrix is SkylineMatrix casted) LinearCombinationIntoThis(thisCoefficient, casted, otherCoefficient);
            else throw new SparsityPatternModifiedException(
                "This operation is legal only if the other matrix has the same sparsity pattern");
        }

        /// <summary>
        /// Performs the following operation for th non-zero entries (i, j), such that 0 &lt;= j &lt; <see cref="NumColumns"/> 
        /// 0 &lt;= i &lt;= j: this[i, j] = <paramref name="thisCoefficient"/> * this[i, j] 
        ///     + <paramref name="otherCoefficient"/> * <paramref name="otherMatrix"/>[i, j].  
        /// The resulting matrix overwrites the entries of this <see cref="CscMatrix"/> instance.
        /// </summary>
        /// <param name="thisCoefficient">A scalar that multiplies each entry of this <see cref="Matrix"/>.</param>
        /// <param name="otherMatrix">A matrix for whose active columns (non-zero part of the column) are shorter or equal to the
        ///     active columns of this <see cref="SkylineMatrix"/> instance.</param>
        /// <param name="otherCoefficient">A scalar that multiplies each entry of <paramref name="otherMatrix"/>.</param>
        /// <exception cref="SparsityPatternModifiedException">Thrown if <paramref name="otherMatrix"/> has taller active columns 
        ///     than this <see cref="SkylineMatrix"/> instance.</exception>
        public void LinearCombinationIntoThis(double thisCoefficient, SkylineMatrix otherMatrix, double otherCoefficient)
        {
            if (HasSameIndexer(otherMatrix)) // no need to check dimensions if the indexing arrays are the same
            {
                if (thisCoefficient == 1.0)
                {
                    Blas.Daxpy(values.Length, otherCoefficient, otherMatrix.values, 0, 1, this.values, 0, 1);
                }
                else
                {
                    BlasExtensions.Daxpby(values.Length, otherCoefficient, otherMatrix.values, 0, 1,
                        thisCoefficient, this.values, 0, 1);
                }
            }
            else
            {
                Preconditions.CheckSameMatrixDimensions(this, otherMatrix);
                for (int j = 0; j < NumColumns; ++j)
                {
                    // Check if the current column of this matrix is tall enough
                    int thisDiagOffset = this.diagOffsets[j];
                    int otherDiagOffset = otherMatrix.diagOffsets[j];
                    int thisColTop = j - this.diagOffsets[j + 1] + thisDiagOffset + 1;
                    int otherColTop = j - this.diagOffsets[j + 1] + otherDiagOffset + 1;
                    if (thisColTop < otherColTop) throw new SparsityPatternModifiedException(
                        $"Column {j} of this matrix is shorter, which would result in overflow");

                    // Do the operation between the two columns. The column of this matrix is taller or equal to the other one.
                    for (int i = j; i >= otherColTop; --i) // non zero entries of shortest column, including diagonal
                    {
                        int thisIndex = thisDiagOffset + j - i;
                        this.values[thisIndex] = thisCoefficient * this.values[thisIndex] 
                            + otherCoefficient * otherMatrix.values[otherDiagOffset + j - i];
                    }
                    for (int i = otherColTop - 1; i >= thisColTop; --i) // non zero entries of the above the shortest column
                    {
                        int thisIndex = thisDiagOffset + j - i;
                        this.values[thisIndex] = thisCoefficient * this.values[thisIndex];
                    }
                }
            }            
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyLeft(IMatrixView, bool, bool)"/>.
        /// </summary>
        /// <remarks>
        /// <paramref name="transposeThis"/> does not affect the result, as a <see cref="SkylineMatrix"/> is symmetric.
        /// </remarks>
        public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return DenseStrategies.Multiply(other, this, transposeOther, transposeThis);
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyRight(IMatrixView, bool, bool)"/>.
        /// </summary>
        /// <remarks>
        /// <paramref name="transposeThis"/> does not affect the result, as a <see cref="SkylineMatrix"/> is symmetric.
        /// </remarks>
        public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
        {
            return DenseStrategies.Multiply(this, other, transposeThis, transposeOther);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Multiply(IVectorView, bool)"/>.
        /// </summary>
        /// <remarks>
        /// <paramref name="transposeThis"/> does not affect the result, as a <see cref="SkylineMatrix"/> is symmetric.
        /// </remarks>
        public IVector Multiply(IVectorView vector, bool transposeThis = false)
        {
            if (vector is Vector casted) return Multiply(casted);
            else throw new NotImplementedException();
        }

        /// <summary>
        /// Performs the matrix-vector multiplication: this * <paramref name="vector"/> = <paramref name="vector"/> * this.
        /// </summary>
        /// <param name="vector">A vector with <see cref="IIndexable1D.Length"/> being equal to 
        ///     this.<see cref="NumColumns"/>.</param>
        /// <exception cref="NonMatchingDimensionsException">Thrown if the <see cref="IIndexable1D.Length"/> of
        ///     <paramref name="vector"/> is different than the <see cref="NumColumns"/> of this.</exception>
        public Vector Multiply(Vector vector)
        {
            //TODO: this performs redundant dimension checks
            var result = Vector.CreateZero(NumColumns);
            MultiplyIntoResult(vector, result);
            return result;
        }

        /// <summary>
        /// See <see cref="IMatrixView.MultiplyIntoResult(IVectorView, IVector, bool)"/>.
        /// </summary>
        public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
        {
            if ((lhsVector is Vector lhsDense) && (rhsVector is Vector rhsDense))
            {
                MultiplyIntoResult(lhsDense, rhsDense);
            }
            else throw new NotImplementedException();
        }

        /// <summary>
        /// Performs the matrix-vector multiplication: <paramref name="rhsVector"/> = this * <paramref name="vector"/>.
        /// To multiply this * columnVector, set <paramref name="transposeThis"/> to false.
        /// To multiply rowVector * this, set <paramref name="transposeThis"/> to true.
        /// The resulting vector will overwrite the entries of <paramref name="rhsVector"/>.
        /// </summary>
        /// <param name="lhsVector">
        /// The vector that will be multiplied by this matrix. It sits on the left hand side of the equation y = A * x.
        /// Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
        /// == this.<see cref="IIndexable2D.NumColumns"/>.
        /// </param>
        /// <param name="rhsVector">
        /// The vector that will be overwritten by the result of the multiplication. It sits on the right hand side of the 
        /// equation y = A * x. Constraints: <paramref name="lhsVector"/>.<see cref="IIndexable1D.Length"/> 
        /// == this.<see cref="IIndexable2D.NumRows"/>.
        /// </param>
        /// <exception cref="NonMatchingDimensionsException">
        /// Thrown if the <see cref="IIndexable1D.Length"/> of <paramref name="lhsVector"/> or <paramref name="rhsVector"/> 
        /// violate the described contraints.
        /// </exception>
        public void MultiplyIntoResult(Vector lhsVector, Vector rhsVector)
        {
            Preconditions.CheckMultiplicationDimensions(NumColumns, lhsVector.Length);
            Preconditions.CheckSystemSolutionDimensions(NumRows, rhsVector.Length);
            ManagedSparseBlasProvider.UniqueInstance.Dskymv(
                NumColumns, values, diagOffsets, lhsVector.RawData, rhsVector.RawData);
        }

        /// <summary>
        /// See <see cref="IReducible.Reduce(double, ProcessEntry, ProcessZeros, Reduction.Finalize)"/>.
        /// </summary>
        public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
        {
            double aggregator = identityValue;
            int nnz = values.Length;
            for (int i = 0; i < nnz; ++i) aggregator = processEntry(values[i], aggregator);
            aggregator = processZeros(NumColumns * NumColumns - nnz, aggregator);
            return finalize(aggregator);
        }

        /// <summary>
        /// See <see cref="IMatrixView.Scale(double)"/>.
        /// </summary>
        IMatrix IMatrixView.Scale(double scalar) => Scale(scalar);

        /// <summary>
        /// Performs the following operation for the non-zero entries (i, j), such that 0 &lt;= j &lt; <see cref="NumColumns"/>,
        /// 0 &lt;= i &lt;= j: result[i, j] = <paramref name="scalar"/> * this[i, j].
        /// The resulting matrix is written to a new <see cref="SkylineMatrix"/> and then returned.
        /// </summary>
        /// <param name="scalar">A scalar that multiplies each entry of this matrix.</param>
        public SkylineMatrix Scale(double scalar)
        {
            int nnz = this.values.Length;
            double[] resultValues = new double[nnz];
            Array.Copy(this.values, resultValues, nnz); //TODO: perhaps I should also copy the indexers
            Blas.Dscal(nnz, scalar, resultValues, 0, 1);
            return new SkylineMatrix(this.NumColumns, resultValues, this.diagOffsets);
        }

        /// <summary>
        /// See <see cref="IMatrix.ScaleIntoThis(double)"/>.
        /// </summary>
        public void ScaleIntoThis(double scalar) => Blas.Dscal(values.Length, scalar, values, 0, 1);

        /// <summary>
        /// See <see cref="IMatrix.SetEntryRespectingPattern(int, int, double)"/>.
        /// </summary>
        public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value) => this[rowIdx, colIdx] = value;

        /// <summary>
        /// See <see cref="IMatrixView.Transpose"/>.
        /// </summary>
        public IMatrix Transpose() => Copy(false);

        /// <summary>
        /// Perhaps this should be manually inlined. Testing needed.
        /// </summary>
        /// <param name="rowIdx"></param>
        /// <param name="colIdx"></param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        internal static void ProcessIndices(ref int rowIdx, ref int colIdx)
        {
            if (rowIdx > colIdx)
            {
                int swap = rowIdx;
                rowIdx = colIdx;
                colIdx = swap;
            }
        }

        private bool HasSameIndexer(SkylineMatrix other) => this.diagOffsets == other.diagOffsets;
    }
}
