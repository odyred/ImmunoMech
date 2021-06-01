using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices.Builders;
using ISAAR.MSolve.Solvers.Commons;

//TODO: Instead of storing the raw CSC arrays, use a reusable DOK or CscIndexer class. That class should provide methods to 
//      assemble the values part of the global matrix more efficiently than the general purpose DOK. The general purpose DOK 
//      should only be used to assemble the first global matrix and whenever the dof ordering changes. Now it is used everytime 
//      and the indexing arrays are discarded.
namespace ISAAR.MSolve.Solvers.Assemblers
{
    /// <summary>
    /// Builds the global matrix of the linear system that will be solved. This matrix is square and stored in CSC format, but
    /// both triangles are explicitly stored. This format is suitable for matrix/vector multiplications, therefore it can be 
    /// combined with many iterative solvers. 
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class CscAssembler : IGlobalMatrixAssembler<CscMatrix>
    {
        private const string name = "CscAssembler"; // for error messages
        private readonly bool sortRowsOfEachCol;
        private ConstrainedMatricesAssembler constrainedAssembler = new ConstrainedMatricesAssembler();

        bool isIndexerCached = false;
        private int[] cachedRowIndices, cachedColOffsets;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="sortRowsOfEachCol">
        /// Sorting the rows of each column in the CSC storage format may increase performance of the matrix vector 
        /// multiplications. It is recommended to set it to true, especially for iterative linear system solvers.
        /// </param>
        public CscAssembler(bool sortRowsOfEachCol = true)
        {
            this.sortRowsOfEachCol = sortRowsOfEachCol;
        }

        public CscMatrix BuildGlobalMatrix(ISubdomainFreeDofOrdering dofOrdering, IEnumerable<IElement> elements, 
            IElementMatrixProvider matrixProvider)
        {
            int numFreeDofs = dofOrdering.NumFreeDofs;
            var subdomainMatrix = DokColMajor.CreateEmpty(numFreeDofs, numFreeDofs);

            foreach (IElement element in elements)
            {
                (int[] elementDofIndices, int[] subdomainDofIndices) = dofOrdering.MapFreeDofsElementToSubdomain(element);
                IMatrix elementMatrix = matrixProvider.Matrix(element);
                subdomainMatrix.AddSubmatrix(elementMatrix, elementDofIndices, subdomainDofIndices, elementDofIndices, subdomainDofIndices);
            }

            (double[] values, int[] rowIndices, int[] colOffsets) = subdomainMatrix.BuildCscArrays(sortRowsOfEachCol);
            if (!isIndexerCached)
            {
                cachedRowIndices = rowIndices;
                cachedColOffsets = colOffsets;
                isIndexerCached = true;
            }
            else
            {
                Debug.Assert(Utilities.AreEqual(cachedRowIndices, rowIndices));
                Debug.Assert(Utilities.AreEqual(cachedColOffsets, colOffsets));
            }
            return CscMatrix.CreateFromArrays(numFreeDofs, numFreeDofs, values, cachedRowIndices, cachedColOffsets, false);
        }

        public (CscMatrix matrixFreeFree, IMatrixView matrixFreeConstr, IMatrixView matrixConstrFree,
            IMatrixView matrixConstrConstr) BuildGlobalSubmatrices(
            ISubdomainFreeDofOrdering freeDofOrdering, ISubdomainConstrainedDofOrdering constrainedDofOrdering,
            IEnumerable<IElement> elements, IElementMatrixProvider matrixProvider)
        {
            int numFreeDofs = freeDofOrdering.NumFreeDofs;
            var subdomainMatrix = DokColMajor.CreateEmpty(numFreeDofs, numFreeDofs);

            //TODO: also reuse the indexers of the constrained matrices.
            constrainedAssembler.InitializeNewMatrices(freeDofOrdering.NumFreeDofs, constrainedDofOrdering.NumConstrainedDofs);

            // Process the stiffness of each element
            foreach (IElement element in elements)
            {
                (int[] elementDofsFree, int[] subdomainDofsFree) = freeDofOrdering.MapFreeDofsElementToSubdomain(element);
                (int[] elementDofsConstrained, int[] subdomainDofsConstrained) =
                    constrainedDofOrdering.MapConstrainedDofsElementToSubdomain(element);

                IMatrix elementMatrix = matrixProvider.Matrix(element);
                subdomainMatrix.AddSubmatrixSymmetric(elementMatrix, elementDofsFree, subdomainDofsFree);
                constrainedAssembler.AddElementMatrix(elementMatrix, elementDofsFree, subdomainDofsFree,
                    elementDofsConstrained, subdomainDofsConstrained);
            }

            // Create and cache the CSC arrays for the free dofs.
            (double[] values, int[] rowIndices, int[] colOffsets) = subdomainMatrix.BuildCscArrays(sortRowsOfEachCol);
            if (!isIndexerCached)
            {
                cachedRowIndices = rowIndices;
                cachedColOffsets = colOffsets;
                isIndexerCached = true;
            }
            else
            {
                Debug.Assert(Utilities.AreEqual(cachedRowIndices, rowIndices));
                Debug.Assert(Utilities.AreEqual(cachedColOffsets, colOffsets));
            }

            // Create the free and constrained matrices. 
            subdomainMatrix = null; // Let the DOK be garbaged collected early, in case there isn't sufficient memory.
            var matrixFreeFree = 
                CscMatrix.CreateFromArrays(numFreeDofs, numFreeDofs, values, cachedRowIndices, cachedColOffsets, false);
            (CsrMatrix matrixConstrFree, CsrMatrix matrixConstrConstr) = constrainedAssembler.BuildMatrices();
            return (matrixFreeFree, matrixConstrFree, matrixConstrFree.TransposeToCSC(false), matrixConstrConstr);
        }

        public void HandleDofOrderingWillBeModified()
        {
            //TODO: perhaps the indexer should be disposed altogether. Then again it could be in use by other matrices.
            cachedRowIndices = null;
            cachedColOffsets = null;
            isIndexerCached = false;
        }
    }
}
