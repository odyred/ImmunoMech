﻿using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Reordering;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    public class SubdomainFreeRowDofOrderingGeneral : ISubdomainFreeDofOrdering
    {
        public SubdomainFreeRowDofOrderingGeneral(int numFreeDofs, DofTable subdomainFreeDofs)
        {
            this.NumFreeDofs = numFreeDofs;
            this.FreeDofs = subdomainFreeDofs;
        }

        public DofTable FreeDofs { get; }
        public int NumFreeDofs { get; }

        public void AddVectorElementToSubdomain(IElement element, double[] elementVector, IVector subdomainVector)
        {
            IReadOnlyList<INode> elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            IReadOnlyList<IReadOnlyList<IDofType>> elementDofs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);

            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx],
                        out int subdomainDofIdx);
                    if (isFree)
                    {
                        subdomainVector.Set(subdomainDofIdx, subdomainVector[subdomainDofIdx] + elementVector[elementDofIdx]);
                    }
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
        }

        public int CountElementDofs(IElement element)
        {
            IReadOnlyList<INode> elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            IReadOnlyList<IReadOnlyList<IDofType>> elementDofs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            int numElementDofs = 0;
            for (int nodeIdx = 0; nodeIdx < elementDofs.Count; ++nodeIdx) numElementDofs += elementDofs[nodeIdx].Count;
            return numElementDofs;
        }

        public double[] ExtractVectorElementFromSubdomain(IElement element, IVectorView subdomainVector)
        {
            IReadOnlyList<INode> elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            IReadOnlyList<IReadOnlyList<IDofType>> elementDofs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);

            int numElementDofs = 0;
            for (int nodeIdx = 0; nodeIdx < elementDofs.Count; ++nodeIdx) numElementDofs += elementDofs[nodeIdx].Count;
            var elementVector = new double[numElementDofs];

            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx],
                        out int subdomainDofIdx);
                    if (isFree) elementVector[elementDofIdx] = subdomainVector[subdomainDofIdx];
                    // Else, the quantity of interest is 0.0 at all constrained dofs.
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }

            return elementVector;
        }

        public void ExtractVectorElementFromSubdomain(IElement element, IVectorView subdomainVector, IVector elementVector)
        {
            IReadOnlyList<INode> elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            IReadOnlyList<IReadOnlyList<IDofType>> elementDofs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            //int numElementDofs = 0;
            //for (int nodeIdx = 0; nodeIdx < elementDofs.Count; ++nodeIdx) numElementDofs += elementDofs[nodeIdx].Count;

            int elementDofIdx = 0;
            for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
            {
                for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx],
                        out int subdomainDofIdx);
                    if (isFree) elementVector.Set(elementDofIdx, subdomainVector[subdomainDofIdx]);
                    // Else, the quantity of interest is 0.0 at all constrained dofs.
                    ++elementDofIdx; // This must be incremented for constrained dofs as well
                }
            }
        }

        public (int[] elementDofIndices, int[] subdomainDofIndices) MapFreeDofsElementToSubdomain(IElement element)
        {
            var collocationElement = ((ICollocationElement)element);

            //IList<INode> elementNodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            //IList<IList<DOFType>> elementDofs = element.ElementType.DofEnumerator.GetDOFTypes(element);

            // Count the dof superset (free and constrained) to allocate enough memory and avoid resizing
            var rowDofs = collocationElement.GetDOFTypesForDOFEnumeration(element);
            int allElementDofs = rowDofs.Count;
            var elementDofIndices = new List<int>(allElementDofs);
            var subdomainDofIndices = new List<int>(allElementDofs);
            var elementDofs = new IDofType[1][];
            elementDofs[0] = rowDofs.ToArray();

            int elementDofIdx = 0;
                for (int dofIdx = 0; dofIdx < elementDofs[0].Length; ++dofIdx)
                {
                    bool isFree = FreeDofs.TryGetValue(collocationElement.CollocationPoint, elementDofs[0][dofIdx],
                        out int subdomainDofIdx);
                    if (isFree)
                    {
                        elementDofIndices.Add(elementDofIdx);
                        subdomainDofIndices.Add(subdomainDofIdx);
                    }
                    ++elementDofIdx; 
                }
            
            return (elementDofIndices.ToArray(), subdomainDofIndices.ToArray());
        }

        public void Reorder(IReorderingAlgorithm reorderingAlgorithm, ISubdomain subdomain)
        {
            var pattern = SparsityPatternSymmetric.CreateEmpty(NumFreeDofs);
            foreach (var element in subdomain.Elements)
            {
                (int[] elementDofIndices, int[] subdomainDofIndices) = MapFreeDofsElementToSubdomain(element);

                //TODO: ISubdomainFreeDofOrdering could perhaps return whether the subdomainDofIndices are sorted or not.
                pattern.ConnectIndices(subdomainDofIndices, false);
            }
            (int[] permutation, bool oldToNew) = reorderingAlgorithm.FindPermutation(pattern);
            FreeDofs.Reorder(permutation, oldToNew);
        }

        public void ReorderNodeMajor(IReadOnlyList<INode> sortedNodes) => FreeDofs.ReorderNodeMajor(sortedNodes);

        //public IReadOnlyDictionary<int, int> MapFreeDofsElementToSubdomain(IElement element)
        //{
        //    IList<INode> elementNodes = element.IElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
        //    IList<IList<DOFType>> elementDofs = element.IElementType.DOFEnumerator.GetDOFTypes(element);

        //    var dofMap = new Dictionary<int, int>();
        //    int elementDofIdx = 0;
        //    for (int nodeIdx = 0; nodeIdx < elementNodes.Count; ++nodeIdx)
        //    {
        //        for (int dofIdx = 0; dofIdx < elementDofs[nodeIdx].Count; ++dofIdx)
        //        {
        //            bool isFree = FreeDofs.TryGetValue(elementNodes[nodeIdx], elementDofs[nodeIdx][dofIdx],
        //                out int subdomainDofIdx);
        //            if (isFree) dofMap[elementDofIdx] = subdomainDofIdx;
        //            ++elementDofIdx; // This must be incremented for constrained dofs as well
        //        }
        //    }
        //    return dofMap;
        //}
    }
}
