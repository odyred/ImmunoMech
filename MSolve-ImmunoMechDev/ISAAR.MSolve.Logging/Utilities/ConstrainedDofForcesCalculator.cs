﻿using System.Collections.Generic;
using System.Diagnostics;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: finding the contributing elements and the corresponding local dof indices can be done only once in the constructor.
namespace ISAAR.MSolve.Logging.Utilities
{
    /// <summary>
    /// This does not work if the requested node belongs to an element that contains embedded elements.
    /// </summary>
    internal class ConstrainedDofForcesCalculator
    {
        private readonly Subdomain subdomain;

        internal ConstrainedDofForcesCalculator(Subdomain subdomain)
        {
            this.subdomain = subdomain;
        }

        internal double CalculateForceAt(Node node, IDofType dofType, IVectorView totalDisplacements)
        {
            double totalForce = 0.0;

            foreach (Element element in node.ElementsDictionary.Values)
            {
                // It is possible that one of the elements at this node does not engage this dof type, in which case -1 will be returned.
                // We will not have any contribution from them. If none of the elements engage this dof type, the total force will always be 0.
                int monitorDofIdx = FindLocalDofIndex(element, node, dofType);
                if (monitorDofIdx == -1) continue;

                //TODO: if an element has embedded elements, then we must also take into account their forces.
                double[] totalElementDisplacements = subdomain.CalculateElementDisplacements(element, totalDisplacements);
                double[] elementForces = element.ElementType.CalculateForcesForLogging(element, totalElementDisplacements);

                totalForce += elementForces[monitorDofIdx];
            }

            return totalForce;
        }

        /// <summary>
        /// Returns -1 if the element does not engage the requested <see cref="IDofType"/>
        /// </summary>
        private int FindLocalDofIndex(Element element, Node node, IDofType dofType)
        {
            int localNodeIdx = element.Nodes.IndexOf(node);
            Debug.Assert(localNodeIdx != -1, "The element does not contain this node.");
            IReadOnlyList<IReadOnlyList<IDofType>> elementDofs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            int localDofIdx = elementDofs[localNodeIdx].FindFirstIndex(dofType);
            int multNum = elementDofs[localNodeIdx].Count;
            int dofIdx = multNum * (localNodeIdx + 1) - (localDofIdx + 1);
            return dofIdx;
        }
    }
}
