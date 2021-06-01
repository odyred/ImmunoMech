﻿using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Output
{
    class StressRecovery
    {
        private readonly Model2D_old model;
        private readonly IDofOrderer dofOrderer;

        public StressRecovery(Model2D_old model, IDofOrderer dofOrderer)
        {
            this.model = model;
            this.dofOrderer = dofOrderer;
        }

        public IReadOnlyList<Tensor2D> ComputeSmoothedNodalStresses(Vector solution)
        {
            var stressesFromAllElements = new Dictionary<XNode, List<Tensor2D>>();
            foreach (var node in model.Nodes) stressesFromAllElements[node] = new List<Tensor2D>();
            Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(dofOrderer);

            foreach (var element in model.Elements)
            {
                IReadOnlyDictionary<XNode, Tensor2D> elementStresses =
                    ComputeNodalStressesOfElement(element, solution, constrainedDisplacements);
                foreach (var nodeStressPair in elementStresses)
                {
                    stressesFromAllElements[nodeStressPair.Key].Add(nodeStressPair.Value);
                }
            }

            // Average with equal weights for all elements. TODO: perhaps vary the weights depending on the element type/area
            var nodalStresses = new Tensor2D[model.Nodes.Count];
            for (int i = 0; i < model.Nodes.Count; ++i)
            {
                XNode node = model.Nodes[i];
                double stressXX = 0.0, stressYY = 0.0, stressXY = 0.0;
                foreach (var tensor in stressesFromAllElements[node])
                {
                    stressXX += tensor.XX;
                    stressYY += tensor.YY;
                    stressXY += tensor.XY;
                }
                int contributingElementsCount = stressesFromAllElements[node].Count;
                stressXX /= contributingElementsCount;
                stressYY /= contributingElementsCount;
                stressXY /= contributingElementsCount;
                nodalStresses[i] = new Tensor2D(stressXX, stressYY, stressXY);
            }

            return nodalStresses;
        }

        public IReadOnlyDictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>> ComputeElementWiseStresses(
            Vector solution)
        {
            Vector constrainedDisplacements = model.CalculateConstrainedDisplacements(dofOrderer);
            var allStresses = new Dictionary<XContinuumElement2D, IReadOnlyList<Tensor2D>>();
            foreach (var element in model.Elements)
            {
                IReadOnlyDictionary<XNode, Tensor2D> elementStresses =
                    ComputeNodalStressesOfElement(element, solution, constrainedDisplacements);
                allStresses[element] = elementStresses.Values.ToArray();
            }
            return allStresses;
        }

        // Computes stresses directly at the nodes. The other approach is to compute them at Gauss points and then extrapolate
        private IReadOnlyDictionary<XNode, Tensor2D> ComputeNodalStressesOfElement(XContinuumElement2D element,
            Vector freeDisplacements, Vector constrainedDisplacements)
        {
            Vector standardDisplacements = dofOrderer.ExtractDisplacementVectorOfElementFromGlobal(element, 
                freeDisplacements, constrainedDisplacements);
            Vector enrichedDisplacements = 
                dofOrderer.ExtractEnrichedDisplacementsOfElementFromGlobal(element, freeDisplacements);

            IReadOnlyList<NaturalPoint> naturalNodes = element.Interpolation.NodalNaturalCoordinates;
            var nodalStresses = new Dictionary<XNode, Tensor2D>();
            for (int i = 0; i < element.Nodes.Count; ++i)
            {
                EvalInterpolation2D evaluatedInterpolation =
                    element.Interpolation.EvaluateAllAt(element.Nodes, naturalNodes[i]);
                Matrix2by2 displacementGradient = element.CalculateDisplacementFieldGradient(
                    naturalNodes[i], evaluatedInterpolation, standardDisplacements, enrichedDisplacements);
                Matrix constitutive = 
                    element.Material.CalculateConstitutiveMatrixAt(naturalNodes[i], evaluatedInterpolation);
                nodalStresses[element.Nodes[i]] = element.CalculateStressTensor(displacementGradient, constitutive);
            }

            return nodalStresses;
        }
    }
}
