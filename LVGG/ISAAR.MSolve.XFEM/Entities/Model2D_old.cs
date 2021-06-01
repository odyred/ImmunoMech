﻿using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Enrichments.Items;
using ISAAR.MSolve.XFEM.FreedomDegrees.Ordering;

namespace ISAAR.MSolve.XFEM.Entities
{
    /// <summary>
    /// This class both manages and assembles the FEM entities. TODO: Those 2 should be split. I could have a nested 
    /// Builder or better yet a UI class for assembling the model.
    /// </summary>
    public class Model2D_old
    {
        private readonly List<XNode> nodes;
        private readonly List<XContinuumElement2D> elements;
        private readonly List<IEnrichmentItem2D> enrichments;
        private readonly Table<XNode, StructuralDof, double> constraints;

        /// <summary>
        /// Multiple loads for the same dof are not allowed. Any attempt to input them will result in an exception 
        /// thrown by the table object.
        /// </summary>
        private readonly Table<XNode, StructuralDof, double> loads;

        public IReadOnlyList<XNode> Nodes { get { return nodes; } }
        public IReadOnlyList<XContinuumElement2D> Elements { get { return elements; } }
        public IReadOnlyList<IEnrichmentItem2D> Enrichments { get { return enrichments; } }
        public ITable<XNode, StructuralDof, double> Constraints { get { return constraints; } }

        public Model2D_old()
        {
            this.nodes = new List<XNode>();
            this.elements = new List<XContinuumElement2D>();
            this.enrichments = new List<IEnrichmentItem2D>();
            this.constraints = new Table<XNode, StructuralDof, double>();
            this.loads = new Table<XNode, StructuralDof, double>();
        }

        public void AddNode(XNode node)
        {
            if (nodes.Contains(node))
                throw new ArgumentException("There is already a node with id = " + node.ID);
            nodes.Add(node);
        }

        public void AddElement(XContinuumElement2D element)
        {
            if (elements.Contains(element))
                throw new ArgumentException("This element is already inserted");
            elements.Add(element);
        }

        public void AddEnrichment(IEnrichmentItem2D enrichment)
        {
            if (enrichments.Contains(enrichment))
                throw new ArgumentException("This enrichment is already inserted");
            enrichments.Add(enrichment);
        }

        public void AddConstraint(XNode node, StructuralDof dof, double displacement)
        {
            // TODO: This should be done more efficiently than O(N)
            if (!nodes.Contains(node)) throw new ArgumentException("There is no such node");
            constraints[node, dof] = displacement;
        }

        //TODO: Should I use the node's id instead? In a UI class, I probably should.
        public void AddNodalLoad(XNode node, StructuralDof dof, double magnitude)
        {
            // TODO: This should be done more efficiently than O(N)
            if (!nodes.Contains(node)) throw new ArgumentException("There is no such node");
            loads[node, dof] = magnitude;
        }

        public Vector CalculateFreeForces(IDofOrderer dofOrderer)
        {
            double[] rhs = new double[dofOrderer.NumStandardDofs + dofOrderer.NumEnrichedDofs];
            foreach ((XNode node, StructuralDof dofType, double load) in loads)
            {
                try
                {
                    int dof = dofOrderer.GetStandardDofOf(node, dofType);
                    rhs[dof] += load; // This supports multiple loads on the same dof, which isn't implemented yet
                }
                catch (KeyNotFoundException ex)
                {
                    throw new NotImplementedException("Load on a constrained dof at node "
                    + node.ID + ", axis " + dofType, ex);
                }
            }
            return Vector.CreateFromArray(rhs);
        }

        public Vector CalculateStandardForces(IDofOrderer dofOrderer)
        {
            double[] rhs = new double[dofOrderer.NumStandardDofs];
            foreach ((XNode node, StructuralDof dofType, double load) in loads)
            {
                try
                {
                    int dof = dofOrderer.GetStandardDofOf(node, dofType);
                    rhs[dof] += load; // This supports multiple loads on the same dof, which isn't implemented yet
                }
                catch (KeyNotFoundException ex)
                {
                    throw new NotImplementedException("Load on a constrained dof at node "
                    + node.ID + ", axis " + dofType, ex);
                }
            }
            return Vector.CreateFromArray(rhs);
        }

        public Vector CalculateConstrainedDisplacements(IDofOrderer dofOrderer)
        {
            double[] uc = new double[dofOrderer.NumConstrainedDofs];
            foreach ((XNode node, StructuralDof dofType, double displacement) in constraints)
            {
                int dof = dofOrderer.GetConstrainedDofOf(node, dofType);
                uc[dof] = displacement;
            }
            return Vector.CreateFromArray(uc);
        }
    }
}
