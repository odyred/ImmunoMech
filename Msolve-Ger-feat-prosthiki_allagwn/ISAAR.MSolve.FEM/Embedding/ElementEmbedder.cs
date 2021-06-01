﻿using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM.Embedding
{
    public class SuperElementDOF
    {
        public Element Element { get; set; }
        public Node HostNode { get; set; }
        public Node EmbeddedNode { get; set; }
        public DOFType DOF { get; set; }

        public override bool Equals(object obj)
        {
            if (obj is SuperElementDOF == false) return false;
            var e = obj as SuperElementDOF;
            if (e == null) return false;

            if ((e.Element == null && this.Element != null) || (e.Element != null && this.Element == null)) return false;
            if ((e.HostNode == null && this.HostNode != null) || (e.HostNode != null && this.HostNode == null)) return false;
            if ((e.EmbeddedNode == null && this.EmbeddedNode != null) || (e.EmbeddedNode != null && this.EmbeddedNode == null)) return false;

            return (e.DOF == this.DOF && 
                ((e.Element == null && this.Element == null) || (e.Element.ID == this.Element.ID)) && 
                ((e.HostNode == null && this.HostNode == null) || (e.HostNode.ID == this.HostNode.ID)) &&
                ((e.EmbeddedNode == null && this.EmbeddedNode == null) || (e.EmbeddedNode.ID == this.EmbeddedNode.ID)));
        }

        public override int GetHashCode()
        {
            int elementID = Element == null ? 0 : Element.ID;
            int hostNodeID = HostNode == null ? 0 : HostNode.ID;
            return String.Format("H{3}E{0}N{1}{2}", elementID, hostNodeID, DOF, EmbeddedNode.ID).GetHashCode();
        }

        public override string ToString()
        {
            int elementID = Element == null ? 0 : Element.ID;
            int hostNodeID = HostNode == null ? 0 : HostNode.ID;
            return String.Format("N:{3} -> E:{0}, N:{1}, {2}", elementID, hostNodeID, DOF, EmbeddedNode.ID);
        }
    }

    public class ElementEmbedder : IFiniteElementDOFEnumerator
    {
        private readonly Model model;
        private readonly Element embeddedElement;
        private readonly IEmbeddedDOFInHostTransformationVector transformation;
        private readonly Dictionary<SuperElementDOF, int> superElementMap = new Dictionary<SuperElementDOF, int>();
        private readonly Dictionary<EmbeddedNode, Dictionary<DOFType, int>> dofToHostMapping = new Dictionary<EmbeddedNode, Dictionary<DOFType, int>>();
        private Matrix2D transformationMatrix;
        //private bool isElementEmbedded = false;

        public ElementEmbedder(Model model, Element embeddedElement, IEmbeddedDOFInHostTransformationVector transformation)
        {
            this.model = model;
            this.embeddedElement = embeddedElement;
            this.transformation = transformation;
            Initialize();
        }

        private void InitializeMappings()
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            superElementMap.Clear();
            int index = 0;
            foreach (var embeddedNode in e.EmbeddedNodes)
            {
                int nodeOrderInEmbeddedElement = embeddedElement.Nodes.IndexOf(embeddedNode.Node);
                var currentEmbeddedNodeDOFs = embeddedElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement];
                //var currentNodeDOFs = currentEmbeddedNodeDOFs.Intersect(embeddedNode.DependentDOFs);
                var independentEmbeddedDOFs = currentEmbeddedNodeDOFs.Except(embeddedNode.DependentDOFs);

                // TODO: Optimization to exclude host DOFs that embedded node does not depend on.
                for (int i = 0; i < embeddedNode.EmbeddedInElement.Nodes.Count; i++)
                {
                    var currentNodeDOFs = embeddedNode.EmbeddedInElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedNode.EmbeddedInElement)[i];
                    foreach (var dof in currentNodeDOFs)
                    {
                        var superElementDOF = new SuperElementDOF() { DOF = dof, EmbeddedNode = embeddedNode.Node, HostNode = embeddedNode.EmbeddedInElement.Nodes[i], Element = embeddedNode.EmbeddedInElement };
                        if (!superElementMap.ContainsKey(superElementDOF))
                        {
                            superElementMap.Add(superElementDOF, index);
                            index++;
                        }
                    }
                }

                //var independentEmbeddedDOFs = model.NodalDOFsDictionary[embeddedNode.Node.ID].Select(x => x.Key).Except(embeddedNode.DependentDOFs);
                //int nodeOrderInEmbeddedElement = embeddedElement.Nodes.IndexOf(embeddedNode.Node);

                //var independentEmbeddedDOFs = embeddedElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement].Except(embeddedNode.DependentDOFs);

                foreach (var dof in independentEmbeddedDOFs)
                {
                    var superElementDOF = new SuperElementDOF() { DOF = dof, EmbeddedNode = embeddedNode.Node, HostNode = null, Element = null };
                    if (!superElementMap.ContainsKey(superElementDOF))
                    {
                        superElementMap.Add(superElementDOF, index);
                        index++;
                    }
                }
            }

            foreach (var node in embeddedElement.Nodes.Except(e.EmbeddedNodes.Select(x => x.Node)))
            {
                int nodeOrderInEmbeddedElement = embeddedElement.Nodes.IndexOf(node);
                var currentNodeDOFs = embeddedElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement];
                foreach (var dof in currentNodeDOFs)
                {
                    var superElementDOF = new SuperElementDOF() { DOF = dof, EmbeddedNode = node, HostNode = null, Element = null };
                    if (!superElementMap.ContainsKey(superElementDOF))
                    {
                        superElementMap.Add(superElementDOF, index);
                        index++;
                    }
                }
            }
        }

        private void CalculateTransformationMatrix()
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            int row = 0;
            int col = 0;
            int totalRows = embeddedElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedElement).SelectMany(x => x).Count();
            var matrix = new double[totalRows, superElementMap.Count];

            foreach (var embeddedNode in e.EmbeddedNodes)
            {
                var localTransformationMatrix = transformation.GetTransformationVector(embeddedNode);
                var localHostDOFs = transformation.GetDOFTypesOfHost(embeddedNode);
                int nodeOrderInEmbeddedElement = embeddedElement.Nodes.IndexOf(embeddedNode.Node);
                var embeddedNodeDOFQuantity = embeddedElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement].Count;
                int dependentDOFs = transformation.GetDependentDOFTypes.Count;

                for (int i = 0; i < dependentDOFs; i++)
                {
                    col = 0;
                    for (int j = 0; j < localHostDOFs.Count; j++)
                    {
                        for (int k = 0; k < localHostDOFs[j].Count; k++)
                        {
                            var superelement = new SuperElementDOF() { DOF = localHostDOFs[j][k], Element = embeddedNode.EmbeddedInElement, EmbeddedNode = embeddedNode.Node, HostNode = embeddedNode.EmbeddedInElement.Nodes[j] };
                            matrix[row + i, superElementMap[superelement]] = localTransformationMatrix[i][col];
                            col++;
                        }
                    }
                }
                row += dependentDOFs;

                var independentEmbeddedDOFs = embeddedElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement].Except(embeddedNode.DependentDOFs).ToArray();
                for (int j = 0; j < independentEmbeddedDOFs.Length; j++)
                {
                    var superelement = new SuperElementDOF() { DOF = independentEmbeddedDOFs[j], Element = null, HostNode = null, EmbeddedNode = embeddedNode.Node };
                    matrix[row, superElementMap[superelement]] = 1;
                    row++;
                }
            }

            foreach (var node in embeddedElement.Nodes.Except(e.EmbeddedNodes.Select(x => x.Node)))
            {
                int nodeOrderInEmbeddedElement = embeddedElement.Nodes.IndexOf(node);
                var currentNodeDOFs = embeddedElement.ElementType.DOFEnumerator.GetDOFTypes(embeddedElement)[nodeOrderInEmbeddedElement];
                for (int j = 0; j < currentNodeDOFs.Count; j++)
                {
                    var superelement = new SuperElementDOF() { DOF = currentNodeDOFs[j], Element = null, HostNode = null, EmbeddedNode = node };
                    matrix[row, superElementMap[superelement]] = 1;
                    row++;
                }
            }

            //StreamWriter sw = File.CreateText(String.Format("TransformationMatrix{0}.txt", embeddedElement.ID));
            //for (int i = 0; i < totalRows; i++)
            //{
            //    var line = String.Empty;
            //    for (int j = 0; j < superElementMap.Count; j++)
            //        line += matrix[i,j].ToString() + ";";
            //    sw.WriteLine(line);
            //}
            //sw.Close();
            transformationMatrix = new Matrix2D(matrix);
        }

        private void Initialize()
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            if (e == null) return;
            if (e.EmbeddedNodes.Count == 0) return;

            InitializeMappings();
            CalculateTransformationMatrix();
        }

        public IMatrix2D GetTransformedMatrix(IMatrix2D matrix)
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            //if (e == null || !isElementEmbedded) return matrix;
            if (e == null) return matrix;
            if (e.EmbeddedNodes.Count == 0) return matrix;

            return transformationMatrix.Transpose() * ((SymmetricMatrix2D)matrix).ToMatrix2D() * transformationMatrix;
        }

        public double[] GetTransformedVector(double[] vector)
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            //if (e == null || !isElementEmbedded) return matrix;
            if (e == null) return vector;
            if (e.EmbeddedNodes.Count == 0) return vector;

            return (transformationMatrix * new Vector(vector)).Data;
        }

        public double[] GetTransformedForcesVector(double[] vector) //compa prosthiki msolve
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            //if (e == null || !isElementEmbedded) return matrix;
            if (e == null) return vector;
            if (e.EmbeddedNodes.Count == 0) return vector;

            return (transformationMatrix.Transpose() * new Vector(vector)).Data; // compa Vector<double>
        }

        public IList<IList<DOFType>> GetDOFTypes(Element element)
        {
            //return element.ElementType.GetElementDOFTypes(element);

            var dofs = new List<IList<DOFType>>();
            Node currentNode = null;
            List<DOFType> nodeDOFs = null;

            foreach (var superElement in superElementMap)
            {
                Node node = superElement.Key.HostNode == null ? superElement.Key.EmbeddedNode : superElement.Key.HostNode;

                if (currentNode != node)
                {
                    if (nodeDOFs != null)
                        dofs.Add(nodeDOFs);
                    currentNode = node;
                    nodeDOFs = new List<DOFType>();
                }
                nodeDOFs.Add(superElement.Key.DOF);
            }
            if (nodeDOFs != null)
                dofs.Add(nodeDOFs);

            return dofs;
        }

        public IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(Element element)
        {
            //if (embeddedElement != element) throw new ArgumentException();

            var nodesDictionary = new Dictionary<Node, int>();
            int index = 0;
            foreach (var node in element.Nodes)
            {
                nodesDictionary.Add(node, index);
                index++;
            }
            
            var dofs = new List<IList<DOFType>>();
            for (int i = 0; i < element.Nodes.Count; i++)
                dofs.Add(new List<DOFType>());

            Node currentNode = null;
            List<DOFType> nodeDOFs = null;

            foreach (var superElement in superElementMap)
            {
                if (superElement.Key.HostNode != null) continue;
                Node node = superElement.Key.EmbeddedNode;
                //Node node = superElement.Key.HostNode == null ? superElement.Key.EmbeddedNode : superElement.Key.HostNode;

                if (currentNode != node)
                {
                    if (nodeDOFs != null)
                        dofs[nodesDictionary[currentNode]] = nodeDOFs;
                    currentNode = node;
                    nodeDOFs = new List<DOFType>();
                }
                nodeDOFs.Add(superElement.Key.DOF);
            }
            if (nodeDOFs != null)
                dofs[nodesDictionary[currentNode]] = nodeDOFs;
            //dofs.Add(nodeDOFs);

            return dofs;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            var nodes = new List<Node>();
            Node currentNode = null;
            foreach (var superElement in superElementMap)
            {
                Node node = superElement.Key.HostNode == null ? superElement.Key.EmbeddedNode : superElement.Key.HostNode;
                if (currentNode != node)
                {
                    if (currentNode != null) 
                        nodes.Add(currentNode);
                    currentNode = node;
                }
                //if (nodes.IndexOf(node) < 0)
                //    nodes.Add(node);
            }
            if (currentNode != null)
                nodes.Add(currentNode);

            return nodes;
        }
    }
}
