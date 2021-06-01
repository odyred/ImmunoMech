﻿using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra;
using ISAAR.MSolve.LinearAlgebra.Matrices;

namespace ISAAR.MSolve.FEM.Embedding
{
    public class ElementEmbedder : IElementDofEnumerator
    {
        private readonly Model model;
        private readonly Element embeddedElement;
        private readonly IEmbeddedDOFInHostTransformationVector transformation;
        private readonly Dictionary<SuperElementDof, int> superElementMap = new Dictionary<SuperElementDof, int>();
        private readonly Dictionary<EmbeddedNode, Dictionary<IDofType, int>> dofToHostMapping = new Dictionary<EmbeddedNode, Dictionary<IDofType, int>>();
        private Matrix transformationMatrix;
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
            var e = (IEmbeddedElement)(embeddedElement.ElementType);
            superElementMap.Clear();
            int index = 0;
            foreach (var embeddedNode in e.EmbeddedNodes)
            {
                int nodeOrderInEmbeddedElement = embeddedElement.Nodes.IndexOf(embeddedNode.Node);
                var currentEmbeddedNodeDOFs = embeddedElement.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(embeddedElement)[nodeOrderInEmbeddedElement];
                //var currentNodeDOFs = currentEmbeddedNodeDOFs.Intersect(embeddedNode.DependentDOFs);
                var independentEmbeddedDOFs = currentEmbeddedNodeDOFs.Except(embeddedNode.DependentDOFs);

                // TODO: Optimization to exclude host DOFs that embedded node does not depend on.
                for (int i = 0; i < embeddedNode.EmbeddedInElement.Nodes.Count; i++)
                {
                    var currentNodeDOFs = embeddedNode.EmbeddedInElement.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(embeddedNode.EmbeddedInElement)[i];
                    foreach (var dof in currentNodeDOFs)
                    {
                        var superElementDOF = new SuperElementDof() { DOF = dof, EmbeddedNode = embeddedNode.Node, HostNode = embeddedNode.EmbeddedInElement.Nodes[i], Element = embeddedNode.EmbeddedInElement };
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
                    var superElementDOF = new SuperElementDof() { DOF = dof, EmbeddedNode = embeddedNode.Node, HostNode = null, Element = null };
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
                var currentNodeDOFs = embeddedElement.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(embeddedElement)[nodeOrderInEmbeddedElement];
                foreach (var dof in currentNodeDOFs)
                {
                    var superElementDOF = new SuperElementDof() { DOF = dof, EmbeddedNode = node, HostNode = null, Element = null };
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
            var e = (IEmbeddedElement)(embeddedElement.ElementType);
            int row = 0;
            int col = 0;
            int totalRows = embeddedElement.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(embeddedElement).SelectMany(x => x).Count();
            var matrix = new double[totalRows, superElementMap.Count];

            foreach (var embeddedNode in e.EmbeddedNodes)
            {
                var localTransformationMatrix = transformation.GetTransformationVector(embeddedNode);
                var localHostDOFs = transformation.GetDOFTypesOfHost(embeddedNode);
                int nodeOrderInEmbeddedElement = embeddedElement.Nodes.IndexOf(embeddedNode.Node);
                var embeddedNodeDOFQuantity = embeddedElement.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(embeddedElement)[nodeOrderInEmbeddedElement].Count;
                int dependentDOFs = transformation.GetDependentDOFTypes.Count;

                for (int i = 0; i < dependentDOFs; i++)
                {
                    col = 0;
                    for (int j = 0; j < localHostDOFs.Count; j++)
                    {
                        for (int k = 0; k < localHostDOFs[j].Count; k++)
                        {
                            var superelement = new SuperElementDof() { DOF = localHostDOFs[j][k], Element = embeddedNode.EmbeddedInElement, EmbeddedNode = embeddedNode.Node, HostNode = embeddedNode.EmbeddedInElement.Nodes[j] };
                            matrix[row + i, superElementMap[superelement]] = localTransformationMatrix[i][col];
                            col++;
                        }
                    }
                }
                row += dependentDOFs;

                var independentEmbeddedDOFs = embeddedElement.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(embeddedElement)[nodeOrderInEmbeddedElement].Except(embeddedNode.DependentDOFs).ToArray();
                for (int j = 0; j < independentEmbeddedDOFs.Length; j++)
                {
                    var superelement = new SuperElementDof() { DOF = independentEmbeddedDOFs[j], Element = null, HostNode = null, EmbeddedNode = embeddedNode.Node };
                    matrix[row, superElementMap[superelement]] = 1;
                    row++;
                }
            }

            foreach (var node in embeddedElement.Nodes.Except(e.EmbeddedNodes.Select(x => x.Node)))
            {
                int nodeOrderInEmbeddedElement = embeddedElement.Nodes.IndexOf(node);
                var currentNodeDOFs = embeddedElement.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(embeddedElement)[nodeOrderInEmbeddedElement];
                for (int j = 0; j < currentNodeDOFs.Count; j++)
                {
                    var superelement = new SuperElementDof() { DOF = currentNodeDOFs[j], Element = null, HostNode = null, EmbeddedNode = node };
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
            transformationMatrix = Matrix.CreateFromArray(matrix);
        }

        private void Initialize()
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            if (e == null) return;
            if (e.EmbeddedNodes.Count == 0) return;

            InitializeMappings();
            CalculateTransformationMatrix();
        }

        public IMatrix GetTransformedMatrix(IMatrix matrix)
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            //if (e == null || !isElementEmbedded) return matrix;
            if (e == null) return matrix;
            if (e.EmbeddedNodes.Count == 0) return matrix;

            return transformationMatrix.ThisTransposeTimesOtherTimesThis(matrix);
        }

        public double[] GetTransformedDisplacementsVector(double[] vector)
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            //if (e == null || !isElementEmbedded) return matrix;
            if (e == null) return vector;
            if (e.EmbeddedNodes.Count == 0) return vector;

            return transformationMatrix.Multiply(vector, false);
        }

        public double[] GetTransformedForcesVector(double[] vector) //compa prosthiki msolve
        {
            var e = embeddedElement.ElementType as IEmbeddedElement;
            //if (e == null || !isElementEmbedded) return matrix;
            if (e == null) return vector;
            if (e.EmbeddedNodes.Count == 0) return vector;

            return transformationMatrix.Multiply(vector, true);
        }


        public IReadOnlyList<IReadOnlyList<IDofType>> GetDofTypesForMatrixAssembly(IElement element)
        {
            //return element.ElementType.GetElementDOFTypes(element);

            var dofs = new List<IReadOnlyList<IDofType>>();
            INode currentNode = null;
            List<IDofType> nodeDOFs = null;

            foreach (var superElement in superElementMap)
            {
                INode node = superElement.Key.HostNode == null ? superElement.Key.EmbeddedNode : superElement.Key.HostNode;

                if (currentNode != node)
                {
                    if (nodeDOFs != null)
                        dofs.Add(nodeDOFs);
                    currentNode = node;
                    nodeDOFs = new List<IDofType>();
                }
                nodeDOFs.Add(superElement.Key.DOF);
            }
            if (nodeDOFs != null)
                dofs.Add(nodeDOFs);

            return dofs;
        }

        public IReadOnlyList<IReadOnlyList<IDofType>> GetDofTypesForDofEnumeration(IElement element)
        {
            //if (embeddedElement != element) throw new ArgumentException();

            var nodesDictionary = new Dictionary<INode, int>();
            int index = 0;
            foreach (var node in element.Nodes)
            {
                nodesDictionary.Add(node, index);
                index++;
            }
            
            var dofs = new List<IReadOnlyList<IDofType>>();
            for (int i = 0; i < element.Nodes.Count; i++)
                dofs.Add(new List<IDofType>());

            INode currentNode = null;
            List<IDofType> nodeDOFs = null;

            foreach (var superElement in superElementMap)
            {
                if (superElement.Key.HostNode != null) continue;
                INode node = superElement.Key.EmbeddedNode;
                //Node node = superElement.Key.HostNode == null ? superElement.Key.EmbeddedNode : superElement.Key.HostNode;

                if (currentNode != node)
                {
                    if (nodeDOFs != null)
                        dofs[nodesDictionary[currentNode]] = nodeDOFs;
                    currentNode = node;
                    nodeDOFs = new List<IDofType>();
                }
                nodeDOFs.Add(superElement.Key.DOF);
            }
            if (nodeDOFs != null)
                dofs[nodesDictionary[currentNode]] = nodeDOFs;
            //dofs.Add(nodeDOFs);

            return dofs;
        }

        public IReadOnlyList<INode> GetNodesForMatrixAssembly(IElement element)
        {
            var nodes = new List<INode>();
            INode currentNode = null;
            foreach (var superElement in superElementMap)
            {
                INode node = superElement.Key.HostNode == null ? superElement.Key.EmbeddedNode : superElement.Key.HostNode;
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
