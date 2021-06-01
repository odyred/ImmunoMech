﻿using ISAAR.MSolve.FEM.Entities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;

namespace ISAAR.MSolve.FEM
{
    public class AutomaticDomainDecomposer2 // exei allaxei mono to onoma den exei ginei update se model xrhsh klp.
    {
        private readonly Model Model;
        private readonly int NumberOfProcessors;
        private int numberOfElementsPerSubdomain;
        Dictionary<Element, List<Element>> ElementAdjacency;
        List<Element> ElementsRenumbered = new List<Element>();
        Dictionary<int,List<Node>> SubdomainInterfaceNodes =new Dictionary<int, List<Node>>();

        public AutomaticDomainDecomposer2(Model model, int numberOfProcessors)
        {
            this.NumberOfProcessors = numberOfProcessors;
            this.Model = model;
        }


        public void UpdateModel(bool isColoringEnabled = false)
        {
            Adjacency();

            CreateSubdomains();
            AssignElementsRenumberedToSubdomains();

            if (isColoringEnabled)
            {
                var purgedElements = Purge();
                ColorDisconnectedElements(purgedElements);
            }
        }

        private void AssignElementsRenumberedToSubdomains()
        {
            Model.SubdomainsDictionary.Clear();
            var indexElement = 0;
            for (int i = 0; i < NumberOfProcessors; i++)
            {
                if (indexElement >= ElementsRenumbered.Count) break;
                Model.SubdomainsDictionary.Add(i, new Subdomain(i)); //new Subdomain() { ID = i }
                for (int j = 0; j < numberOfElementsPerSubdomain; j++)
                {
                    if (indexElement >= ElementsRenumbered.Count) break;
                    Model.SubdomainsDictionary[i].Elements.Add( ElementsRenumbered[indexElement++]);
                }
            }

            UpdateModelDataStructures();
        }

        private void UpdateModelDataStructures()
        {
            foreach (Subdomain subdomain in Model.SubdomainsDictionary.Values)
            {
                foreach (Element element in subdomain.Elements)
                    element.Subdomain = subdomain;
            }

            foreach (Node node in Model.NodesDictionary.Values)
            {
                node.SubdomainsDictionary.Clear();
                node.BuildSubdomainDictionary();
            }
            //TEMP comment
            //foreach (Subdomain subdomain in Model.SubdomainsDictionary.Values)
            //{
            //    subdomain.NodesDictionary.Clear();
            //    subdomain.BuildNodesDictionary();
            //}
        }

        private void Adjacency()
        {
            ElementAdjacency = new Dictionary<Element, List<Element>>();
            // mask is an integer that shows if the element is used
            

            foreach (var element in Model.ElementsDictionary.Values)
            {
                var usedElement = new Dictionary<Element, bool>(Model.ElementsDictionary.Count);//bool[] usedElement = new bool[Model.ElementsDictionary.Count];//mask
                foreach(Element e1 in Model.ElementsDictionary.Values) { usedElement.Add(e1, false); }

                ElementAdjacency.Add(element,new List<Element>());
                usedElement[element] = true;
                foreach (var node in element.NodesDictionary.Values)
                {
                    foreach (var nodeElement in node.ElementsDictionary.Values)
                    {
                        if (!usedElement.ContainsKey(nodeElement)) { usedElement.Add(nodeElement, false); };

                        if (usedElement[nodeElement]) continue;
                        ElementAdjacency[element].Add(nodeElement);
                        usedElement[nodeElement] = true;
                    }
                }
            }
        }
        
        private void CreateSubdomains()
        {
            //TODO: return IndexOutOfRangeException if nodes,elements or subdomains numbering does not start with 0.
            var isInteriorBoundaryElement = new Dictionary<Element, bool>(Model.ElementsDictionary.Count); //bool[] isInteriorBoundaryElement= new bool[Model.ElementsDictionary.Count];
            foreach(Element element in Model.ElementsDictionary.Values) { isInteriorBoundaryElement.Add(element, false); }
            var isInteriorBoundaryNode = new Dictionary<Node, bool>(Model.NodesDictionary.Count); //bool[] isInteriorBoundaryNode = new bool[Model.NodesDictionary.Count];
            foreach(Node node in Model.NodesDictionary.Values) { isInteriorBoundaryNode.Add(node, false); }

            // Number of Elements per subdomain
            numberOfElementsPerSubdomain =(Model.ElementsDictionary.Count % NumberOfProcessors==0)?
                Model.ElementsDictionary.Count / NumberOfProcessors: Model.ElementsDictionary.Count / NumberOfProcessors+1;

            Dictionary<Node, int> nodeWeight = new Dictionary<Node, int>();
            foreach (Node node in Model.NodesDictionary.Values)
                nodeWeight.Add(node, node.ElementsDictionary.Count);

            var usedElementsCounter = 0;
            var mlabel = 0;
            int counterSubdomain = 0;
            
            do
            {
                var flag = true;
                var flagStop=true;
                #region Find Node with next minimum weight
                var finalSubdomainElement = usedElementsCounter;
                var minimumNodeWeight = int.MaxValue;
                int nodeID = Model.NodesDictionary.ElementAt(0).Value.ID;// int nodeID = 0; ;
                for (int i = 0; i < Model.NodesDictionary.Count; i++)
                {
                    if (nodeWeight[Model.NodesDictionary.ElementAt(i).Value] == 0) continue;
                    if (nodeWeight[Model.NodesDictionary.ElementAt(i).Value] < minimumNodeWeight)
                    {
                        minimumNodeWeight = nodeWeight[Model.NodesDictionary.ElementAt(i).Value];
                        nodeID = Model.NodesDictionary.ElementAt(i).Value.ID;
                    }
                }
                #endregion

                // Start fill list with elements connected to node with minimum weight
                var counterSubdomainElements = 0;
                foreach (var element in Model.NodesDictionary[nodeID].ElementsDictionary.Values)
                {
                    var elementID = element.ID;
                    if (isInteriorBoundaryElement[element]) continue;
                    counterSubdomainElements++;
                    isInteriorBoundaryElement[element] = true;
                    ElementsRenumbered.Add(element);

                    #region nomask
                    //Reduce nodeWeight for all nodes connected to this element
                    foreach (Node node in Model.ElementsDictionary[elementID].NodesDictionary.Values)
                        nodeWeight[node]--;
                    #endregion

                    if (counterSubdomainElements == numberOfElementsPerSubdomain)
                    {
                        flag = false;
                        break;
                    }
                }

                if (flag)
                {
                    // Recursively add adjacent elements to list
                    do
                    {
                        var initialSubdomainElement = finalSubdomainElement;
                        finalSubdomainElement = usedElementsCounter + counterSubdomainElements;
                        var nnstart = initialSubdomainElement + 1;
                        
                        for (int i = initialSubdomainElement; i <= finalSubdomainElement-1; i++)
                        {
                            int lc = 0;
                            for (int j=0; j< ElementAdjacency.First(x=>x.Key.ID==ElementsRenumbered[i].ID).Value.Count; j++)
                            {
                                var element = ElementAdjacency[ElementsRenumbered[i]][j]; int elementID = element.ID;  //int elementID = ElementAdjacency[ElementsRenumbered[i]][j].ID;
                                if (isInteriorBoundaryElement[element]) continue;
                                lc++;
                                counterSubdomainElements++;
                                isInteriorBoundaryElement[element] = true;
                                ElementsRenumbered.Add(Model.ElementsDictionary.First(x =>x.Value.ID== elementID).Value);

                                #region nomask
                                foreach (Node node in Model.ElementsDictionary[elementID].NodesDictionary.Values)
                                    nodeWeight[node]--;
                                #endregion

                                if (counterSubdomainElements== numberOfElementsPerSubdomain)
                                {
                                    flag = false;
                                    break;
                                }
                            } // 800

                            if (flag)
                            {
                                if (lc == 0 && (usedElementsCounter + counterSubdomainElements) == finalSubdomainElement && i == finalSubdomainElement-1)
                                {
                                    usedElementsCounter = usedElementsCounter + counterSubdomainElements;
                                    flagStop = false;
                                    flag = false;
                                }
                            }

                            if (!flag) break;
                        } 
                        if (!flag) break;
                    } while (counterSubdomainElements < numberOfElementsPerSubdomain);
                    if (!flagStop) continue;
                }
                SubdomainInterfaceNodes.Add(counterSubdomain++,CalculateInterface(nodeWeight, isInteriorBoundaryNode, ElementsRenumbered, usedElementsCounter, counterSubdomainElements, counterSubdomain));
                usedElementsCounter = usedElementsCounter + counterSubdomainElements;
                mlabel = usedElementsCounter;
            } while (usedElementsCounter < Model.ElementsDictionary.Count);

        }

        //isInteriorBoundaryNode-> if node is on the interior interface
        private List<Node>  CalculateInterface(Dictionary<Node,int> nodeWeight, Dictionary<Node,bool> isInteriorBoundaryNode, List<Element>  ElementsRenumbered, int usedElementsCounter, int counterSubdomainElements, int counterSubdomain)
        {
            var locmask = new Dictionary<Node, bool>(Model.NodesDictionary.Count);// new bool[Model.NodesDictionary.Count];
            foreach(Node node in Model.NodesDictionary.Values) { locmask.Add(node, false); }

            List<Node> SubdomainInterfaceNodes = new List<Node>();            

            for (int i = usedElementsCounter; i < usedElementsCounter + counterSubdomainElements; i++)
            {
                int elementID = ElementsRenumbered[i].ID;
                foreach (Node node in Model.ElementsDictionary[elementID].NodesDictionary.Values)
                {
                    if ((nodeWeight[node] !=0||isInteriorBoundaryNode[node])&&!locmask[node])
                    {
                        isInteriorBoundaryNode[node] = true;
                        locmask[node] = true;
                        SubdomainInterfaceNodes.Add(node);
                    }
                }
            }
            return SubdomainInterfaceNodes;
        }

        private void ColorDisconnectedElements(List<Element> purgedElements)
        {
            int numberOfColors = 1;
            Dictionary<int, List<Element>> elementsPerColor = new Dictionary<int, List<Element>>();
            int indexColor = 0;
            int counterElement = 0;

            var isElementUsed = new Dictionary<Element, bool>(Model.ElementsDictionary.Count);//bool[] isElementUsed = new bool[Model.ElementsDictionary.Count];

            // Form the list of distinct colors
            do
            {
                elementsPerColor.Add(indexColor,new List<Element>());
                var usedNodes = new Dictionary<Node, bool>(Model.NodesDictionary.Count);// bool[] usedNodes = new bool[Model.NodesDictionary.Count];
                foreach (var element in purgedElements)
                {
                    if (isElementUsed[element]) continue;

                    bool disjoint = CheckIfElementDisjoint(usedNodes, element.ID);

                    if (disjoint)
                    {
                        counterElement++;
                        elementsPerColor[indexColor].Add(Model.ElementsDictionary[element.ID]);

                        #region marker
                        foreach (Node node in Model.ElementsDictionary[element.ID].NodesDictionary.Values)
                            usedNodes[node] = true;
                        #endregion

                    }

                }
                indexColor++;
            } while (counterElement < purgedElements.Count);
            numberOfColors = indexColor;
        }

        private bool CheckIfElementDisjoint(Dictionary<Node, bool> usedNodes, int elementID)
        {
            foreach (Node node in Model.ElementsDictionary[elementID].NodesDictionary.Values)
                if (usedNodes[node])
                    return false;
            return true;
        }


        private List<Element> Purge()
        {
            var isElementUsed = new Dictionary<Element, bool>(Model.ElementsDictionary.Count);// bool[] isElementUsed = new bool[Model.ElementsDictionary.Count];

            List<Node> interfaceNodes = new List<Node>();
            foreach (var nodeList in SubdomainInterfaceNodes.Values)
                interfaceNodes.AddRange(nodeList);

            interfaceNodes=interfaceNodes.Distinct().ToList();
            int numberOfInterfaceNodes = interfaceNodes.Count;

            List<Element> purgedElements=new List<Element>();
            int numberOfPurgedElements = 0;
            for (int indexInterfaceNode = 0; indexInterfaceNode < numberOfInterfaceNodes; indexInterfaceNode++)
            {
                foreach (Element element in interfaceNodes[indexInterfaceNode].ElementsDictionary.Values)
                {
                    if (!isElementUsed[element])
                    {
                        isElementUsed[element] = true;
                        purgedElements.Add(element);
                    }
                }
            }

            return purgedElements;
        }

    }
}
