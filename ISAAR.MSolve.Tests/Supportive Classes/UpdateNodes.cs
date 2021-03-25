using System;
using System.Linq;
using System.Text;
using System.Collections.Generic;
using System.Globalization;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.Discretization.Mesh;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.Materials;
using ISAAR.MSolve.FEM.Elements.BoundaryConditionElements;
using ISAAR.MSolve.FEM.Readers.Interfaces;

namespace ISAAR.MSolve.Tests.Supportive_Classes
{
    public class UpdateNodes : IModelReader
    {
        public Model Model { get; private set; }
        public IList<IList<Node>> nodeBoundaries;
        public IList<IList<Element>> elementBoundaries;
        public IList<IList<Element>> elementDomains { get; private set; }
        public IList<IList<Node>> nodeDomains { get; private set; }
        public int[] DomainIDs { get; set; }
        public IList<IList<IList<Node>>> quadBoundaries { get; private set; }
        public IList<IList<IList<Node>>> triBoundaries { get; private set; }
        private double[] diffusionCoeff;
        private double[][] convectionCoeff;
        private double[] loadFromUnknownCoeff;
        enum Attributes
        {
            sdim = 1001,
            Meshpointcoordinates = 1002,
            numberofmeshpoints = 1003,
            numberofelementtypes = 1004,
            vtx = 1006,
            edg = 1007,
            tri = 1008,
            quad = 1009,
            tet = 1010,
            hex = 1011,
        }

        public string Filename { get; private set; }

        public UpdateNodes(string filename, double[] k, double[][] U, double[] L)
        {
            Filename = filename;
            diffusionCoeff = k;
            convectionCoeff = U;
            loadFromUnknownCoeff = L;
        }

        int NumberOfNodes;
        int NumberOfTriElements;
        int NumberOfTetElements;
        int NumberOfHexElements;
        int NumberOfQuadElements;

        public UpdateNodes UpdateModelReader(double[] k, double[][] U, double[] L)
        {
            diffusionCoeff = k;
            convectionCoeff = U;
            loadFromUnknownCoeff = L;
            return this;
        }

        public Model UpdateModel()
        {//only works for tet4 now
            ConvectionDiffusionMaterial[] CDMaterial = new ConvectionDiffusionMaterial[diffusionCoeff.Length];
            ConvectionDiffusionElement3DFactory[] elementFactory3D = new ConvectionDiffusionElement3DFactory[diffusionCoeff.Length];
            for (int i = 0; i < diffusionCoeff.Length; i++)
            {
                CDMaterial[i] = new ConvectionDiffusionMaterial(1, diffusionCoeff[i], convectionCoeff[i], loadFromUnknownCoeff[i]);
                elementFactory3D[i] = new ConvectionDiffusionElement3DFactory(CDMaterial[i]);
            }
            for (int domain = 0; domain < elementDomains.Count; domain++)
                foreach (Element elem in elementDomains[domain])
                    elem.ElementType = elementFactory3D[domain].CreateElement(CellType.Tet4, elem.Nodes.ToList());
            Model.BodyLoads.Clear();
            return Model;
        }

        public Model CreateModelFromFile(int[] domainIDs = null, int[] boundaryIDs = null)
        {
            ConvectionDiffusionMaterial[] CDMaterial = new ConvectionDiffusionMaterial[diffusionCoeff.Length];
            ConvectionDiffusionElement3DFactory[] elementFactory3D = new ConvectionDiffusionElement3DFactory[diffusionCoeff.Length];
            for (int i = 0; i < diffusionCoeff.Length; i++)
            {
                CDMaterial[i] = new ConvectionDiffusionMaterial(1, diffusionCoeff[i], convectionCoeff[i], loadFromUnknownCoeff[i]);
                elementFactory3D[i] = new ConvectionDiffusionElement3DFactory(CDMaterial[i]);
            }
            //var boundaryFactory3D = new SurfaceBoundaryFactory3D(0, new ConvectionDiffusionMaterial(diffusionCoeff, new double[] {0,0,0}, 0));
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(0);
            // Material

            char[] delimeters = { ' ', '=', '\t' };
            Attributes? name = null;

            String[] text = System.IO.File.ReadAllLines(Filename);
            elementBoundaries = new List<IList<Element>>();
            elementDomains = new List<IList<Element>>();
            nodeDomains = new List<IList<Node>>();
            nodeBoundaries = new List<IList<Node>>();
            quadBoundaries = new List<IList<IList<Node>>>();
            triBoundaries = new List<IList<IList<Node>>>();
            for (int i = 0; i < 10; i++)
            {
                elementBoundaries.Add(new List<Element>());
                nodeBoundaries.Add(new List<Node>());
                quadBoundaries.Add(new List<IList<Node>>());
                triBoundaries.Add(new List<IList<Node>>());
            }
            for (int i = 0; i < 2; i++)
            {
                elementDomains.Add(new List<Element>());
                nodeDomains.Add(new List<Node>());
            }

            for (int i = 0; i < text.Length; i++)
            {
                String[] line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
                StringBuilder comparisonString = new StringBuilder(null, 50);
                if (line.Length == 0)
                {
                    continue;
                }
                if ((line[0] == "3" || line[0] == "4") && line.Length > 1)
                {
                    try
                    {
                        name = (Attributes)Enum.Parse(typeof(Attributes), line[1]);
                    }
                    catch (Exception exception)
                    {
                        name = null;
                    }
                }
                else
                {
                    for (int linePosition = 0; linePosition < line.Length; linePosition++)
                    {
                        if (line[linePosition] == "#")
                        {
                            for (int ij = linePosition + 1; ij < line.Length; ij++)
                            {
                                comparisonString.Append(line[ij]);
                            }
                            try
                            {
                                name = (Attributes)Enum.Parse(typeof(Attributes), comparisonString.ToString());
                            }
                            catch (Exception exception)
                            {
                                name = null;
                            }
                        }
                    }
                }
                switch (name)
                {
                    case Attributes.sdim:
                        int NumberOfDimensions = Int32.Parse(line[0]);
                        break;
                    case Attributes.numberofmeshpoints:
                        NumberOfNodes = Int32.Parse(line[0]);
                        break;
                    case Attributes.Meshpointcoordinates:
                        IList<Node> nodelist = new List<Node> { null };
                        for (int j = 0; j < NumberOfNodes; j++)
                        {
                            i++;
                            line = text[i].Split(delimeters);
                            int nodeGlobalID = j;
                            double x = Double.Parse(line[0], CultureInfo.InvariantCulture);
                            double y = Double.Parse(line[1], CultureInfo.InvariantCulture);
                            double z = Double.Parse(line[2], CultureInfo.InvariantCulture);
                            Node node = new Node(nodeGlobalID, x, y, z);
                            model.NodesDictionary.Add(nodeGlobalID, node);
                            nodelist.Add(node);
                        }

                        break;
                    case Attributes.tri:
                        do
                        {
                            i++;
                            line = text[i].Split(delimeters);
                        } while (line[0] == "");
                        i++;
                        line = text[i].Split(delimeters);
                        NumberOfTriElements = Int32.Parse(line[0]);
                        i++;
                        IList<IList<Node>> triNodesCollection = new List<IList<Node>>();
                        for (int TriID = 0; TriID < NumberOfTriElements; TriID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);

                            triNodesCollection.Add(new List<Node>
                            {
                                model.NodesDictionary[Int32.Parse(line[2])],
                                model.NodesDictionary[Int32.Parse(line[1])],
                                model.NodesDictionary[Int32.Parse(line[0])]
                            });

                        }
                        i = i + 3;
                        for (int TriID = 0; TriID < NumberOfTriElements; TriID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);
                            int boundaryID = Int32.Parse(line[0]);
                            bool res = Array.Exists(boundaryIDs, x => x == boundaryID - 1);
                            if (res)
                            {
                                triBoundaries[boundaryID].Add(triNodesCollection[TriID]);
                                foreach (Node node in triNodesCollection[TriID])
                                {
                                    nodeBoundaries[boundaryID].Add(node);
                                }
                            }
                            //int elementBoundaryID = Int32.Parse(line[0]);

                            //elementBoundaries[elementBoundaryID].Add(model.ElementsDictionary[QuadID]);
                        }

                        foreach (int boundaryID in boundaryIDs)
                        {
                            nodeBoundaries[boundaryID] = nodeBoundaries[boundaryID].Distinct().ToList();
                        }
                        break;
                    case Attributes.quad:
                        do
                        {
                            i++;
                            line = text[i].Split(delimeters);
                        } while (line[0] == "");
                        i++;
                        line = text[i].Split(delimeters);
                        NumberOfQuadElements = Int32.Parse(line[0]);
                        i++;
                        IList<IList<Node>> quadNodesCollection = new List<IList<Node>>();
                        for (int QuadID = 0; QuadID < NumberOfQuadElements; QuadID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);

                            quadNodesCollection.Add(new List<Node>
                            {
                                model.NodesDictionary[Int32.Parse(line[1])], //0
                                model.NodesDictionary[Int32.Parse(line[3])], //1
                                model.NodesDictionary[Int32.Parse(line[2])], //2
                                model.NodesDictionary[Int32.Parse(line[0])]  //3
                            });
                        }
                        i = i + 3;
                        for (int QuadID = 0; QuadID < NumberOfQuadElements; QuadID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);
                            int boundaryID = Int32.Parse(line[0]);
                            bool res = Array.Exists(boundaryIDs, x => x == boundaryID - 1);
                            if (res)
                            {
                                quadBoundaries[boundaryID].Add(quadNodesCollection[QuadID]);
                                foreach (Node node in quadNodesCollection[QuadID])
                                {
                                    nodeBoundaries[boundaryID].Add(node);
                                }
                            }
                        }

                        foreach (int boundaryID in boundaryIDs)
                        {
                            nodeBoundaries[boundaryID] = nodeBoundaries[boundaryID].Distinct().ToList();
                        }
                        break;
                    case Attributes.tet:
                        do
                        {
                            i++;
                            line = text[i].Split(delimeters);
                        } while (line[0] == "");
                        i++;
                        line = text[i].Split(delimeters);
                        NumberOfTetElements = Int32.Parse(line[0]);
                        i++;
                        List<IReadOnlyList<Node>> TetraNodes = new List<IReadOnlyList<Node>>();
                        for (int TetID = 0; TetID < NumberOfTetElements; TetID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);

                            TetraNodes.Add(new List<Node>
                            {
                                model.NodesDictionary[Int32.Parse(line[3])],
                                model.NodesDictionary[Int32.Parse(line[2])],
                                model.NodesDictionary[Int32.Parse(line[1])],
                                model.NodesDictionary[Int32.Parse(line[0])]
                            });
                        }
                        i = i + 3;
                        for (int TetID = 0; TetID < NumberOfTetElements; TetID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);
                            int elementDomainID = Int32.Parse(line[0]);
                            bool res = Array.Exists(domainIDs, x => x == elementDomainID - 1);
                            if (res)
                            {
                                IReadOnlyList<Node> nodesTet = TetraNodes[TetID];
                                var Tet4 = elementFactory3D[elementDomainID - 1].CreateElement(CellType.Tet4, nodesTet);
                                var element = new Element();
                                element.ID = TetID;
                                element.ElementType = Tet4;
                                foreach (Node node in nodesTet)
                                {
                                    element.AddNode(node);
                                    nodeDomains[elementDomainID - 1].Add(node);
                                }
                                model.SubdomainsDictionary[0].Elements.Add(element);
                                model.ElementsDictionary.Add(TetID, element);
                                elementDomains[elementDomainID - 1].Add(model.ElementsDictionary[TetID]);
                            }
                        }
                        foreach (int domainID in domainIDs)
                        {
                            nodeDomains[domainID] = nodeDomains[domainID].Distinct().ToList();
                        }
                        break;
                    case Attributes.hex:
                        do
                        {
                            i++;
                            line = text[i].Split(delimeters);
                        } while (line[0] == "");
                        i++;
                        line = text[i].Split(delimeters);
                        NumberOfHexElements = Int32.Parse(line[0]);
                        i++;
                        List<IReadOnlyList<Node>> HexaNodes = new List<IReadOnlyList<Node>>();
                        for (int HexID = 0; HexID < NumberOfHexElements; HexID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);

                            HexaNodes.Add(new List<Node>
                            {
                                model.NodesDictionary[Int32.Parse(line[4])],
                                model.NodesDictionary[Int32.Parse(line[6])],
                                model.NodesDictionary[Int32.Parse(line[7])],
                                model.NodesDictionary[Int32.Parse(line[5])],
                                model.NodesDictionary[Int32.Parse(line[0])],
                                model.NodesDictionary[Int32.Parse(line[2])],
                                model.NodesDictionary[Int32.Parse(line[3])],
                                model.NodesDictionary[Int32.Parse(line[1])],
                            });
                        }
                        i = i + 3;
                        for (int HexID = 0; HexID < NumberOfHexElements; HexID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);
                            int elementDomainID = Int32.Parse(line[0]);
                            bool res = Array.Exists(domainIDs, x => x == elementDomainID - 1);
                            if (res)
                            {
                                IReadOnlyList<Node> nodesHex = HexaNodes[HexID];
                                var Hexa8 = elementFactory3D[elementDomainID - 1].CreateElement(CellType.Hexa8, nodesHex);
                                var element = new Element();
                                element.ID = HexID;
                                element.ElementType = Hexa8;
                                foreach (Node node in nodesHex)
                                {
                                    element.AddNode(node);
                                    nodeDomains[elementDomainID - 1].Add(node);
                                }
                                model.SubdomainsDictionary[0].Elements.Add(element);
                                model.ElementsDictionary.Add(HexID, element);
                                elementDomains[elementDomainID - 1].Add(model.ElementsDictionary[HexID]);
                            }
                        }
                        foreach (int domainID in domainIDs)
                        {
                            nodeDomains[domainID] = nodeDomains[domainID].Distinct().ToList();
                        }
                        break;
                }
            }
            IList<Node> modelNodes = nodeDomains.SelectMany(x => x).ToList();
            IList<int> nodesToRemove = new List<int>();
            foreach (Node node in model.NodesDictionary.Values)
            {
                if (!(modelNodes.Contains(node))) nodesToRemove.Add(node.ID);
            }
            foreach (int nodeID in nodesToRemove)
            {
                model.NodesDictionary.Remove(nodeID);
            }
            this.Model = model;
            return model;
        }

    }
}
