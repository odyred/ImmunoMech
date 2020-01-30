using System;
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

namespace ISAAR.MSolve.FEM.Readers
{
    public class ComsolMeshReader
    {
        enum Attributes
        {
            sdim = 1001,
            Meshpointcoordinates = 1002,
            numberofmeshpoints = 1003,
            numberofelementtypes = 1004,
            //typename = 1005,
            vtx = 1006,
            edg = 1007,
            tri = 1008,
            tet = 1009,
            //numberofnodesperelement = 1010,
            //numberofelements = 1011,
            //Elements = 1012
        }

        public string Filename { get; private set; }

        public ComsolMeshReader(string filename)
        {
            Filename = filename;
        }

        private static readonly ElasticMaterial3D Material0 = new ElasticMaterial3D
        {
            YoungModulus = 1000,
            PoissonRatio = 0.3
        };

        private static readonly DynamicMaterial DynamicMaterial = new DynamicMaterial(1, 0, 0);



        int NumberOfNodes;
        int NumberOfVtxElements;
        int NumberOfEdgElements;
        int NumberOfTriElements;
        int NumberOfElements;

        public void CreateModelFromFile()
        {
            double density = 1.0;
            double k = 1.0;
            double c = 1.0;
            var elementFactory = new ThermalElement3DFactory(new ThermalMaterial(density, c, k));
            var elements = new ThermalElement3D[4];
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(0);
            // Material

            char[] delimeters = { ' ', '=', '\t' };
            Attributes? name = null;

            String[] text = System.IO.File.ReadAllLines(Filename);

            for (int i = 0; i < text.Length; i++)
            {
                String[] line = text[i].Split(delimeters, StringSplitOptions.RemoveEmptyEntries);
                StringBuilder comparisonString = new StringBuilder(null, 50);
                if (line.Length == 0)
                {
                    continue;
                }
                if (line[0] == "3" & line.Length > 1)
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
                        //IList<Node> boundary1 = new List<Node> { null };
                        //IList<Node> boundary2 = new List<Node> { null };
                        //IList<Node> boundary3 = new List<Node> { null };
                        for (int j = 0; j < NumberOfNodes; j++)
                        {
                            i++;
                            line = text[i].Split(delimeters);
                            int nodeGlobalID = j;
                            double x = Double.Parse(line[0], CultureInfo.InvariantCulture);
                            double y = Double.Parse(line[1], CultureInfo.InvariantCulture);
                            double z = Double.Parse(line[2], CultureInfo.InvariantCulture);
                            Node node = new Node (nodeGlobalID, x, y, z);
                            model.NodesDictionary.Add(nodeGlobalID, node);
                            nodelist.Add(node);
                        }

                        //foreach (Node node in nodelist)
                        //{
                        //    if (node.X == 0)
                        //        model.NodesDictionary[node.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0});
                        //    if (node.Y == 0)
                        //        model.NodesDictionary[node.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0 });
                        //    if (node.Z == 0)
                        //        model.NodesDictionary[node.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0 });
                        //}
                        break;
                    //case Attributes.vtx:
                    //    do
                    //    {
                    //        i++;
                    //        line = text[i].Split(delimeters);
                    //    } while (line[0] == "");
                    //    i++;
                    //    line = text[i].Split(delimeters);
                    //    NumberOfTetElements = Int32.Parse(line[0]);
                    //    i++;
                    //    for (int VtxID = 0; VtxID < NumberOfTetElements; VtxID++)
                    //    {
                    //        i++;
                    //        line = text[i].Split(delimeters);

                    //        IReadOnlyList<Node> nodelist = new List<Node>
                    //        {
                    //            model.NodesDictionary[Int32.Parse(line[0])],
                    //        };
                    //        var factory = new ContinuumElement3DFactory(Material0, DynamicMaterial);
                    //        var Tet4 = factory.CreateElement(CellType.Tet4, nodelist);
                    //        var element = new Element();
                    //        element.ID = VtxID;
                    //        element.ElementType = Tet4;
                    //        foreach (Node node in nodelist)
                    //        {
                    //            element.AddNode(node);
                    //        }
                    //        model.SubdomainsDictionary[0].Elements.Add(element);
                    //        model.ElementsDictionary.Add(VtxID, element);
                    //    }
                    //    break;
                    //case Attributes.edg:
                    //    do
                    //    {
                    //        i++;
                    //        line = text[i].Split(delimeters);
                    //    } while (line[0] == "");
                    //    i++;
                    //    line = text[i].Split(delimeters);
                    //    NumberOfEdgElements = Int32.Parse(line[0]);
                    //    i++;
                    //    for (int EdgID = 0; EdgID < NumberOfEdgElements; EdgID++)
                    //    {
                    //        i++;
                    //        line = text[i].Split(delimeters);

                    //        IReadOnlyList<Node> nodesCollection = new List<Node>
                    //        {
                    //            model.NodesDictionary[Int32.Parse(line[0])],
                    //            model.NodesDictionary[Int32.Parse(line[1])]
                    //        };
                    //        var factory = new ContinuumElement3DFactory(Material0, DynamicMaterial);
                    //        var Line = factory.CreateElement(CellType.Line, nodesCollection);
                    //        var element = new Element();
                    //        element.ID = NumberOfVtxElements + EdgID;
                    //        element.ElementType = Line;
                    //        foreach (Node node in nodesCollection)
                    //        {
                    //            element.AddNode(node);
                    //        }
                    //        model.SubdomainsDictionary[0].Elements.Add(element);
                    //        model.ElementsDictionary.Add(EdgID, element);
                    //    }
                    //    break;
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
                        for (int TriID = 0; TriID < NumberOfTriElements; TriID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);

                            IReadOnlyList<Node> nodesCollection = new List<Node>
                            {
                                model.NodesDictionary[Int32.Parse(line[0])],
                                model.NodesDictionary[Int32.Parse(line[1])],
                                model.NodesDictionary[Int32.Parse(line[2])]
                            };

                            foreach (Node node in nodesCollection)
                            {
                                if (node.X < 5e-4 | node.Y < 5e-4 | node.Z < 5e-4)
                                    model.NodesDictionary[node.ID].Constraints.Add(new Constraint() { DOF = ThermalDof.Temperature, Amount = 0 });
                            }
                            //var factory = new ContinuumElement3DFactory(Material0, DynamicMaterial);
                            //var Tri3 = factory.CreateElement(CellType.Tri3, nodesCollection);
                            //var element = new Element();
                            //element.ID = NumberOfVtxElements + NumberOfEdgElements + TriID;
                            //element.ElementType = Tri3;
                            //foreach (Node node in nodesCollection)
                            //{
                            //    element.AddNode(node);
                            //}
                            //model.SubdomainsDictionary[0].Elements.Add(element);
                            //model.ElementsDictionary.Add(TriID, element);
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
                        NumberOfElements = Int32.Parse(line[0]);
                        i++;
                        for (int TetID = 0; TetID < NumberOfElements; TetID++)
                        {
                            i++;
                            line = text[i].Split(delimeters);

                            IReadOnlyList<Node> nodesCollection = new List<Node>
                            {
                                model.NodesDictionary[Int32.Parse(line[0])],
                                model.NodesDictionary[Int32.Parse(line[1])],
                                model.NodesDictionary[Int32.Parse(line[2])],
                                model.NodesDictionary[Int32.Parse(line[3])]
                            };
//                            var factory = new ContinuumElement3DFactory(Material0, DynamicMaterial);
                            var Tet4 = elementFactory.CreateElement(CellType.Tet4, nodesCollection);
                            var element = new Element();
                            element.ID = TetID;
                            element.ElementType = Tet4;
                            foreach (Node node in nodesCollection)
                            {
                                element.AddNode(node);
                            }
                            model.SubdomainsDictionary[0].Elements.Add(element);
                            model.ElementsDictionary.Add(TetID, element);
                        }
                        break;
                }
            }
        }
    }
}
