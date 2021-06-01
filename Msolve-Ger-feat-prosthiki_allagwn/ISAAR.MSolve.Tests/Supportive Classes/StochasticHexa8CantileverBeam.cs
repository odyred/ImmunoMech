using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;


namespace ISAAR.MSolve.Tests
{
    public static class StochasticHexaSimpleCantileverBeam
    {
        public static void MakeCantileverBeam(Model model, double startX, double startY, double startZ, int startNodeID, int startElementID, int subdomainID, IStochasticMaterialCoefficientsProvider coefficientsProvider)

        {

            int nodeID = startNodeID;

            for (int j = 0; j < 4; j++)
            {
                if (nodeID % 2 == 0)
                {
                    model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY + 0.25, Z = startZ + 0.25 * (j / 2) });
                }
                else
                {
                    model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX, Y = startY, Z = startZ + 0.25 * (j / 2) });
                }
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.X);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Y);
                model.NodesDictionary[nodeID].Constraints.Add(DOFType.Z);

                nodeID++;
            }

            for (int i = 0; i < 4; i++)
            {
                for (int k = 0; k < 4; k++)
                {
                    if (nodeID % 2 == 0)
                    {
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + 0.25 * (i + 1), Y = startY + 0.25, Z = startZ + 0.25 * (k / 2) });
                    }
                    else
                    {
                        model.NodesDictionary.Add(nodeID, new Node() { ID = nodeID, X = startX + 0.25 * (i + 1), Y = startY, Z = startZ + 0.25 * (k / 2) });
                    }
                    nodeID++;
                }
            }

            int elementID = startElementID;
            Element e;
            StochasticElasticMaterial3D material = new StochasticElasticMaterial3D(coefficientsProvider)
            {
                YoungModulus = 2.0e7,
                PoissonRatio = 0.3
            };

            var conMat = material.ConstitutiveMatrix;

            for (int i = 0; i < 4; i++)
            {

                e = new Element()
                {
                    ID = elementID,
                    ElementType = new Hexa8WithStochasticMaterial(material)

                };

                e.NodesDictionary.Add(startNodeID + 4 * i, model.NodesDictionary[startNodeID + 4 * i]);
                e.NodesDictionary.Add(startNodeID + 4 * i + 4, model.NodesDictionary[startNodeID + 4 * i + 4]);
                e.NodesDictionary.Add(startNodeID + 4 * i + 5, model.NodesDictionary[startNodeID + 4 * i + 5]);
                e.NodesDictionary.Add(startNodeID + 4 * i + 1, model.NodesDictionary[startNodeID + 4 * i + 1]);

                e.NodesDictionary.Add(startNodeID + 4 * i + 2, model.NodesDictionary[startNodeID + 4 * i + 2]);
                e.NodesDictionary.Add(startNodeID + 4 * i + 6, model.NodesDictionary[startNodeID + 4 * i + 6]);
                e.NodesDictionary.Add(startNodeID + 4 * i + 7, model.NodesDictionary[startNodeID + 4 * i + 7]);
                e.NodesDictionary.Add(startNodeID + 4 * i + 3, model.NodesDictionary[startNodeID + 4 * i + 3]);

                model.ElementsDictionary.Add(e.ID, e);
                model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e.ID, e);


                elementID++;
            }
        }


    }

}