using ISAAR.MSolve.PreProcessor;
//using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.Materials;//using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.FEM.Materials;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;//using ISAAR.MSolve.Matrices.Interfaces;
using System.Runtime.InteropServices;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
// compa
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Elements;

namespace ISAAR.MSolve.SamplesConsole
{
    class DddmExamplesBuilder
    {
        public static void Reference1RVEExample10000(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 1;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp, Dq);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters(model, gp, ekk_xyz[j], model_o_x_parameteroi[j]);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\Reference_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF1_10000\Fxk_p_komvoi_rve.txt");
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa3);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void Reference1RVEExample10000_Hexaonly(Model model)
        {
            double[,] Dq = new double[1, 1];
            Tuple<rveMatrixParameters, grapheneSheetParameters> mpgp;
            rveMatrixParameters mp;
            grapheneSheetParameters gp;
            int subdiscr1 = 4;
            int discr1 = 4;
            // int discr2 dn xrhsimopoieitai
            int discr3 = 10;
            int subdiscr1_shell = 7;
            int discr1_shell = 1;
            mpgp = RVEExamplesBuilder.GetReferenceRveExampleParameters(subdiscr1, discr1, discr3, subdiscr1_shell, discr1_shell);
            mp = mpgp.Item1;
            gp = mpgp.Item2;
            double[][] ekk_xyz = new double[2][] { new double[] { 0, 0, 0 }, new double[] { 0.25 * 105, 0, 0.25 * 40 } };

            int graphene_sheets_number = 0;
            o_x_parameters[] model_o_x_parameteroi = new o_x_parameters[graphene_sheets_number];

            RVEExamplesBuilder.HexaElementsOnlyRVE(model, mp, Dq);
            int hexaElementsNumber = model.ElementsDictionary.Count();

            IEnumerable<Element> hostGroup = model.ElementsDictionary.Where(x => (x.Key < hexaElementsNumber + 1)).Select(kv => kv.Value);
            List<int> EmbeddedElementsIDs = new List<int>();
            int element_counter_after_Adding_sheet;
            element_counter_after_Adding_sheet = hexaElementsNumber; // initial value before adding first graphene sheet
            int shellElementsNumber;

            for (int j = 0; j < graphene_sheets_number; j++)
            {
                RVEExamplesBuilder.AddGrapheneSheet_with_o_x_parameters(model, gp, ekk_xyz[j], model_o_x_parameteroi[j]);
                shellElementsNumber = (model.ElementsDictionary.Count() - element_counter_after_Adding_sheet) / 3; //tha xrhsimefsei
                //embdeddedGroup_adittion= model.ElementsDictionary.Where(x => (x.Key >= shellElementsNumber + element_counter_after_Adding_sheet + 1)).Select(kv => kv.Value);
                //embdeddedGroup.Concat(embdeddedGroup_adittion);
                for (int k = shellElementsNumber + element_counter_after_Adding_sheet + 1; k < model.ElementsDictionary.Count() + 1; k++)
                {
                    EmbeddedElementsIDs.Add(model.ElementsDictionary[k].ID);
                }
                element_counter_after_Adding_sheet = model.ElementsDictionary.Count();

            }

            // model: add loads
            RVEExamplesBuilder.AddLoadsOnRveFromFile(model, mp.hexa1, mp.hexa2, mp.hexa3, @"C:\Users\turbo-x\Desktop\notes_elegxoi\Reference_Examples_Dokimi\fe2_tax_me1_arxiko_chol_dixws_me1_OriginalRVEExampleChol_me_a1_REF1_10000\Fxk_p_komvoi_rve.txt");
            //RVEExamplesBuilder.AddXLoadsOnYZConstrainedOneElementRVE(model);
            // model: add constraints
            RVEExamplesBuilder.AddConstraintsForNonSingularStiffnessMatrix(model, mp.hexa1, mp.hexa2, mp.hexa3);
            //RVEExamplesBuilder.AddConstraintsForYZConstraindeOneElementRVE(model);

            int[] EmbElementsIds = EmbeddedElementsIDs.ToArray();
            IEnumerable<Element> embdeddedGroup = model.ElementsDictionary.Where(x => (Array.IndexOf(EmbElementsIds, x.Key) > -1)).Select(kv => kv.Value); // dld einai null afth th stigmh
            var embeddedGrouping = new EmbeddedCohesiveGrouping(model, hostGroup, embdeddedGroup);
        }

        public static void MakeModelDictionariesZeroBasedForDecomposer(Model model)
        {
            model.SubdomainsDictionary[1].ID = 0;
            Subdomain subdomain_ini = model.SubdomainsDictionary[1];
            model.SubdomainsDictionary.Remove(1);
            model.SubdomainsDictionary.Add(0, subdomain_ini);

            Dictionary<int, Element> ElementsDictionary_2 = new Dictionary<int, Element>(model.ElementsDictionary.Count);
            for (int i1 = 0; i1 < model.ElementsDictionary.Count; i1++)
            {
                ElementsDictionary_2.Add(model.ElementsDictionary[i1 + 1].ID - 1, model.ElementsDictionary[i1 + 1]);
                ElementsDictionary_2[model.ElementsDictionary[i1 + 1].ID - 1].ID += -1;
            }
            int nElement = model.ElementsDictionary.Count;
            for (int i1 = 0; i1 < nElement; i1++)
            {
                model.ElementsDictionary.Remove(i1 + 1);
                model.SubdomainsDictionary[0].ElementsDictionary.Remove(i1 + 1);
            }
            for (int i1 = 0; i1 < nElement; i1++)
            {
                model.ElementsDictionary.Add(ElementsDictionary_2[i1].ID, ElementsDictionary_2[i1]);
                model.SubdomainsDictionary[0].ElementsDictionary.Add(ElementsDictionary_2[i1].ID, ElementsDictionary_2[i1]);
            }

            Dictionary<int, Node> NodesDictionary_2 = new Dictionary<int, Node>(model.NodesDictionary.Count);
            for (int i1 = 0; i1 < model.NodesDictionary.Count; i1++)
            {
                NodesDictionary_2.Add(model.NodesDictionary[i1 + 1].ID - 1, model.NodesDictionary[i1 + 1]);
                NodesDictionary_2[model.NodesDictionary[i1 + 1].ID - 1].ID += -1;
            }

            int nNode = model.NodesDictionary.Count;
            for (int i1 = 0; i1 < nNode; i1++)
            {
                model.NodesDictionary.Remove(i1 + 1);
                //model.SubdomainsDictionary[0].NodesDictionary.Remove(i1 + 1); // den peirazoume to subdomain nodesDictionary, ftiahnetai mono tou (pithanws apo to connect data Structures)
            }
            for (int i1 = 0; i1 < nNode; i1++)
            {
                model.NodesDictionary.Add(NodesDictionary_2[i1].ID, NodesDictionary_2[i1]);
                //model.SubdomainsDictionary[0].NodesDictionary.Add(NodesDictionary_2[i1].ID, NodesDictionary_2[i1]); // den peirazoume to subdomain nodesDictionary, ftiahnetai mono tou (pithanws apo to connect data Structures)
            }
        }
    }
}
