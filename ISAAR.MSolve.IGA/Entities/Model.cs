﻿using System;
using System.Collections.Generic;
using System.Linq;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Vectors;

//TODO: find what is going on with the dynamic loads and refactor them. That 564000000 in AssignMassAccelerationHistoryLoads()
//      cannot be correct.
namespace ISAAR.MSolve.IGA.Entities
{
	public class Model : IStructuralModel
	{
		private int numberOfPatches = 0;
		private int numberOfInterfaces = 0;

		private readonly Dictionary<int, ControlPoint> controlPointsDictionary = new Dictionary<int, ControlPoint>();
		private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();
		private readonly Dictionary<int, Patch> patchesDictionary = new Dictionary<int, Patch>();

		private readonly IList<Load> loads = new List<Load>();

		private readonly IList<IMassAccelerationHistoryLoad> massAccelerationHistoryLoads =
			new List<IMassAccelerationHistoryLoad>();

		private IGlobalFreeDofOrdering globalDofOrdering;

		//public IList<EmbeddedNode> EmbeddedNodes { get; } = new List<EmbeddedNode>();

		public IList<Cluster> Clusters => ClustersDictionary.Values.ToList();
		public Dictionary<int, Cluster> ClustersDictionary { get; } = new Dictionary<int, Cluster>();

		IReadOnlyList<IElement> IStructuralModel.Elements => ElementsDictionary.Values.ToList();
		public IList<Element> Elements => ElementsDictionary.Values.ToList();
		public Dictionary<int, Element> ElementsDictionary => elementsDictionary;

		public IList<Load> Loads { get; private set; } = new List<Load>();

		public IList<IMassAccelerationHistoryLoad> MassAccelerationHistoryLoads { get; } =
			new List<IMassAccelerationHistoryLoad>();

		public IList<ControlPoint> ControlPoints => controlPointsDictionary.Values.ToList();
		IReadOnlyList<INode> IStructuralModel.Nodes => controlPointsDictionary.Values.ToList();

		public Dictionary<int, ControlPoint> ControlPointsDictionary
		{
			get => controlPointsDictionary;
		}

		public int NumberOfPatches
		{
			get { return numberOfPatches; }
			set { numberOfPatches = value; }
		}

		public int NumberOfInterfaces
		{
			get { return numberOfInterfaces; }
			set { numberOfInterfaces = value; }
		}

		IReadOnlyList<ISubdomain> IStructuralModel.Subdomains => patchesDictionary.Values.ToList();
		public IList<Patch> Patches => patchesDictionary.Values.ToList();

		public Dictionary<int, Patch> PatchesDictionary
		{
			get => patchesDictionary;
		}

		public Table<INode, IDofType, double> Constraints { get; private set; } =
			new Table<INode, IDofType, double>(); //TODOMaria: maybe it's useless in model class

		public IGlobalFreeDofOrdering GlobalDofOrdering
		{
			get => globalDofOrdering;
			set
			{
				globalDofOrdering = value;
				foreach (Patch patch in Patches)
				{
					patch.FreeDofRowOrdering = GlobalDofOrdering.SubdomainDofOrderings[patch];
				}

				//EnumerateSubdomainLagranges();
				//EnumerateDOFMultiplicity();
			}
		}

        private IGlobalFreeDofOrdering globalRowDofOrdering;
        public IGlobalFreeDofOrdering GlobalRowDofOrdering 
        {
            get => globalRowDofOrdering;
            set
            {
                globalRowDofOrdering = value;
                foreach (Patch patch in Patches)
                {
                    patch.FreeDofColOrdering = GlobalRowDofOrdering.SubdomainDofOrderings[patch];
                    patch.Forces = Vector.CreateZero(patch.FreeDofColOrdering.NumFreeDofs);
                }

                //EnumerateSubdomainLagranges();
                //EnumerateDOFMultiplicity();
            }
        }

        private IGlobalFreeDofOrdering globalColDofOrdering;
        public IGlobalFreeDofOrdering GlobalColDofOrdering
        {
            get => globalColDofOrdering;
            set
            {
                globalColDofOrdering = value;
                foreach (Patch patch in Patches)
                {
                    patch.FreeDofColOrdering = GlobalColDofOrdering.SubdomainDofOrderings[patch];
                    patch.Forces = Vector.CreateZero(patch.FreeDofColOrdering.NumFreeDofs);
                }

                //EnumerateSubdomainLagranges();
                //EnumerateDOFMultiplicity();
            }
        }
        
        public void AssignLoads(NodalLoadsToSubdomainsDistributor distributeNodalLoads)
		{
			foreach (Patch patch in PatchesDictionary.Values) patch.Forces.Clear();
			AssignControlPointLoads(distributeNodalLoads);
			AssignBoundaryLoads();
		}

		public void AssignMassAccelerationHistoryLoads(int timeStep)
		{
			throw new NotImplementedException();
		}

		private void AssignBoundaryLoads()
		{
			foreach (Patch patch in patchesDictionary.Values)
			{
				foreach (Edge edge in patch.EdgesDictionary.Values)
				{
					Dictionary<int, double> edgeLoadDictionary = edge.CalculateLoads();
					foreach (int dof in edgeLoadDictionary.Keys)
						if (dof != -1)
							patch.Forces[dof] += edgeLoadDictionary[dof];
				}

				foreach (Face face in patch.FacesDictionary.Values)
				{
					Dictionary<int, double> faceLoadDictionary = face.CalculateLoads();
					foreach (int dof in faceLoadDictionary.Keys)
						if (dof != -1)
							patch.Forces[dof] += faceLoadDictionary[dof];
				}
			}
		}

		private void AssignControlPointLoads(NodalLoadsToSubdomainsDistributor distributeControlPointLoads)
		{
            var globalPointLoads = new Table<INode, IDofType, double>();
            foreach (Load load in Loads) globalPointLoads.TryAdd(load.ControlPoint, load.DOF, load.Amount);

            Dictionary<int, SparseVector> patchPointLoads = distributeControlPointLoads(globalPointLoads);
            foreach (var idPatchLoads in patchPointLoads)
            {
                PatchesDictionary[idPatchLoads.Key].Forces.AddIntoThis(idPatchLoads.Value);
            }
      

            // Old code. It should probably be deleted.
            //foreach (Patch patch in PatchesDictionary.Values)
			//{
			//	patch.ControlPointLoads = new Table<ControlPoint, DOFType, double>();
			//}

			//foreach (Load load in Loads)
			//{
			//	var cp = ((ControlPoint) load.ControlPoint);
			//	double amountPerPatch = load.Amount / cp.PatchesDictionary.Count;
			//	foreach (Patch patch in cp.PatchesDictionary.Values)
			//	{
			//		bool wasNotContained = patch.ControlPointLoads.TryAdd(cp, load.DOF, amountPerPatch);
			//	}
			//}

			////TODO: this should be done by the subdomain when the analyzer decides.
			//foreach (Patch patch in PatchesDictionary.Values)
			//{
			//	foreach ((ControlPoint node, DOFType dofType, double amount) in patch.ControlPointLoads)
			//	{
			//		if (!patch.DofOrdering.FreeDofs.Contains(node, dofType)) continue;
			//		int patchDofIdx = patch.DofOrdering.FreeDofs[node, dofType];
			//		patch.Forces[patchDofIdx] = amount;
			//	}
			//}
		}


		//What is the purpose of this method? If someone wanted to clear the Model, they could just create a new one.
		public void Clear()
		{
			Loads.Clear();
			ClustersDictionary.Clear();
			PatchesDictionary.Clear();
			ElementsDictionary.Clear();
			ControlPointsDictionary.Clear();
			globalDofOrdering = null;
			Constraints.Clear();
			MassAccelerationHistoryLoads.Clear();
		}

		public void ConnectDataStructures()
		{
			BuildInterconnectionData();
			AssignConstraints();
            RemoveInactiveNodalLoads();
        }

		//TODO: constraints should not be saved inside the nodes. As it is right now (22/11/2018) the same constraint 
		//      is saved in the node, the model constraints table and the subdomain constraints table. Furthermore,
		//      displacement control analyzer updates the subdomain constraints table only (another bad design decision).  
		//      It is too easy to access the wrong instance of the constraint. 
		private void AssignConstraints()
		{
			foreach (ControlPoint controlPoint in ControlPointsDictionary.Values)
			{
				if (controlPoint.Constraints == null) continue;
				foreach (Constraint constraint in controlPoint.Constraints)
					Constraints[controlPoint, constraint.DOF] = constraint.Amount;
			}

			foreach (Patch patch in PatchesDictionary.Values) patch.ExtractConstraintsFromGlobal(Constraints);
		}

		private void BuildElementDictionaryOfEachControlPoint()
		{
			foreach (Element element in elementsDictionary.Values)
			foreach (ControlPoint controlPoint in element.ControlPoints)
				controlPoint.ElementsDictionary.Add(element.ID, element);
		}

		private void
			BuildInterconnectionData() //TODOMaria: maybe I have to generate the constraints dictionary for each subdomain here
		{
			BuildPatchOfEachElement();
			BuildElementDictionaryOfEachControlPoint();
			foreach (ControlPoint controlPoint in ControlPointsDictionary.Values) controlPoint.BuildPatchesDictionary();

			//foreach (Patch patch in PatchesDictionary.Values) patch.DefineControlPointsFromElements();
		}

		private void BuildPatchOfEachElement()
		{
			foreach (Patch patch in patchesDictionary.Values)
			foreach (Element element in patch.Elements)
			{
				element.Patch = patch;
				element.Model = this;
			}
		}

        private void RemoveInactiveNodalLoads()
        {
            var activeLoads = new List<Load>(Loads.Count);
            foreach (Load load in Loads)
            {
                bool isConstrained = Constraints.Contains(load.ControlPoint, load.DOF);
                if (!isConstrained) activeLoads.Add(load);
            }
            Loads = activeLoads;
        }
    }
}