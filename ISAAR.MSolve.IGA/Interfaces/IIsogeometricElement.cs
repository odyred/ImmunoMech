﻿using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Entities.Loads;

namespace ISAAR.MSolve.IGA.Interfaces
{

    public interface IIsogeometricElement: IElementType
	{
        int ID { get; }
        ElementDimensions ElementDimensions { get; }
        IElementDofEnumerator DofEnumerator { get; set; }
        //IList<IList<DOFType>> GetElementDOFTypes(IElement element);
        bool MaterialModified { get; }
        //IMatrix2D StiffnessMatrix(Element element);
        //IMatrix2D MassMatrix(Element element);
        //IMatrix2D DampingMatrix(Element element);
        Dictionary<int, double> CalculateLoadingCondition(Element element,Edge edge, NeumannBoundaryCondition neumann);
        Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann);
        Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure);
        Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure);
        void ResetMaterialModified();
        //Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements);
        //double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements);
        //double[] CalculateForcesForLogging(Element element, double[] localDisplacements);
		double[,] CalculateDisplacementsForPostProcessing(Element element, double[,] localDisplacements);
        //double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads);
        //void SaveMaterialState();
        void ClearMaterialState();

        //void ClearMaterialStresses();
    }
}
