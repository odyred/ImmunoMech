﻿using System.Collections.Generic;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;

namespace ISAAR.MSolve.FEM
{
    public class GenericDOFEnumerator : IFiniteElementDOFEnumerator
    {
        public IList<IList<DOFType>> GetDOFTypes(Element element)
        {
            return element.ElementType.GetElementDOFTypes(element);
        }

        public IList<IList<DOFType>> GetDOFTypesForDOFEnumeration(Element element)
        {
            return element.ElementType.GetElementDOFTypes(element);
        }

        public IMatrix2D GetTransformedMatrix(IMatrix2D matrix)
        {
            return matrix;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }

        public double[] GetTransformedVector(double[] vector)
        {
            return vector;
        }

        public double[] GetTransformedForcesVector(double[] vector)
        {
            return vector;
        }
    }
}
