﻿using ISAAR.MSolve.FEM.Entities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Embedding
{
    public class EmbeddedNode
    {
        private readonly Node node;
        private readonly Element embeddedInElement;
        private readonly IList<double> coordinates = new List<double>();
        private readonly IList<DOFType> dependentDOFs;

        public Node Node { get { return node; } }
        public Element EmbeddedInElement { get { return embeddedInElement; } }
        public IList<DOFType> DependentDOFs { get { return dependentDOFs; } }
        public IList<double> Coordinates { get { return coordinates; } }

        public EmbeddedNode(Node node, Element hostElement, IList<DOFType> dependentDOFs)
        {
            this.node = node;
            this.embeddedInElement = hostElement;
            this.dependentDOFs = dependentDOFs;
        }

        public override bool Equals(object obj)
        {
            EmbeddedNode e = obj as EmbeddedNode;
            if (e == null) return false;

            return (e.Node == this.node && e.EmbeddedInElement == embeddedInElement);
        }
            
        public override string ToString()
        {
            if (node == null) return "(Null inside node)";
            string hostElementID = embeddedInElement == null ? "(N/A)" : embeddedInElement.ID.ToString();
            string coordinateString = String.Empty;
            foreach (var c in coordinates)
                coordinateString += c.ToString() + " ";

            var header = String.Format("{0} (host: {4}): ({1}, {2}, {3}) [{5}]", node.ID, node.X, node.Y, node.Z, hostElementID, coordinateString);
            string constraintsDescripton = string.Empty;
            foreach (var c in node.Constraints)
            {
                string con = string.Empty;
                switch (c)
                {
                    case DOFType.Pore:
                        con = "Pore";
                        break;
                    case DOFType.RotX:
                        con = "rX";
                        break;
                    case DOFType.RotY:
                        con = "rY";
                        break;
                    case DOFType.RotZ:
                        con = "rZ";
                        break;
                    case DOFType.Unknown:
                        con = "?";
                        break;
                    case DOFType.X:
                        con = "X";
                        break;
                    case DOFType.Y:
                        con = "Y";
                        break;
                    case DOFType.Z:
                        con = "Z";
                        break;
                }
                constraintsDescripton += c.ToString() + ", ";
            }
            constraintsDescripton = constraintsDescripton.Length > 1 ? constraintsDescripton.Substring(0, constraintsDescripton.Length - 2) : constraintsDescripton;

            return String.Format("{0} - Con ({1})", header, constraintsDescripton);
        }

    }
}
