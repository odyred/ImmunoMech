using ISAAR.MSolve.FEM.Entities;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Readers.Interfaces
{
    public interface IModelReader
    {
        Model Model { get; }
        IList<IList<IList<Node>>> quadBoundaries { get; }
        IList<IList<IList<Node>>> triBoundaries { get; }
        IList<IList<Element>> elementDomains { get; }
        IList<IList<Node>> nodeDomains { get; }
    }
}
