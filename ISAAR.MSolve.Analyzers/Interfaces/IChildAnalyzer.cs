﻿using ISAAR.MSolve.LinearAlgebra.Vectors;
using System.Collections.Generic;

namespace ISAAR.MSolve.Analyzers.Interfaces
{
    public interface IChildAnalyzer: IAnalyzer
    {
        //TODO: the setter could be avoided by passing a child analyzer builder in the parent analyzer's constructor and the
        //      parent analyzer to the child analyzers constructor, before the parent analyzer constructor ends. However that is
        //      too much trouble without a specific need for that degree of immutability.
        IParentAnalyzer ParentAnalyzer { get; set; }
        Dictionary<int, IVector> Responses { get; set; }
    }
}
