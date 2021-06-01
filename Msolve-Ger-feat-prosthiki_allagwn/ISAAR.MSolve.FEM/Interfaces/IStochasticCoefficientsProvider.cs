using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Interfaces
{
    public interface IStochasticCoefficientsProvider
    {
        double[] RandomVariables { get; set; }
        double GetCoefficient(double meanValue, double[] coordinates);
    }
}
