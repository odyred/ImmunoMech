using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;
namespace ISAAR.MSolve.Analyzers.Interfaces
{

    public interface IStochasticModelBuilder
    {
        //IStochasticCoefficientsProvider coefficientsProvider { get; set; }
        void BuildStochasticModel();
        void UpdateStochasticCoefficientsProvider(IStochasticCoefficientsProvider coefficientsProvider);
    }
}
