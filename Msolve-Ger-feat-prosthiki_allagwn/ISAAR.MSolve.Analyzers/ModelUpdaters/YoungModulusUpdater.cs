using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Analyzers.Interfaces;
using ISAAR.MSolve.FEM.Elements;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.Analyzers.ModelUpdaters
{
    public class YoungModulusUpdater : IModelUpdater
    {
        private Model model;
        private double youngModulus;
        private IStochasticCoefficientsProvider coefficientsProvider;

        IList<IFiniteElement> listofFiniteElements;
        public YoungModulusUpdater(Model model, double youngModulus, IStochasticCoefficientsProvider coefficientsProvider)
        {
            this.model = model;
            this.youngModulus = youngModulus;
            this.coefficientsProvider = coefficientsProvider;
        }

        public void UpdateModel()
        {
            throw new NotImplementedException();
        }
    }
}

