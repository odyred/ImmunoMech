using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials
{
    public class ConvectionDiffusionMaterial
    {
        public ConvectionDiffusionMaterial(double diffusionCoeff, double[] convectionCoeff, double absorptionRate)
        {
            this.DiffusionCoeff = diffusionCoeff;
            this.ConvectionCoeff = convectionCoeff;
            this.AbsorptionRate = absorptionRate;
        }

        public double DiffusionCoeff { get; }
        public double[] ConvectionCoeff { get; }
        public double AbsorptionRate { get; }

        public ConvectionDiffusionMaterial Clone() => new ConvectionDiffusionMaterial(DiffusionCoeff, ConvectionCoeff, AbsorptionRate);
    }
}
