using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials
{
    public class ConvectionDiffusionMaterial
    {
        public ConvectionDiffusionMaterial(double capacityCoeff, double diffusionCoeff, double[] convectionCoeff, double absorptionRate)
        {
            this.CapacityCoeff = capacityCoeff; 
            this.DiffusionCoeff = diffusionCoeff;
            this.ConvectionCoeff = convectionCoeff;
            this.AbsorptionRate = absorptionRate;
        }

        public double CapacityCoeff { get; }
        public double DiffusionCoeff { get; }
        public double[] ConvectionCoeff { get; }
        public double AbsorptionRate { get; }

        public ConvectionDiffusionMaterial Clone() => new ConvectionDiffusionMaterial(CapacityCoeff, DiffusionCoeff, ConvectionCoeff, AbsorptionRate);
    }
}
