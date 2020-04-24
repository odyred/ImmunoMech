using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Materials
{
    public class ConvectionDiffusionMaterial
    {
        public ConvectionDiffusionMaterial(double diffusionCoeff, double convectionCoeff)
        {
            this.DiffusionCoeff = diffusionCoeff;
            this.ConvectionCoeff = convectionCoeff;
        }

        public double DiffusionCoeff { get; }
        public double ConvectionCoeff { get; }

        public ConvectionDiffusionMaterial Clone() => new ConvectionDiffusionMaterial(DiffusionCoeff, ConvectionCoeff);
    }
}
