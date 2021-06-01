﻿using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISSAR.MSolve.Discretization.Loads;

//TODO: This is probably covered by Load.cs
namespace ISAAR.MSolve.FEM.Entities
{
    /// <summary>
    /// You should use <see cref="Load"/> instead.
    /// </summary>
    public class SteadyNodalLoad: ITimeDependentNodalLoad
    {
        private readonly double constantloadAmount;

        public SteadyNodalLoad(double constantloadAmount)
        {
            this.constantloadAmount = constantloadAmount;
        }

        public INode Node { get; set; }
        public IDofType DOF { get; set; }

        public double GetLoadAmount(int timeStep)
        {
            return constantloadAmount;
        }
    }
}
