﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public class MassAccelerationLoad
    {
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}
