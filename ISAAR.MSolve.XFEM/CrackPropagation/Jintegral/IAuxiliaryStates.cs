﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.CrackGeometry.CrackTip;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Jintegral
{
    interface IAuxiliaryStates
    {
        AuxiliaryStatesTensors ComputeTensorsAt(CartesianPoint2D globalIntegrationPoint, 
            TipCoordinateSystem tipCoordinateSystem);
    }
}
