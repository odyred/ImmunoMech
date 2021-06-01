﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.Commons
{
    internal static class Utilities
    {
        internal static bool AreEqual(int[] array1, int[] array2)
        {
            if (array1.Length != array2.Length) return false;
            for (int i = 0; i < array1.Length; ++i)
            {
                if (array1[i] != array2[i]) return false;
            }
            return true;
        }
    }
}
