﻿using ISAAR.MSolve.LinearAlgebra.Matrices;


namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IContinuumMaterial3DDefGrad : IContinuumMaterial3D//, IFiniteElementMaterial
    {
        //double[] Stresses { get; }
        //IMatrixView ConstitutiveMatrix { get; }
        //void UpdateMaterial(double[] strains);
        void ClearState();
        void SaveState();
        void ClearStresses();
    }
}
