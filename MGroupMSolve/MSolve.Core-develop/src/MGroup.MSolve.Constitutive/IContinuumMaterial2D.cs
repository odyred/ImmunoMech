using MGroup.LinearAlgebra.Matrices;

namespace MGroup.MSolve.Constitutive
{
	/// <summary>
	/// Interface for materials to be used in 2D finite elements
	/// </summary>
	public interface IContinuumMaterial2D : IFiniteElementMaterial
	{
		double[] Stresses { get; }
		IMatrixView ConstitutiveMatrix { get; }
		void UpdateMaterial(double[] strains);
	}
}
