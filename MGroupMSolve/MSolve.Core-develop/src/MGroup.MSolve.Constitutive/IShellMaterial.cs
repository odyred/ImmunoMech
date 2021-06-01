using MGroup.LinearAlgebra.Matrices;

namespace MGroup.MSolve.Constitutive
{
	/// <summary>
	/// Interface for materials laws implementations to be used with shell finite elements
	/// </summary>
	public interface IShellMaterial : IFiniteElementMaterial
	{
		double[] Stresses { get; }
		IMatrixView ConstitutiveMatrix { get; }
		new IShellMaterial Clone();
		double[] NormalVectorV3 { get; set; }
		double[] TangentVectorV2 { get; set; }
		double[] TangentVectorV1 { get; set; }
		void UpdateMaterial(double[] greenLagrangeStrains);
	}
}
