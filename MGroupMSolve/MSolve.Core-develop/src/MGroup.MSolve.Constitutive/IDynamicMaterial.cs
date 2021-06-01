namespace MGroup.MSolve.Constitutive
{
	/// <summary>
	/// Contains material properties used for dynamic or modal analysis. These are uniform for the whole element and immutable.
	/// Authors: Serafeim Bakalakos
	/// </summary>
	public interface IDynamicMaterial
	{
		double Density { get; }
		double RayleighCoeffMass { get; }
		double RayleighCoeffStiffness { get; }
	}
}
