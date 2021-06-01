using System;

namespace MGroup.MSolve.Constitutive
{
	/// <summary>
	/// Interface for materials laws implementations to be used in finite elements 
	/// </summary>
	public interface IFiniteElementMaterial: ICloneable
	{
		//this is only for Structural elements with strain and stress tensors.
		int ID { get; }
		bool Modified { get; }
		void ResetModified();
		double[] Coordinates { get; set; }
		double YoungModulus { get; }
		double PoissonRatio { get; }

		void SaveState();
		void ClearState();
		void ClearStresses();
	}
}
