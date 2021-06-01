using MGroup.LinearAlgebra.Vectors;

namespace MGroup.MSolve.AnalysisWorkflow
{
	public interface INonLinearSubdomainUpdater
	{
		void ScaleConstraints(double scalingFactor);
		IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution);
		void UpdateState();
		void ResetState();
	}
}
