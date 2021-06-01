

//TODO: Add Gray scale filter, as presented in "An efficient 3D topology optimization code written in Matlab (2014)", which
//      requires calculations during the OC optimizer
//TODO: Add filters presented in other publications with Matlab codes.

using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Optimization.Structural.Topology.SIMP.Filtering
{
    public interface IDensityFilter
    {
        void FilterSensitivities(Vector densities, ref Vector sensitivities); //TODO: perhaps use IVectorView?
    }
}
