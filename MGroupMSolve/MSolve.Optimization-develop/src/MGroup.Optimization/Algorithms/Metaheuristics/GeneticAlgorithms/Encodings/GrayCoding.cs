using MGroup.Optimization.Commons;
using MGroup.Optimization.Problems;

namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings
{
    public class GrayCoding : AbstractBinaryCoding
    {
        public GrayCoding(OptimizationProblem problem, int bitsPerContinuousVariable, int bitsPerIntegerVariable):
                        base(problem, bitsPerContinuousVariable, bitsPerIntegerVariable)
        {
        }

        protected sealed override int BitstringToDecimalInteger(bool[] bits, int start, int length)
        {
            return BinaryUtilities.GrayCodeToDecimal(bits, start, length);
        }

        protected override void DecimalIntegerToBitstring(int dec, bool[] bits, int start, int length)
        {
            BinaryUtilities.DecimalToGrayCode(dec, bits, start, length);
        }
    }
}
