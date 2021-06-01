using MGroup.Optimization.Commons;
using MGroup.Optimization.Problems;

namespace MGroup.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings
{
    public class StandardBinaryCoding : AbstractBinaryCoding
    {
        public StandardBinaryCoding(OptimizationProblem problem, int bitsPerContinuousVariable, int bitsPerIntegerVariable) :
                        base(problem, bitsPerContinuousVariable, bitsPerIntegerVariable)
        {
        }

        protected sealed override int BitstringToDecimalInteger(bool[] bits, int start, int length)
        {
            return BinaryUtilities.StandardBinaryToDecimal(bits, start, length);
        }

        protected override void DecimalIntegerToBitstring(int dec, bool[] bits, int start, int length)
        {
            BinaryUtilities.DecimalToStandardBinary(dec, bits, start, length);
        }
    }
}
