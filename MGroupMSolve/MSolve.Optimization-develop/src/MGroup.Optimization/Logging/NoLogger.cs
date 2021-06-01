namespace MGroup.Optimization.Logging
{
    /// <summary>
    /// An <see cref="IOptimizationLogger"/> that does nothing.
    /// </summary>
    public class NoLogger : IOptimizationLogger
    {
        public void Log(IOptimizationAlgorithm algorithm)
        {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        public void PrintToConsole()
        {
        }
    }
}
