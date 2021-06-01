namespace MGroup.Optimization.Initialization.HaltonPoints
{
    // Thread safe
    class HaltonSequence
    {
        private readonly int prime;

        public HaltonSequence(int prime)                                                           
        {
            this.prime = prime;
        }

        public double ElementAt(int index)
        {
            double result = 0.0;
            double f = 1.0;
            while (index > 0)
            {
                f /= prime;
                result += f * (index % prime);
                index = index / prime;
            }
            return result;
        }
    }
}
