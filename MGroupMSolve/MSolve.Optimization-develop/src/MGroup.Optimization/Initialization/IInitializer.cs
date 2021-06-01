namespace MGroup.Optimization.Initialization
{
    public interface IInitializer<T>
    {
        T[] Generate();
    }
}
