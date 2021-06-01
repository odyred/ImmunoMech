using System.Reflection;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using Xunit;

// General Information about an assembly is controlled through the following 
// set of attributes. Change these attribute values to modify the information
// associated with an assembly.

// Setting ComVisible to false makes the types in this assembly not visible 
// to COM components.  If you need to access a type in this assembly from 
// COM, set the ComVisible attribute to true on that type.
[assembly: ComVisible(false)]

// The following GUID is for the ID of the typelib if this project is exposed to COM
[assembly: Guid("ad892d2a-7436-4e18-8375-c9c8dd4d0eae")]

// Version information for an assembly consists of the following four values:
//
//      Major Version
//      Minor Version 
//      Build Number
//      Revision
//

// The following will prevent XUnit from running tests in parallel. Normally it should not be necessary, but up until 10/12/2018
// it was not worth the effort to make the classes used in some tests thread-safe.
[assembly: CollectionBehavior(DisableTestParallelization = true)]