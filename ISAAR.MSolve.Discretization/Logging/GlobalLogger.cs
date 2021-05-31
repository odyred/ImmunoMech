using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;

namespace ISAAR.MSolve.Discretization.Logging
{
    public static class GlobalLogger
    {
        private static StreamWriter writer;

        public static void OpenOutputFile(string path)
        {
            try
            {
                writer = new StreamWriter(path);
                writer.AutoFlush = true;
                Debug.WriteLine("GlobalLogger opened requested file successfully.");

            }
            catch (Exception)
            {
                Debug.WriteLine("GlobalLogger cannot open requested file.");
                writer = null;
            }
        }

        public static void CloseCurrentOutputFile()
        {
            if (writer != null)
            {
                writer.Flush();
                writer.Dispose();
                writer = null;
                Debug.WriteLine("GlobalLogger closed output file.");
            }
        }

        public static void WriteLine(string msg)
        {
            writer.WriteLine(msg);
        }
    }
}
