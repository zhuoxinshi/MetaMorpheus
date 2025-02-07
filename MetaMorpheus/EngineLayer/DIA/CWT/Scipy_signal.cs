using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Python.Runtime;

namespace EngineLayer.DIA.CWT
{
    public class Scipy_signal
    {

        public static double[] SavgolFilter(double[] data, int windowSize, int polyOrder)
        {
            // Ensure Python runtime is initialized
            if (!PythonEngine.IsInitialized)
            {
                Runtime.PythonDLL = @"C:\Users\Zhuoxin Shi\AppData\Local\Programs\Python\Python313\python313.dll";
                PythonEngine.PythonHome = @"C:\Users\Zhuoxin Shi\AppData\Local\Programs\Python\Python313"; // Adjust your Python path
                PythonEngine.Initialize();
            }

            using (Py.GIL()) // Acquire the Python Global Interpreter Lock
            {
                // Import scipy.signal
                dynamic scipySignal = Py.Import("scipy.signal");

                // Convert C# instance data to Python lists
                dynamic pySignalData = PyObject.FromManagedObject(data);

                // Call savgol_filter
                dynamic smoothedData = scipySignal.savgol_filter(pySignalData, windowSize, polyOrder);

                double[] result = smoothedData.AsManagedObject(typeof(double[])) as double[];

                return result;
            }
        }

        public static List<int> FindPeaks_cwt(double[] data, double[] widths)
        {
            // Ensure Python runtime is initialized
            if (!PythonEngine.IsInitialized)
            {
                Runtime.PythonDLL = @"C:\Users\Zhuoxin Shi\AppData\Local\Programs\Python\Python313\python313.dll";
                PythonEngine.PythonHome = @"C:\Users\Zhuoxin Shi\AppData\Local\Programs\Python\Python313"; // Adjust your Python path
                PythonEngine.Initialize();
            }

            using (Py.GIL()) // Acquire the Python Global Interpreter Lock
            {
                // Import scipy.signal
                dynamic scipySignal = Py.Import("scipy.signal");

                // Convert C# instance data to Python lists
                dynamic pySignalData = PyObject.FromManagedObject(data);
                dynamic pyWidths = PyObject.FromManagedObject(widths);

                // Call find_peaks_cwt
                dynamic peaks = scipySignal.find_peaks_cwt(pySignalData, pyWidths);

                // Convert Python result to C# list
                var result = new List<int>();
                foreach (var peak in peaks)
                {
                    result.Add((int)peak);
                }

                return result;
            }
        }

       
    }
}
