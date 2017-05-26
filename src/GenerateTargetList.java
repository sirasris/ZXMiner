import java.text.DecimalFormat;
import java.util.*;
import java.io.*;

// create target list for mass spectrometer
public class GenerateTargetList
{
	public GenerateTargetList()
	{}
	
	public static void generate(ParamStruct param)
	{
		TreeMap<String, TreeSet<Double>> globalMap = new TreeMap<String, TreeSet<Double>>();
		TreeMap<String, Double> globalMassError = new TreeMap<String, Double>();
		
		try
		{
			// check if file exists
			File target = new File(".\\temp\\" + param.outputFileName + "_targetList.temp");
			
			if (!target.exists()) // terminate with error message
			{
				System.out.println("Error! Temporary target list file <\\temp\\" + param.outputFileName + "_targetList.temp> does not exist!");
				System.out.println("No target list is generated!");
				return;
			}
			
			BufferedReader infile = new BufferedReader(new FileReader(".\\temp\\" + param.outputFileName + "_targetList.temp"));
			String st = infile.readLine();
			String[] tempst;
			String precursor;
			Double RT;
			TreeSet<Double> tempTreeSet;
			
			while (st != null) // for each entry
			{
				tempst = st.split("\t");
				precursor = (new DecimalFormat("#.###")).format(Double.valueOf(tempst[0]));
				RT = new Double(Double.valueOf(tempst[1]));
				
				if (!globalMap.containsKey(precursor)) // new precursor
				{
					tempTreeSet = new TreeSet<Double>();
					tempTreeSet.add(RT);
					globalMap.put(precursor, tempTreeSet);
					globalMassError.put(precursor, Double.valueOf(tempst[2]));
				}
				
				else // existing precursor
				{
					tempTreeSet = globalMap.get(precursor);
					tempTreeSet.add(RT); // TreeSet handles sorting and removing duplicates
					
					if (Double.valueOf(tempst[2]) < globalMassError.get(precursor).doubleValue()) // update mass error
						globalMassError.put(precursor, Double.valueOf(tempst[2]));
				}
				
				st = infile.readLine();
			}
			
			infile.close();
			
			BufferedWriter outfile = new BufferedWriter(new FileWriter(".\\output\\" + param.outputFileName + "_targetList.csv"));
			ArrayList<Double> tempResult;
			
			for (Iterator<String> iter = globalMap.keySet().iterator(); iter.hasNext();) // for each precursor
			{
				precursor = iter.next();
				tempResult = mergeRTWindows(globalMap.get(precursor), param.targetListRTHalfWindow);
				
				for (int i = 0; i < tempResult.size(); i += 2)
					outfile.write(precursor + "," + tempResult.get(i) + "," + tempResult.get(i + 1) + ",,35\n");
				
				// System.out.println(globalMassError.get(precursor).toString());
			}
			
			outfile.flush();
			outfile.close();
		}
		
		catch (Exception e)
		{ HelperFunctions.debug("GenerateTargetList", HelperFunctions.getStackTrace(e)); }
	}
	
	// process a sorted set of peak RT by merging overlapping windows
	public static ArrayList<Double> mergeRTWindows(TreeSet<Double> input, double halfWindow)
	{
		ArrayList<Double> result = new ArrayList<Double>();
		Double current = input.pollFirst(); // remove lowest RT
		result.add(new Double(current.doubleValue() - halfWindow));
		
		mergeRTWindowsHelper(result, input, halfWindow, current.doubleValue() + halfWindow);	
		
		return result;
	}
	
	public static void mergeRTWindowsHelper(ArrayList<Double> result, TreeSet<Double> input, double halfWindow, double currentUpperBound)
	{
		if (input.size() == 0) // end of list, add the upper bound and be done
		{
			result.add(new Double(currentUpperBound));
			return;
		}
		
		else
		{
			Double current = input.pollFirst(); // remove next RT
			
			if (current.doubleValue() - halfWindow > currentUpperBound) // next window is disjoint
			{
				result.add(new Double(currentUpperBound)); // close the current window
				result.add(new Double(current.doubleValue() - halfWindow)); // start a new one
				mergeRTWindowsHelper(result, input, halfWindow, current.doubleValue() + halfWindow);
			}
			
			else // next window overlap the current one
				mergeRTWindowsHelper(result, input, halfWindow, current.doubleValue() + halfWindow);
		}
	}
}