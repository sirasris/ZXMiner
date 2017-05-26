import java.util.Arrays;
import java.util.ArrayList;
import java.io.*;

// contain misc. helper functions for other classes
public final class HelperFunctions
{
	public static final String dummyPeptideReport = "n/a\tn/a\tn/a\tn/a\tn/a\tn/a\tn/a";
	public static final String[] tempScoreFileHeader = {"RawFile", "ScanNumber", "ChargeState", "PrecursorM/Z", "PrecursorM+H",
		"ParentScanNum", "RetentionTime", "BasePeak", "PrecursorInt", "MS/MSPeakCount", 
		"ProteinName-A", "ProteinDescription-A", "ProteinIsForward-A", "StartResidue-A", "EndResidue-A", "Sequence-A", "Modification-A", "XlinkSite-A",
		"ProteinName-B", "ProteinDescription-B", "ProteinIsForward-B", "StartResidue-B", "EndResidue-B", "Sequence-B", "Modification-B", "XlinkSite-B",
		"PrecursorPPMError", "PeakCoverage", "IntCoverage", "IonCoverage-A-b", "IonCoverage-A-y", "IonCoverage-A", 
		"IonCoverage-B-b", "IonCoverage-B-y", "IonCoverage-B", "IonCoverage", "GMScore"};
	
	private HelperFunctions()
	{}

	// output to file
	// synchronized
	public static synchronized void appendToFile(String fpath, String[] lines)
	{
		try
		{
			BufferedWriter outfile = new BufferedWriter(new FileWriter(fpath, true));

			for (int i = 0; i < lines.length; i++)
				outfile.write(lines[i] + '\n');

			outfile.flush();
			outfile.close();
		}

		catch (Exception e)
		{ debug("HelperFunctions::appendToFile", HelperFunctions.getStackTrace(e)); }
	}
	
	// output to file
	// synchronized
	public static synchronized void appendToFile(String fpath, String line)
	{
		try
		{
			BufferedWriter outfile = new BufferedWriter(new FileWriter(fpath, true));
			outfile.write(line + "\n");
			
			outfile.flush();
			outfile.close();
		}

		catch (Exception e)
		{ debug("HelperFunctions::appendToFile", HelperFunctions.getStackTrace(e)); }
	}

	// report progress and update elapsed time
	public static long reportStatus(String message, long startTime)
	{
		long currentTime = System.currentTimeMillis();
		String padding = "                                                         ";
		System.out.println("Status: " + message + padding.substring(0, padding.length() - message.length()) + "Elapsed Time: " + ((currentTime - startTime) / 1000.0) + " seconds.");
		
		return currentTime;
	}

	public static String arrayToString(double[] data)
	{
		String result = data[0] + "";
		
		for (int i = 1; i < data.length; i++)
			result += ", " + data[i];
		
		return result;
	}
	
	// compute mass errors
	public static double getAbsppmError(double ref, double target)
	{ return Math.abs(getppmError(ref, target)); }

	public static double getppmError(double ref, double target)
	{ return (target - ref) * 1000000 / ref; }

	public static double getAbsError(double ref, double target)
	{ return Math.abs(getError(ref, target)); }

	public static double getError(double ref, double target)
	{ return target - ref; }

	// get possible charge state range for a fragmented ion
	public static int[] getChargeRange(FragmentIonStruct ion, int precursorChargeState)
	{
		int[] range = new int[2];
		
		if (ion.ionType == 'p') // precursor
		{
			range[0] = precursorChargeState;
			range[1] = precursorChargeState;
			
			return range;
		}

		if (ion.isCrosslink) // crosslinked ion
		{
			range[0] = 2;
			range[1] = precursorChargeState;
		}

		else // linear ion
		{
			range[0] = 1;
			range[1] = precursorChargeState - 1;

			if (ion.ionID < 6)
				range[1] = 1;
			if (ion.ionID > 12)
				range[0] = 2;
		}

		return range;
	}

	// get mass array out of ArrayList<ChargedPeakStruct>at specified charge state
	// take reference to output variable
	public static void getMassArray(ArrayList<ChargedPeakStruct> sortedInput, int chargeState, ArrayList<Integer[]> outputID, ArrayList<Double[]> output)
	{
		ArrayList<Integer> tempID = new ArrayList<Integer>();
		ArrayList<Double> temp = new ArrayList<Double>();

		for (int i = 0; i < sortedInput.size(); i++)
		{
			if (sortedInput.get(i).chargeState == chargeState)
			{
				tempID.add(new Integer(i));
				temp.add(new Double(sortedInput.get(i).massWithOneCharge));
			}
		}

		outputID.add(tempID.toArray(new Integer[0]));
		output.add(temp.toArray(new Double[0]));
	}

	// get mass array out of ArrayList<ChargedPeakStruct>
	public static Double[] getMassArray(ArrayList<ChargedPeakStruct> sortedInput)
	{
		Double[] sortedOutput = new Double[sortedInput.size()];

		for (int i = 0; i < sortedInput.size(); i++)
			sortedOutput[i] = sortedInput.get(i).massWithOneCharge;

		return sortedOutput;
	}

	// sort and return indexed array
	public static Integer[] getIndexArray(Object[] input)
	{
		IndexedComparator idcomp = new IndexedComparator(input);
		Integer[] index = idcomp.initIndexArray();
		Arrays.sort(index, idcomp);
		Arrays.sort(input);

		return index;
	}

	// print out error message for debugging purposes
	public static void debug(String classname, String message)
	{ System.out.println("Error in \"" + classname + "\": " + message); }

	// print out variable value for debugging purposes
	public static void debug(String varname, Object var)
	{
		if (var instanceof Object[])
		{
			Object[] tempcast = (Object[]) var;
			System.out.print(varname + " values: ");

			for (int i = 0; i < tempcast.length; i++)
				System.out.print(tempcast[i].toString() + ", ");

			System.out.println();
			return;
		}
		
		if (var instanceof int[])
		{
			int[] tempcast = (int[]) var;
			System.out.print(varname + " values: ");

			for (int i = 0; i < tempcast.length; i++)
				System.out.print(tempcast[i]+ ", ");

			System.out.println();
			return;
		}
		
		if (var instanceof double[])
		{
			double[] tempcast = (double[]) var;
			System.out.print(varname + " values: ");

			for (int i = 0; i < tempcast.length; i++)
				System.out.print(tempcast[i]+ ", ");

			System.out.println();
			return;
		}
		
		if (var instanceof boolean[])
		{
			boolean[] tempcast = (boolean[]) var;
			System.out.print(varname + " values: ");

			for (int i = 0; i < tempcast.length; i++)
				System.out.print(tempcast[i]+ ", ");

			System.out.println();
			return;
		}
		
		System.out.println(varname + " value: " + var.toString());
	}

	// print out message
	public static void debug(String message)
	{ System.out.println(message); }

	// place-holder function
	public static void toBeImplemented(String note)
	{ System.out.println("To be implemente: " + note); }
	
	public static String getStackTrace(Exception e)
	{
		StringWriter error = new StringWriter();
		e.printStackTrace(new PrintWriter(error));
		
		return error.toString();
	}
}