import java.util.Iterator;

// collection of functions involved in base64 encoding and decoding
public class Base64Parser
{
	public static final String base64characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	
	public Base64Parser()
	{}

	// decode base64 into array of 32-bit integers
	public static int[] base64toInteger(String base64)
	{
		int[] result = new int[base64.length() * 6 / 32]; // ignore the tail
		int current = 0;
		String binaryString = "";

		for (int i = 0; i < base64.length() && current < result.length; i++)
		{
			binaryString += decodeBase64(base64.charAt(i));
			// System.out.println(binaryString.length());

			if (binaryString.length() >= 32)
			{
				result[current] = (int) Long.parseLong(binaryString.substring(0, 32), 2);
				binaryString = binaryString.substring(32, binaryString.length());
				current++;
			}
		}

		return result;
	}

	// convert 32-bit or 64-bit binary to floating point
	public static double binarytoDouble(String binary)
	{
		if (binary.length() == 32) // single-precision
		{
			double sign = binary.charAt(0) - 48.0;
			double exponent = 0.0;
			double factor = 1.0;

			for (int i = 0; i < 8; i++)
			{
				exponent += factor * (binary.charAt(8 - i) - 48.0);
				factor *= 2;
			}

			double significand = 0;
			factor = 0.5;

			for (int i = 0; i < 23; i++)
			{
				significand += factor * (binary.charAt(9 + i) - 48.0);
				factor /= 2;
			}

			return (1 - (2 * sign)) * (1 + significand) * (Math.pow(2.0, exponent - 127.0));
		}

		if (binary.length() == 64) // double-precision
		{
			double sign = binary.charAt(0) - 48.0;
			double exponent = 0.0;
			double factor = 1.0;

			for (int i = 0; i < 11; i++)
			{
				exponent += factor * (binary.charAt(11 - i) - 48.0);
				factor *= 2;
			}

			double significand = 0;
			factor = 0.5;

			for (int i = 0; i < 52; i++)
			{
				significand += factor * (binary.charAt(12 + i) - 48.0);
				factor /= 2;
			}

			return (1 - (2 * sign)) * (1 + significand) * (Math.pow(2.0, exponent - 1023.0));
		}

		HelperFunctions.debug("Unexpected binary string length");
		return -1;
	}

	// decode base64 into 2D array of double, e.g. spectrum
	// work for either 32-bit (float) or 64-bit (double)
	public static double[][] base64toDouble(String base64, int precision)
	{
		// System.out.println(base64.length());
		// System.out.println(base64.length() * 3 / precision);
		double[][] result = new double[base64.length() * 3 / precision][2]; // m/z and intensity
		int current = 0;
		String binaryString = "";

		if (precision == 32) // 32-bit float
		{
			for (int i = 0; i < base64.length() && current < result.length; i++)
			{
				binaryString += decodeBase64(base64.charAt(i));

				if (binaryString.length() >= 64)
				{
					// System.out.println(binaryString);
					result[current][0] = binarytoDouble(binaryString.substring(0, 32));
					result[current][1] = binarytoDouble(binaryString.substring(32, 64));
					current++;
					binaryString = binaryString.substring(64, binaryString.length());
				}
			}
		}

		else // 64-bit double
		{
			for (int i = 0; i < base64.length() && current < result.length; i++)
			{
				binaryString += decodeBase64(base64.charAt(i));

				if (binaryString.length() >= 128)
				{
					// System.out.println(binaryString);
					result[current][0] = binarytoDouble(binaryString.substring(0, 64));
					result[current][1] = binarytoDouble(binaryString.substring(64, 128));
					current++;
					binaryString = binaryString.substring(128, binaryString.length());
				}
			}
		}

		return result;
	}

	// convert CrosslinkStruct to base64
	public static String encodeBase64(CrosslinkStruct crosslink, int maxVarModPerPeptide)
	{ return encodeBase64(toBinaryString(crosslink.peptideA, maxVarModPerPeptide) + toBinaryString(crosslink.peptideB, maxVarModPerPeptide)); }

	// convert PeptideStruct to binary string
	public static String toBinaryString(PeptideStruct peptide, int maxVarModPerPeptide)
	{
		if (peptide == null) // for loop crosslink
		{}

		String binary = "";
		String padding = "00000000000000000000000000000000"; // padding
		String tempbinary = Integer.toBinaryString(peptide.parent.entryID);

		binary += padding.substring(0, 32 - tempbinary.length()) + tempbinary;
		tempbinary = Integer.toBinaryString(peptide.startPos);
		binary += padding.substring(0, 32 - tempbinary.length()) + tempbinary;
		tempbinary = Integer.toBinaryString(peptide.endPos);
		binary += padding.substring(0, 32 - tempbinary.length()) + tempbinary;

		tempbinary = Integer.toBinaryString(-1);
		String varModPadding = padding.substring(0, 32 - tempbinary.length()) + tempbinary;
		Integer tempint;
		int numMod = peptide.varModMap.size();

		if (numMod == 0) // no modification
		{
			for (int i = 0; i < maxVarModPerPeptide; i++) // make every peptide have 'maxVarModPerPeptide' pairs
				binary += varModPadding + varModPadding;
		}

		else // some modification
		{
			for (Iterator iter = peptide.varModMap.keySet().iterator(); iter.hasNext();)
			{
				tempint = (Integer) iter.next();
				tempbinary = Integer.toBinaryString(tempint.intValue());
				binary += padding.substring(0, 32 - tempbinary.length()) + tempbinary;
				tempbinary = Integer.toBinaryString(peptide.varModMap.get(tempint).entryID);
				binary += padding.substring(0, 32 - tempbinary.length()) + tempbinary;
			}

			for (int i = numMod; i < maxVarModPerPeptide; i++) // make every peptide have 'maxVarModPerPeptide' pairs
				binary += varModPadding + varModPadding;
		}

		binary += padding.substring(0, (maxVarModPerPeptide % 3) * 2);
		return binary;
	}

	// convert PeptideStruct to base64
	public static String encodeBase64(PeptideStruct peptide, int maxVarModPerPeptide)
	{ return encodeBase64(toBinaryString(peptide, maxVarModPerPeptide)); }

	// convert binary string to base64
	public static String encodeBase64(String binaryString)
	{
		if (binaryString.length() % 6 > 0)
		{
			HelperFunctions.debug("Binary string length not divisible by 6!");
			return null;
		}

		else
		{
			String encoded = "";

			for (int i = 0; i < binaryString.length() / 6; i++)
				encoded += "" + base64characters.charAt(Integer.parseInt(binaryString.substring(6 * i, 6 * (i + 1)), 2));

			return encoded;
		}
	}

	// convert base64 string to binary
	public static String decodeBase64(String base64)
	{
		String binaryString = "", padding = "000000";

		for (int i = 0; i < base64.length(); i++)
		{
			String tempbinary = Integer.toString(base64characters.indexOf(base64.charAt(i)), 2);
			binaryString += padding.substring(0, 6 - tempbinary.length()) + tempbinary; // make sure the results are 8-bit
		}

		return binaryString;
	}

	// convert base64 character to binary
	public static String decodeBase64(char base64)
	{
		String padding = "000000";
		String binaryString = Integer.toString(base64characters.indexOf(base64), 2);
		binaryString = padding.substring(0, 6 - binaryString.length()) + binaryString; // make sure the results are 6-bit

		return binaryString;
	}
}