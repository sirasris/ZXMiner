import java.util.TreeMap;
import java.util.Iterator;

// class of functions involved in generating neutral losses
public class NeutralLossGenerator
{
	public static final String waterLossResidues = "STDE";
	public static final String ammoniaLossResidues = "RKQN";
	public static final String[] neutralLossName = {"H2O", "NH3", "OHSCH3"};
	
    public NeutralLossGenerator()
    {}
   

	// check if two strings share a common character
	public static boolean intersect(String a, String b)
	{
		for (int i = 0; i < a.length(); i++)
		for (int j = 0; j < b.length(); j++)
			if (a.charAt(i) == b.charAt(j))
				return true;

		return false;
	}
    
    // determine whether a peptide can have neutral losses
	public static boolean[] getNeutralLossFlags(PeptideStruct peptide)
	{
		Integer position;
		boolean[] result = new boolean[3]; // water, ammonia, metox

		TreeMap<Integer, ModificationStruct> varModMap = peptide.varModMap; // look for Met oxidation

		for (Iterator iter = varModMap.keySet().iterator(); iter.hasNext();)
		{
			position = (Integer) iter.next();

			if (peptide.getResidue(position.intValue()) == 'M' && varModMap.get(position).name.equals("oxidation"))
			{
				result[2] = true;
				break;
			}
		}

		String sequence = peptide.getSequence(true);
		result[0] = intersect(sequence, waterLossResidues);
		result[1] = intersect(sequence, ammoniaLossResidues);

		return result;
	}

	// determine whether a crosslink can have neutral losses
	public static boolean[] getNeutralLossFlags(CrosslinkStruct crosslink)
	{
		boolean[] result = getNeutralLossFlags(crosslink.peptideA);
		boolean[] result2 = getNeutralLossFlags(crosslink.peptideB);

		for (int i = 0; i < 3; i++) // combine flags with 'or'
			result[i] = result[i] | result2[i];

		return result;
	}

	// determine whether a fragmented ion can have neutral losses
	public static boolean[] getNeutralLossFlags(PeptideStruct peptide, char ionType, int ionID)
	{
		int startPos = -1, endPos = -1;
		boolean[] result = new boolean[3]; // water, ammonia, metox

		switch (ionType)
		{
			case 'b':
			{
				startPos = 0;
				endPos = ionID - 1;
				break;
			}

			case 'y':
			{
				startPos = peptide.getSequenceLength(false) - ionID;
				endPos = peptide.getSequenceLength(false) - 1;
				break;
			}

			default:
			{
				HelperFunctions.debug("Unexpected fragmented ion types");
				break;
			}
		}

		String sequence = peptide.getSequence(startPos, endPos, false);
		result[0] = intersect(sequence, waterLossResidues);
		result[1] = intersect(sequence, ammoniaLossResidues);
		result[2] = false; // only consider MetOx neutral loss from precursor

		return result;
	}

	// determine whether a fragmented ion can have neutral losses
	public static boolean[] getNeutralLossFlags(CrosslinkStruct crosslink, char ionType, int ionID, boolean isCrosslink, boolean fromSecondPeptide)
	{
		boolean[] result;

		if (fromSecondPeptide) // obtain flags from fragmented peptide
			result = getNeutralLossFlags(crosslink.peptideB, ionType, ionID);
		else
			result = getNeutralLossFlags(crosslink.peptideA, ionType, ionID);

		if (isCrosslink)
		{
			boolean[] result2;

			if (fromSecondPeptide) // obtain flags from full peptide part
				result2 = getNeutralLossFlags(crosslink.peptideA);
			else
				result2 = getNeutralLossFlags(crosslink.peptideB);

			for (int i = 0; i < 3; i++) // combine flags with 'or'
				result[i] = result[i] | result2[i];
		}

		return result;
	}
	
	// determine whether a residue fit any neutral loss rule, then update previous flag
	public static void updateNeutralLossFlags(PeptideStruct peptide, int position, boolean[] previousFlag)
	{
		if (!peptide.isModified(position))
		{
			String sequence = peptide.getResidue(position) + "";
			
			if (!previousFlag[0] && intersect(sequence, waterLossResidues))
				previousFlag[0] = true;
			if (!previousFlag[1] && intersect(sequence, ammoniaLossResidues))
				previousFlag[1] = true;
			// not consider met oxidation for fragment ions
		}
	}
}