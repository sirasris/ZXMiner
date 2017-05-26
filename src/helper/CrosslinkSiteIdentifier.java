import java.util.ArrayList;
 
// collection of functions involved in determining crosslink site
public class CrosslinkSiteIdentifier
{
	public CrosslinkSiteIdentifier()
	{}
	
	// check if two peptides are crosslinkable
	public static boolean isCrosslinkable(PeptideStruct peptide1, PeptideStruct peptide2, CrosslinkerStruct crosslinker)
	{
		boolean[] site1 = containsCrosslinkSites(peptide1, crosslinker);
		boolean[] site2 = containsCrosslinkSites(peptide2, crosslinker);

		return (site1[0] && site2[1]) || (site1[1] && site2[0]);
	}

	// check if two peptides are crosslinkable, first peptide is preprocessed
	public static boolean isCrosslinkable(boolean[] site1, PeptideStruct peptide2, CrosslinkerStruct crosslinker)
	{
		boolean[] site2 = containsCrosslinkSites(peptide2, crosslinker);

		return (site1[0] && site2[1]) || (site1[1] && site2[0]);
	}

	// check if a peptide can form loop
	public static boolean isSelfCrosslinkable(PeptideStruct peptide, CrosslinkerStruct crosslinker)
	{
		return false; // unimplemented for now
	}

	// check if a certain peptide contain crosslink sites
	public static boolean[] containsCrosslinkSites(PeptideStruct peptide, CrosslinkerStruct crosslinker)
	{
		boolean[] result = {false, false};
		String sequence = peptide.getSequence(true);

		for (int i = 0; i < sequence.length() - 1; i++) // generally not consider cleaved residue, this will include the C-terminal residue of the protein
		{
			if (!peptide.isModified(i)) // not consider modified residue
			{
				if (result[0] && result[1])
					break;
				
				if (crosslinker.siteA.contains(sequence.charAt(i) + ""))
					result[0] = true;
				if (crosslinker.siteB.contains(sequence.charAt(i) + ""))
					result[1] = true;
			}
		}
		
		if ((!result[0] || !result[1]) && peptide.getLastResidue() == ']') // continue if C-terminal is present and it's necessary
		{
			if (crosslinker.siteA.contains("]"))
				result[0] = true;
			if (crosslinker.siteB.contains("]"))
				result[1] = true;
		}

  		return result;
	}
	
	// list all possible crosslink site on a peptide
	// adjust residue position to coincide with b- and y- ion convention, without teminus symbols
	// crosslinkRule = 'A' or 'B'
	public static int[] crosslinkSites(PeptideStruct peptide, char crosslinkRule, CrosslinkerStruct crosslinker)
	{
		String sequence = peptide.getSequence(false); // retrieve sequence without terminus symbol
		String rule = crosslinker.getCrosslinkRule(crosslinkRule);
		ArrayList<Integer> positions = new ArrayList<Integer>();
		
		for (int i = 0; i < sequence.length() - 1; i++) // generally not consider cleaved residue at the C-terminal end
		{
			if (!peptide.isModified(i) && rule.contains(sequence.charAt(i) + "")) // unmodified and fit the rule
				positions.add(new Integer(i + 1));
		}
		
		if (peptide.getLastResidue() == ']') // contain C-terminal, need to consider last residue + terminal
		{
			if ((!peptide.isModified(sequence.length() - 1) && rule.contains(sequence.charAt(sequence.length() - 1) + "")) 
					|| rule.contains("]")) // last residue is unmodified and fit the rule or the rule involve C-terminal
				positions.add(new Integer(sequence.length()));
		}
		
		if (peptide.getFirstResidue() == '[' && rule.contains("[") && !positions.contains(new Integer(1))) // N-terminal
		{
			if (positions.size() > 0)
				positions.add(0, new Integer(1)); // add to the front
			else
				positions.add(new Integer(1));
		}
			
		int[] output = new int[positions.size()];
		
		for (int i = 0; i < positions.size(); i++)
			output[i] = positions.get(i);
		
		// if (peptide.getSequence(true).equals("[GSMEQFPKETVVESSGPK") && crosslinkRule == 'A')
		// 	HelperFunctions.debug("crosslink sites", output);
		
		return output;
	}
	
	// return bounary of crosslink sites
	// for computing overall scores
	public static int[][] crosslinkSiteBound(CrosslinkStruct crosslink, CrosslinkerStruct crosslinker)
	{
		int[] siteA_A = crosslinkSites(crosslink.peptideA, 'A', crosslinker);
		int[] siteA_B = crosslinkSites(crosslink.peptideA, 'B', crosslinker);
		int[] siteB_A = crosslinkSites(crosslink.peptideB, 'A', crosslinker);
		int[] siteB_B = crosslinkSites(crosslink.peptideB, 'B', crosslinker);
		
		int[][] boundary = {{crosslink.peptideA.getSequenceLength(false), 1}, {crosslink.peptideB.getSequenceLength(false), 1}}; // minA, maxA, minB, maxB
		
		// HelperFunctions.debug("crosslink sites bound", boundary);
		
		if (siteA_A.length > 0 && siteB_B.length > 0)
		{
			boundary[0][0] = siteA_A[0];
			boundary[0][1] = siteA_A[siteA_A.length - 1];
			boundary[1][0] = siteB_B[0];
			boundary[1][1] = siteB_B[siteB_B.length - 1]; // sites are pre-sorted
		}
		
		if (siteA_B.length > 0 && siteB_A.length > 0) // update
		{
			boundary[0][0] = Math.min(boundary[0][0], siteA_B[0]);
			boundary[0][1] = Math.max(boundary[0][1], siteA_B[siteA_B.length - 1]);
			boundary[1][0] = Math.min(boundary[1][0], siteB_A[0]);
			boundary[1][1] = Math.max(boundary[1][1], siteB_A[siteB_A.length - 1]); // sites are pre-sorted
		}
		
		// HelperFunctions.debug("crosslink sites bound", boundary[0]);
		// HelperFunctions.debug("crosslink sites bound", boundary[1]);
		
		return boundary;
	}
}