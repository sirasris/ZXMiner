import java.util.Stack;

// generate theoretical fragmented ions
public class FragmentGenerator
{
	public FragmentGenerator()
	{}

	public static Stack<FragmentIonStruct> fragment(PeptideStruct peptide)
	{
		Stack<FragmentIonStruct> allIons = new Stack<FragmentIonStruct>();
		allIons.push(new FragmentIonStruct(peptide));

		for (int i = 1; i < peptide.getSequenceLength(false); i++) // add b and y
		{
			allIons.push(new FragmentIonStruct('b', i, peptide));
			allIons.push(new FragmentIonStruct('y', i, peptide));
		}

		return allIons;
	}

	public static Stack<FragmentIonStruct> fragment(CrosslinkStruct crosslink, ParamStruct param)
	{
		Stack<FragmentIonStruct> allIons = new Stack<FragmentIonStruct>();
		allIons.push(new FragmentIonStruct(crosslink, param));

		// generate all b and y for now
		for (int i = 1; i < crosslink.peptideA.getSequenceLength(false); i++) // add b and y from peptideA
		{
			allIons.push(new FragmentIonStruct('b', i, false, false, crosslink, param));
			allIons.push(new FragmentIonStruct('y', i, false, false, crosslink, param));
			allIons.push(new FragmentIonStruct('b', i, true, false, crosslink, param));
			allIons.push(new FragmentIonStruct('y', i, true, false, crosslink, param));
		}

		for (int i = 1; i < crosslink.peptideB.getSequenceLength(false); i++) // add b and y from peptideB
		{
			allIons.push(new FragmentIonStruct('b', i, false, true, crosslink, param));
			allIons.push(new FragmentIonStruct('y', i, false, true, crosslink, param));
			allIons.push(new FragmentIonStruct('b', i, true, true, crosslink, param));
			allIons.push(new FragmentIonStruct('y', i, true, true, crosslink, param));
		}

		return allIons;
	}
	
	// faster, more manual method
	public static Stack<FragmentIonStruct> fragmentManual(PeptideStruct peptide)
	{
		Stack<FragmentIonStruct> allIons = new Stack<FragmentIonStruct>();
		allIons.push(new FragmentIonStruct(peptide));
		
		// keep track of mass and neutral loss flag
		double currentMass = 0, peptideMass = MassInfo.getMass(peptide);
		boolean[] currentFlagB = {false, false, false};
		boolean[] currentFlagY = {false, false, false};
		int offset = 0, length = peptide.getSequenceLength(false);
		
		if (peptide.getFirstResidue() == '[')
			offset = 1; // skip the N-terminal
		
		currentMass = MassInfo.getResidueMass(peptide, offset);
		NeutralLossGenerator.updateNeutralLossFlags(peptide, offset, currentFlagB);
		NeutralLossGenerator.updateNeutralLossFlags(peptide, offset + length - 1, currentFlagY);
		
		for (int i = 1; i < length; i++) // add b and y
		{
			allIons.push(new FragmentIonStruct('b', i, false, false, currentMass, currentFlagB));
			allIons.push(new FragmentIonStruct('y', length - i, false, false, peptideMass - currentMass, currentFlagY));
			
			currentMass += MassInfo.getResidueMass(peptide, offset + i); // update mass and flags
			NeutralLossGenerator.updateNeutralLossFlags(peptide, offset + i, currentFlagB);
			NeutralLossGenerator.updateNeutralLossFlags(peptide, offset - i + length - 1, currentFlagY);
		}

		return allIons;
	}
	
	public static Stack<FragmentIonStruct> fragmentManual(CrosslinkStruct crosslink, ParamStruct param)
	{
		Stack<FragmentIonStruct> allIons = new Stack<FragmentIonStruct>();
		allIons.push(new FragmentIonStruct(crosslink, param));
		
		// keep track of mass and neutral loss flag
		double currentMassA = 0, currentMassB, peptideMassA = MassInfo.getMass(crosslink.peptideA), peptideMassB = MassInfo.getMass(crosslink.peptideB);
		boolean[] currentFlagBA = {false, false, false};
		boolean[] currentFlagYA = {false, false, false};
		boolean[] currentFlagBB = {false, false, false};
		boolean[] currentFlagYB = {false, false, false};
		boolean[] flagA = NeutralLossGenerator.getNeutralLossFlags(crosslink.peptideA);
		boolean[] flagB = NeutralLossGenerator.getNeutralLossFlags(crosslink.peptideB);
		int offsetA = 0, offsetB = 0, lengthA = crosslink.peptideA.getSequenceLength(false), lengthB = crosslink.peptideB.getSequenceLength(false);
		
		if (crosslink.peptideA.getFirstResidue() == '[')
			offsetA = 1; // skip the N-terminal
		if (crosslink.peptideB.getFirstResidue() == '[')
			offsetB = 1; // skip the N-terminal
		
		currentMassA = MassInfo.getResidueMass(crosslink.peptideA, offsetA);
		currentMassB = MassInfo.getResidueMass(crosslink.peptideB, offsetB);
		NeutralLossGenerator.updateNeutralLossFlags(crosslink.peptideA, offsetA, currentFlagBA);
		NeutralLossGenerator.updateNeutralLossFlags(crosslink.peptideA, offsetA + lengthA - 1, currentFlagYA);
		NeutralLossGenerator.updateNeutralLossFlags(crosslink.peptideB, offsetB, currentFlagBB);
		NeutralLossGenerator.updateNeutralLossFlags(crosslink.peptideB, offsetB + lengthB - 1, currentFlagYB);
		
		for (int i = 1; i < lengthA; i++) // add b and y of peptideA
		{
			allIons.push(new FragmentIonStruct('b', i, false, false, currentMassA, currentFlagBA)); // linear ions
			allIons.push(new FragmentIonStruct('y', lengthA - i, false, false, peptideMassA - currentMassA, currentFlagYA));
			
			allIons.push(new FragmentIonStruct('b', i, true, false, currentMassA + peptideMassB + param.crosslinker.deltaMass, or(currentFlagBA, flagB))); // crosslinked ions
			allIons.push(new FragmentIonStruct('y', lengthA - i, true, false, peptideMassA - currentMassA + peptideMassB + param.crosslinker.deltaMass, or(currentFlagYA, flagB)));
			
			currentMassA += MassInfo.getResidueMass(crosslink.peptideA, offsetA + i); // update mass and flags
			NeutralLossGenerator.updateNeutralLossFlags(crosslink.peptideA, offsetA + i, currentFlagBA);
			NeutralLossGenerator.updateNeutralLossFlags(crosslink.peptideA, offsetA - i + lengthA - 1, currentFlagYA);
		}
		
		for (int i = 1; i < lengthB; i++) // add b and y of peptideB
		{
			allIons.push(new FragmentIonStruct('b', i, false, true, currentMassB, currentFlagBB)); // linear ions
			allIons.push(new FragmentIonStruct('y', lengthB - i, false, true, peptideMassB - currentMassB, currentFlagYB));
			
			allIons.push(new FragmentIonStruct('b', i, true, true, currentMassB + peptideMassA + param.crosslinker.deltaMass, or(currentFlagBB, flagA))); // crosslinked ions
			allIons.push(new FragmentIonStruct('y', lengthB - i, true, true, peptideMassB - currentMassB + peptideMassA + param.crosslinker.deltaMass, or(currentFlagYB, flagA)));
			
			currentMassB += MassInfo.getResidueMass(crosslink.peptideB, offsetB + i); // update mass and flags
			NeutralLossGenerator.updateNeutralLossFlags(crosslink.peptideB, offsetB + i, currentFlagBB);
			NeutralLossGenerator.updateNeutralLossFlags(crosslink.peptideB, offsetB - i + lengthB - 1, currentFlagYB);
		}

		return allIons;
	}
	
	// or operation for boolean array
	public static boolean[] or(boolean[] a, boolean[] b)
	{
		boolean[] result = new boolean[a.length];
		
		for (int i = 0; i < a.length; i++)
			result[i] = a[i] || b[i];
		
		return result;
	}

	/*
	public static TreeMap<String, FragmentIonStruct> fragment(PeptideStruct peptide)
	{
		TreeMap<String, FragmentIonStruct> allIons = new TreeMap<String, FragmentIonStruct>();
		FragmentIonStruct tempfragment = new FragmentIonStruct(peptide);
		allIons.put("precursor", tempfragment);

		for (int i = 1; i < peptide.getSequenceLength(false); i++) // add b and y
		{
			tempfragment = new FragmentIonStruct('b', i, peptide);
			allIons.put(tempfragment.getName(), tempfragment);

			tempfragment = new FragmentIonStruct('y', i, peptide);
			allIons.put(tempfragment.getName(), tempfragment);
		}

		return allIons;
	}

	public static TreeMap<String, FragmentIonStruct> fragment(CrosslinkStruct crosslink, ParamStruct param)
	{
		TreeMap<String, FragmentIonStruct> allIons = new TreeMap<String, FragmentIonStruct>();
		FragmentIonStruct tempfragment = new FragmentIonStruct(crosslink, param);
		allIons.put("precursor", tempfragment);

		// generate all b and y for now
		for (int i = 1; i < crosslink.peptideA.getSequenceLength(false); i++) // add b and y from peptideA
		{
			tempfragment = new FragmentIonStruct('b', i, false, false, crosslink, param);
			allIons.put(tempfragment.getName(), tempfragment);

			tempfragment = new FragmentIonStruct('y', i, false, false, crosslink, param);
			allIons.put(tempfragment.getName(), tempfragment);

			tempfragment = new FragmentIonStruct('b', i, true, false, crosslink, param);
			allIons.put(tempfragment.getName(), tempfragment);

			tempfragment = new FragmentIonStruct('y', i, true, false, crosslink, param);
			allIons.put(tempfragment.getName(), tempfragment);
		}

		for (int i = 1; i < crosslink.peptideB.getSequenceLength(false); i++) // add b and y from peptideB
		{
			tempfragment = new FragmentIonStruct('b', i, false, true, crosslink, param);
			allIons.put(tempfragment.getName(), tempfragment);

			tempfragment = new FragmentIonStruct('y', i, false, true, crosslink, param);
			allIons.put(tempfragment.getName(), tempfragment);

			tempfragment = new FragmentIonStruct('b', i, true, true, crosslink, param);
			allIons.put(tempfragment.getName(), tempfragment);

			tempfragment = new FragmentIonStruct('y', i, true, true, crosslink, param);
			allIons.put(tempfragment.getName(), tempfragment);
		}

		return allIons;
	}
	*/
}