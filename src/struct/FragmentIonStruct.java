// data structure for fragmented ion
public class FragmentIonStruct
{
	public final char ionType; // b or y
	public final int ionID; // 1, 2, ...
	public final boolean isCrosslink;
	public final boolean fromSecondPeptide; // true only if derived from second peptide of a crosslink
	public final double mass;
	public final String neutralLossTag;
	public final boolean[] neutralLossFlag; // flag whether this ion can lose [water, ammonia, metox-specific]

	public FragmentIonStruct(char ionType, int ionID, PeptideStruct peptide)
	{
		this.ionType = ionType;
		this.ionID = ionID;
		isCrosslink = false;
		fromSecondPeptide = false;
		neutralLossTag = "";
		mass = MassInfo.getFragmentMass(peptide, ionType, ionID);
		neutralLossFlag = NeutralLossGenerator.getNeutralLossFlags(peptide, ionType, ionID);
	}
	
	// initialize with mass and neutral loss flag pre-computed
	public FragmentIonStruct(char ionType, int ionID, boolean isCrosslink, boolean fromSecondPeptide, double mass, boolean[] neutralLossFlag)
	{
		this.ionType = ionType;
		this.ionID = ionID;
		this.isCrosslink = isCrosslink;
		this.fromSecondPeptide = fromSecondPeptide;
		neutralLossTag = "";
		this.mass = mass;
		this.neutralLossFlag = neutralLossFlag;
	}

	// precursor ion
	public FragmentIonStruct(PeptideStruct peptide)
	{
		this.ionType = 'p';
		this.ionID = -1;
		isCrosslink = false;
		fromSecondPeptide = false;
		neutralLossTag = "";
		mass = MassInfo.getMass(peptide);
		neutralLossFlag = NeutralLossGenerator.getNeutralLossFlags(peptide);
	}

	public FragmentIonStruct(char ionType, int ionID, boolean isCrosslink, boolean fromSecondPeptide, CrosslinkStruct crosslink, ParamStruct param)
	{
		this.ionType = ionType;
		this.ionID = ionID;
		this.isCrosslink = isCrosslink;
		this.fromSecondPeptide = fromSecondPeptide;
		neutralLossTag = "";
		mass = MassInfo.getFragmentMass(crosslink, ionType, ionID, isCrosslink, fromSecondPeptide, param);
		neutralLossFlag = NeutralLossGenerator.getNeutralLossFlags(crosslink, ionType, ionID, isCrosslink, fromSecondPeptide);
	}

	// precursor ion
	public FragmentIonStruct(CrosslinkStruct crosslink, ParamStruct param)
	{
		this.ionType = 'p';
		this.ionID = -1;
		this.isCrosslink = true;
		this.fromSecondPeptide = false;
		neutralLossTag = "";
		mass = MassInfo.getMass(crosslink, param);
		neutralLossFlag = NeutralLossGenerator.getNeutralLossFlags(crosslink);
	}

	// generate neutral loss version
	// add only one neutral loss at a time
	public FragmentIonStruct(FragmentIonStruct source, String neutralLossTag)
	{
		this.ionType = source.ionType;
		this.ionID = source.ionID;
		this.isCrosslink = source.isCrosslink;
		this.fromSecondPeptide = source.fromSecondPeptide;
		this.neutralLossFlag = source.neutralLossFlag;
		this.neutralLossTag = source.neutralLossTag + "-" + neutralLossTag;

		switch (neutralLossTag)
		{
			case "H2O":
				this.mass = source.mass - MassInfo.water;
				break;

			case "NH3":
				this.mass = source.mass - MassInfo.ammonia;
				break;

			case "OHSCH3":
				this.mass = source.mass - MassInfo.metOxLoss;
				break;

			case "CO":
				this.mass = source.mass - MassInfo.aIonLoss;
				break;
				
			case "C13":
				this.mass = source.mass + MassInfo.neutron;
				break;

			default:
				this.mass = -1;
				HelperFunctions.debug("Unexpected neutral loss!");
				break;
		}
	}
	
	// indicator of whether the ion is unmodified b or y ion
	public boolean isMajorIon()
	{ return (ionType != 'p') && neutralLossTag.equals(""); }

	// indicator of whether the ion is unmodified a, b, or y ion
	public boolean isUnmodifiedABYIon()
	{ return (ionType != 'p') && (neutralLossTag.equals("") || neutralLossTag.equals("-CO")); }

	// indicator of whether the ion is unmodified precursor
	public boolean isPrecursor()
	{ return (ionType == 'p') && neutralLossTag.equals(""); }

	// indicator of whether the ion can result in a-ion
	// b -> a has priority over other modification
	public boolean canHaveAion(ParamStruct param)
	{ return (ionType == 'b' && neutralLossTag.equals("") && ionID <= param.numAion); }

	// construct unique name for this fragmented ion
	public String getName()
	{
		String name;

		if (ionType == 'p') // precursor ion
			return "precursor" + neutralLossTag;

		if (fromSecondPeptide)
			name = "B" + ionType + "-" + ionID;
		else
			name = "A" + ionType + "-" + ionID;

		if (isCrosslink)
			name += "-linked";

		name += neutralLossTag;

		return name;
	}
	
	/*
	// equals method for comparing fragments in a set environment, i.e. Ab-1 portion only
	public boolean equals(FragmentIonStruct another)
	{ return (ionType == another.ionType) && (ionID == another.ionID) && (fromSecondPeptide == another.fromSecondPeptide) && (isCrosslink == another.isCrosslink); }
	*/
	
	// return major tag only
	public String getMajorTag()
	{
		String name;

		if (ionType == 'p') // precursor ion
			return "precursor";

		if (fromSecondPeptide)
			name = "B" + ionType + "-" + ionID;
		else
			name = "A" + ionType + "-" + ionID;
		
		return name;
	}

	public String toString()
	{ return "name: " + getName() + ", mass: " + mass; }
}