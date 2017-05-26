import java.util.TreeMap;
import java.util.ArrayList;

public class MassInfo
{
	// frequently used mass values
	static final double neutron = 1.003355; // for heavy isotope mass calculation
	static final double hydrogen = 1.007825;
	static final double proton = 1.007276;
	static final double water = 18.010565;
	static final double ammonia = 17.026549;
	static final double aIonLoss = 27.994915;
	static final double metOxLoss = 63.999978;

	// all modifications
	static ArrayList<ModificationStruct> modifications = new ArrayList<ModificationStruct>();

	// atom mass table
	static TreeMap<Character, Double> atomMasses = new TreeMap<Character, Double>()
	{{
		put(new Character('C'), new Double(12));
		put(new Character('H'), new Double(1.007825));
		put(new Character('O'), new Double(15.994915));
		put(new Character('N'), new Double(14.003074));
		put(new Character('S'), new Double(31.972072));
		put(new Character('P'), new Double(30.973763));
	}};

	// amino acid mass table
	static TreeMap<Character, Double> aminoacidMasses = new TreeMap<Character, Double>()
	{{
		put(new Character('G'), new Double(57.02146));
		put(new Character('A'), new Double(71.03711));
		put(new Character('S'), new Double(87.03203));
		put(new Character('P'), new Double(97.05276));
		put(new Character('V'), new Double(99.06841));
		put(new Character('T'), new Double(101.04768));
		put(new Character('C'), new Double(103.00919));
		put(new Character('L'), new Double(113.08406));
		put(new Character('I'), new Double(113.08406));
		put(new Character('N'), new Double(114.04293));
		put(new Character('D'), new Double(115.02694));
		put(new Character('Q'), new Double(128.05858));
		put(new Character('K'), new Double(128.09496));
		put(new Character('E'), new Double(129.04259));
		put(new Character('M'), new Double(131.04049));
		put(new Character('H'), new Double(137.05891));
		put(new Character('F'), new Double(147.06841));
		put(new Character('R'), new Double(156.10111));
		put(new Character('Y'), new Double(163.06333));
		put(new Character('W'), new Double(186.07931));
		put(new Character('['), new Double(0));
		put(new Character(']'), new Double(0));
	}};

	// amino acid modifications
	static TreeMap<Character, ModificationStruct> fixedModifications = new TreeMap<Character, ModificationStruct>(); // at most one modification
	static TreeMap<Character, ArrayList<ModificationStruct>> varModifications = new TreeMap<Character, ArrayList<ModificationStruct>>();

	public MassInfo()
	{}

	// check if a residue has fixed modification
	public static boolean hasFixedMod(char residue)
	{ return fixedModifications.containsKey(new Character(residue)); }

	// return list of variable modifications for a residue
	public static ArrayList<ModificationStruct> getAllVarMod(char residue)
	{
		if (varModifications.containsKey(new Character(residue)))
			return varModifications.get(new Character(residue));
		else
			return null;
	}

	// return mass of target molecule
	public static double getMass(String molecule, boolean isAminoAcid)
	{
		double mass = 0;

		if (isAminoAcid) // determine which mass table to use
		{
			for (int i = 0; i < molecule.length(); i++)
				mass += aminoacidMasses.get(new Character(molecule.charAt(i))).doubleValue();

			mass += water;
		}

		else
		{
			for (int i = 0; i < molecule.length(); i++)
				mass += atomMasses.get(new Character(molecule.charAt(i))).doubleValue();
		}

		return mass;
	}
	
	// return pure residue mass
	public static double getResidueMass(char residue)
	{ return aminoacidMasses.get(new Character(residue)).doubleValue(); }
	
	// return residue + modification mass
	public static double getResidueMass(PeptideStruct peptide, int position)
	{ return getResidueMass(peptide.getResidue(position)) + peptide.getModDeltaMass(position); }

	// return mass of peptide
	public static double getMass(PeptideStruct peptide)
	{
		if (peptide == null) // trick for loop crosslink
			return 0;
		else
			return getMass(peptide.getSequence(true), true) + peptide.getModDeltaMass();
	}

	// return mass of peptide
	// only consider certain range, inclusive
	public static double getMass(PeptideStruct peptide, int startPos, int endPos, boolean includeTerminal)
	{
		if (peptide == null) // trick for loop crosslink
			return 0;
		else
			return getMass(peptide.getSequence(startPos, endPos, includeTerminal), true) + peptide.getModDeltaMass(startPos, endPos, includeTerminal);
	}

	// return mass of crosslink
	public static double getMass(CrosslinkStruct crosslink, ParamStruct param)
	{ return getMass(crosslink.peptideA) + getMass(crosslink.peptideB) + param.crosslinker.deltaMass; }

	// return mass of fragmented ion of linear peptide
	public static double getFragmentMass(PeptideStruct peptide, char ionType, int ionID)
	{
		double mass = 0;

		switch (ionType) // note that getMass() already add 'water'
		{
			case 'b':
				return getMass(peptide, 0, ionID - 1, false) - water;

			case 'y':
				return getMass(peptide, peptide.getSequenceLength(false) - ionID, peptide.getSequenceLength(false) - 1, false);

			default:
			{
				HelperFunctions.debug("Unexpected fragmented ion types");
				break;
			}
		}

		return mass;
	}

	// return mass of fragmented ion of crosslink
	// do not handle loop crosslink yet
	public static double getFragmentMass(CrosslinkStruct crosslink, char ionType, int ionID, boolean isCrosslink, boolean fromSecondPeptide, ParamStruct param)
	{
		double mass;

		if (fromSecondPeptide)
			mass = getFragmentMass(crosslink.peptideB, ionType, ionID);
		else
			mass = getFragmentMass(crosslink.peptideA, ionType, ionID);

		if (isCrosslink)
		{
			mass += param.crosslinker.deltaMass;

			if (fromSecondPeptide)
				mass += getMass(crosslink.peptideA);
			else
				mass += getMass(crosslink.peptideB);
		}

		return mass;
	}

	// return certain modification
	public static ModificationStruct getModification(int id)
	{ return modifications.get(id); }

	// add new entry modification table
	// input: "residue,name,mass", ex. "C,Carboxyamidomethyl,57.02146"
	public static void addModificationEntry(String input, boolean isVariable)
	{
		String[] content = input.split(",");
		Character aminoAcid = new Character(content[0].toUpperCase().charAt(0));
		ModificationStruct mod = new ModificationStruct(modifications.size(), content[1].toLowerCase(), Double.valueOf(content[2]));
		ArrayList<ModificationStruct> temp;

		if (isVariable) // variable modification
		{
			if (varModifications.containsKey(aminoAcid)) // check for existing amino acid in the map
			{
				temp = varModifications.get(aminoAcid);
				temp.add(mod);
			}

			else // create new entry
			{
				temp = new ArrayList<ModificationStruct>();
				temp.add(mod);
				varModifications.put(aminoAcid, temp);
			}

			modifications.add(mod); // add to global modification table
		}

		else // fixed modification
		{
			if (fixedModifications.containsKey(aminoAcid)) // complain
				HelperFunctions.debug("MassInfo", "Multiple fixed modification assigned to an amino acid!");
			else
			{
				double modMass = aminoacidMasses.get(aminoAcid).doubleValue();
				fixedModifications.put(aminoAcid, mod); // only record the first modification
				aminoacidMasses.put(aminoAcid, new Double(modMass + mod.deltaMass)); // update amino acid mass table
				modifications.add(mod); // add to global modification table
			}
		}
	}

	// add all modification from parameter
	public static void addModifications(ParamStruct param)
	{
		String[] tempst = param.fixedModification.split(";");

		for (int i = 0; i < tempst.length; i++) // add fixed modification
			addModificationEntry(tempst[i], false);

		tempst = param.varModification.split(";");

		for (int i = 0; i < tempst.length; i++) // add variable modification
			addModificationEntry(tempst[i], true);

		if (param.hasDeadEnd()) // consider dead-end products
		{
			for (int i = 0; i < param.crosslinker.siteA.length(); i++) // add as variable modification
				addModificationEntry(param.crosslinker.siteA.charAt(i) + ",Dead-end," + (param.crosslinker.deltaMass + hydrogen), true);
		}

		// HelperFunctions.debug("amino acid mass table", MassInfo.aminoacidMasses);
		// HelperFunctions.debug("all modifications", MassInfo.modifications);
	}
	
	// for clearing modification list
	public static void clearMemory()
	{
		fixedModifications.clear();
		varModifications.clear();
		modifications.clear();
		
		aminoacidMasses = new TreeMap<Character, Double>()
		{{
			put(new Character('G'), new Double(57.02146));
			put(new Character('A'), new Double(71.03711));
			put(new Character('S'), new Double(87.03203));
			put(new Character('P'), new Double(97.05276));
			put(new Character('V'), new Double(99.06841));
			put(new Character('T'), new Double(101.04768));
			put(new Character('C'), new Double(103.00919));
			put(new Character('L'), new Double(113.08406));
			put(new Character('I'), new Double(113.08406));
			put(new Character('N'), new Double(114.04293));
			put(new Character('D'), new Double(115.02694));
			put(new Character('Q'), new Double(128.05858));
			put(new Character('K'), new Double(128.09496));
			put(new Character('E'), new Double(129.04259));
			put(new Character('M'), new Double(131.04049));
			put(new Character('H'), new Double(137.05891));
			put(new Character('F'), new Double(147.06841));
			put(new Character('R'), new Double(156.10111));
			put(new Character('Y'), new Double(163.06333));
			put(new Character('W'), new Double(186.07931));
			put(new Character('['), new Double(0));
			put(new Character(']'), new Double(0));
		}};
	}
}