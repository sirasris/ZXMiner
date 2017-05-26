import java.util.TreeMap;
import java.util.ArrayList;
import java.util.Iterator;

// data structure for linear peptides
public class PeptideStruct
{
	public final ProteinStruct parent;
	public final int startPos, endPos;
	public final TreeMap<Integer, ModificationStruct> varModMap;

	public PeptideStruct(ProteinStruct parent, int startPos, int endPos)
	{
		this.parent = parent;
		this.startPos = startPos;
		this.endPos = endPos;
		varModMap = new TreeMap<Integer, ModificationStruct>();
	}

	public PeptideStruct(ProteinStruct parent, int startPos, int endPos, TreeMap<Integer,ModificationStruct> varModMap)
	{
		this.parent = parent;
		this.startPos = startPos;
		this.endPos = endPos;
		this.varModMap = new TreeMap<Integer, ModificationStruct>(varModMap); // TreeMap can be safely cloned
	}

	public PeptideStruct(ArrayList<ProteinStruct> proteins, int[] data)
	{
		parent = proteins.get(data[0]);
		startPos = data[1];
		endPos = data[2];
		varModMap = new TreeMap<Integer, ModificationStruct>();

		for (int i = 3; i < data.length; i += 2)
		{
			if (data[i] > -1)
				varModMap.put(new Integer(data[i]), MassInfo.getModification(data[i + 1]));
			else
				break;
		}
	}

	public PeptideStruct clone()
	{ return new PeptideStruct(this.parent, this.startPos, this.endPos, this.varModMap); }

	// add variable modification
	public void addModification(int pos, ModificationStruct mod)
	{
		if (varModMap.containsKey(new Integer(pos)))
			HelperFunctions.debug("PeptideStruct", "Multiple variable modifications assigned to same residue!");
		else
			varModMap.put(new Integer(pos), mod);
	}

	// return number of modifications
	public int getNumVarMod()
	{ return varModMap.size(); }

	// return total delta mass of all modifications
	public double getModDeltaMass()
	{
		double mass = 0;

		for (Iterator<ModificationStruct> iter = varModMap.values().iterator(); iter.hasNext();)
			mass += iter.next().deltaMass;

		return mass;
	}
	
	public double getModDeltaMass(int position)
	{
		if (isModified(position))
			return varModMap.get(new Integer(position)).deltaMass;
		else
			return 0;
	}

	// return total delta mass of all modifications
	public double getModDeltaMass(int startPos, int endPos, boolean includeTerminus)
	{
		double mass = 0;
		int offset = 0;
		Integer position;

		if (!includeTerminus && getFirstResidue() == '[')
			offset = 1;

		for (Iterator<Integer> iter = varModMap.keySet().iterator(); iter.hasNext();)
		{
			position = iter.next();

			if (position.intValue() >= startPos + offset && position.intValue() <= endPos + offset)
				mass += varModMap.get(position).deltaMass;
		}

		return mass;
	}

	// return whether a residue is modified
	public boolean isModified(int location)
	{ return varModMap.containsKey(new Integer(location)); }

	// return sequence length
	public int getSequenceLength(boolean includeTerminus)
	{ return parent.getSequenceLength(includeTerminus, startPos, endPos); }

	// return full sequence
	public String getSequence(boolean includeTerminus)
	{ 
		if (includeTerminus)
			return parent.getSequence(startPos, endPos);
		else
			return parent.getSequence(startPos, endPos).replaceAll("[\\W]", ""); // remove non-letter symbols
	}

	// return amino acid sequence
	public String getSequence(int startPos, int endPos, boolean includeTerminus)
	{
		if (!includeTerminus && getFirstResidue() == '[') // N-terminal symbol shift index by 1
			return parent.getSequence(this.startPos + startPos + 1, this.startPos + endPos + 1);
		else
			return parent.getSequence(this.startPos + startPos, this.startPos + endPos);
	}

	// return specific residue
	public char getResidue(int location)
	{ return parent.getResidue(startPos + location); }
	
	// return first residue
	// for checking for '['
	public char getFirstResidue()
	{ return parent.getResidue(startPos); }
	
	// return last residue
	// for checking for ']'
	public char getLastResidue()
	{ return parent.getResidue(endPos); }

	public boolean isForward()
	{ return parent.isForward; }
	
	public String toString()
	{
		String details = "parentProtein: " + parent.toString();
		details += ", startPos: " + startPos;
		details += ", endPos: " + endPos;
		details += ", sequence: " + getSequence(true);
		details += ", varModMap: " + varModMap.toString();

		return details;
	}
	
	public String toReport()
	{
		String posReport;
		
		if (getFirstResidue() == '[')
			posReport = "1\t";
		else
			posReport = startPos + "\t";
		
		if (getLastResidue() == ']')
			posReport += (endPos - 1) + "";
		else
			posReport += endPos + "";
		
		return parent.toReport() + "\t" + posReport + "\t" + getSequence(true) + "\t" + varModMap.toString();
	}
}