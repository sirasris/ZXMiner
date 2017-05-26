import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// data structure for proteins
public class ProteinStruct
{
	public final int entryID;
	public final String name, description, sequence;
	public final int length;
	public final boolean isForward; // flag is TRUE for 'forward' and FALSE for 'reverse' and 'unrelated'
	ArrayList<Integer> cleaveSites = new ArrayList<Integer>(); // contain cleavage sites, starting index = 1

	public ProteinStruct(int entryID, String name, String description, String sequence, boolean isForward)
	{
		this.entryID = entryID;
		this.name = name;
		this.description = description;
		this.sequence = "[" + sequence + "]"; // add terminus symbols
		this.isForward = isForward;
		length = sequence.length(); // ignore terminus symbol
	}

	// compute all possible cleavage site for a given protease
	// convention is to extract sequence from 'cleaveSites[i] + 1' to 'cleaveSites[i + 1]' inclusive
	public void computeCleaveSites(ProteaseStruct protease)
	{
		cleaveSites = new ArrayList<Integer>();
		Pattern p = Pattern.compile(protease.getPatternString());
		Matcher m = p.matcher(sequence);
		int currentPos = 0; // keep track of the amino acid position to start the search

		cleaveSites.add(new Integer(-1)); // add N-terminal

		if (protease.isNterm) // protease cleave N-terminal
		{
			if (m.find(currentPos))
			{
				if (!(m.start() == 1)) // the situation where '-1' will be duplicated
					cleaveSites.add(new Integer(m.start() - 1));

				currentPos = m.start() + 1;
			}

			while (m.find(currentPos))
			{
				cleaveSites.add(new Integer(m.start() - 1));
	      		currentPos = m.start() + 1;
	    	}

	    	cleaveSites.add(new Integer(length + 1)); // add C-terminal
		}

		else // protease cleave C-terminal
		{
			while (m.find(currentPos))
			{
				cleaveSites.add(new Integer(m.start()));
	      		currentPos = m.start() + 1;
	    	}

			if (cleaveSites.get(cleaveSites.size() - 1).intValue() == length) // situation where 'length' will be duplicated
				cleaveSites.set(cleaveSites.size() - 1, new Integer(length + 1));
			else
				cleaveSites.add(new Integer(length + 1));
		}

		// HelperFunctions.debug("cleaveSites", cleaveSites);
	}

	public ArrayList<Integer> getCleaveSites()
	{ return cleaveSites; }

	// return number of starting site for peptides on this protein
	public int getNumStartSites()
	{ return cleaveSites.size() - 1; }

	// return sequence at specified indices [startPos, startPos + 1, ..., endPos]
	public String getSequence(int startPos, int endPos)
	{ return sequence.substring(startPos, endPos + 1); }

	// return sequence length
	public int getSequenceLength(boolean includeTerminus, int startPos, int endPos)
	{
		if (includeTerminus)
			return endPos + 1 - startPos;
		else
		{
			int result = endPos + 1 - startPos;

			if (startPos == 0)
				result--;

			if (endPos == length + 1)
				result--;

			return result;
		}
	}

	// return specific residue
	public char getResidue(int location)
	{ return sequence.charAt(location); }

	public String toString()
	{
		String details = "name: " + name;
		details += ", description: " + description;
		details += ", length: " + length;
		details += ", isForward: " + isForward;

		return details;
	}
	
	public String toReport()
	{ return name + "\t" + description + "\t" + isForward; }
}