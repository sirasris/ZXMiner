import java.util.Comparator;

// comparator for sorting 'MatchedFragmentStruct'
public class MatchedFragmentComparator implements Comparator<MatchedFragmentStruct>
{
	public MatchedFragmentComparator()
	{}

	public int compare(MatchedFragmentStruct m1, MatchedFragmentStruct m2)
	{
		if (m1.isMajorIon())
		{
			if (m2.isMajorIon())
				return compareMassError(m1, m2);
			else
				return -1; // prefer major ion explanation
		}

		else
		{
			if (m2.isMajorIon())
				return 1; // prefer major ion explanation
			else
				return compareMassError(m1, m2);
		}
	}

	// sort by absolute mass error from low to high
	public int compareMassError(MatchedFragmentStruct m1, MatchedFragmentStruct m2)
	{ return (int) Math.signum(Math.abs(m1.massError) - Math.abs(m2.massError)); }
}