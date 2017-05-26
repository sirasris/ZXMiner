// data structure for crosslinkers
public class CrosslinkerStruct
{
	public final boolean hasDeadEnd; // whether deadend and loop products are possible
	public final String siteA, siteB; // residues at the first and second sites
								      // assume that dead-end product only come from siteA
	public final double deltaMass; // mass added by crosslinker in crosslinked product

	public CrosslinkerStruct(boolean hasDeadEnd, String siteA, String siteB, double deltaMass)
	{
		this.hasDeadEnd = hasDeadEnd;
		this.siteA = siteA;
		this.siteB = siteB;
		this.deltaMass = deltaMass;
	}

	public String toString()
	{
		String details = "hasDeadEnd: " + hasDeadEnd;
		details += ", siteA: " + siteA;
		details += ", siteB: " + siteB;
		details += ", deltaMass: " + deltaMass;

		return details;
	}
	
	public String getCrosslinkRule(char rule)
	{
		if (rule == 'A')
			return siteA;
		else
			return siteB;
	}
}