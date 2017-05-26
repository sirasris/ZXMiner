// data structure for proteases
public class ProteaseStruct
{
	public final boolean isNterm; // direction of cleaving
	public final String restrictions; // residues that can prevent cleavage
	public final String sites; // residues at cleavage site

	public ProteaseStruct(boolean isNterm, String restrictions, String sites)
	{
		this.isNterm = isNterm;
		this.restrictions = restrictions;
		this.sites = sites;
	}

	public ProteaseStruct(boolean isNterm, String sites)
	{
		this.isNterm = isNterm;
		this.restrictions = null; // default
		this.sites = sites;
	}

	// return regular expression rule for this protease
	public String getPatternString()
	{
		String pattern = "[" + sites + "]";

		if (restrictions != null)
		{
			if (isNterm)
				return "[^" + restrictions + "]" + pattern;
			else
				return pattern + "[^" + restrictions + "]";
		}

		else
			return pattern;
	}

	public String toString()
	{
		String details = "isNterm: " + isNterm;
		details += ", restrictions: " + restrictions;
		details += ", sites: " + sites;

		return details;
	}
}