// data structure for matched fragmented ion
public class MatchedFragmentStruct
{
	public final double massError; // same type as specificed in ParamStruct
	public final FragmentIonStruct ion;
	public final int chargeState;

	public MatchedFragmentStruct(double massError, FragmentIonStruct ion, int chargeState)
	{
		this.massError = massError;
		this.ion = ion;
		this.chargeState = chargeState;
	}

	public boolean isUnmodifiedABYIon()
	{ return ion.isUnmodifiedABYIon(); }
	
	public boolean isMajorIon()
	{ return ion.isMajorIon(); }
	
	public String toString()
	{
		String details = "massError: " + massError;
		details += ", FragmentIonStruct: " + ion.toString();
		details += ", chargeState: " + chargeState;
		
		return details;
	}
}