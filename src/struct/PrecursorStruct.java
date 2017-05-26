// data structure for precursor information
public class PrecursorStruct
{
	public final int rawFileID;
	public final int scanNumber;
	public final int chargeState;
	public final double precursorMZ; // mass over charge

	public PrecursorStruct(int rawFileID, String[] inputString)
	{
		this.rawFileID = rawFileID;
		scanNumber = Integer.valueOf(inputString[1]);
		chargeState = Integer.valueOf(inputString[2]);
		precursorMZ = Double.valueOf(inputString[3]);
	}

	public String toString()
	{
		String details = "rawFileID: " + rawFileID;
		details += ", scanNumber: " + scanNumber;
		details += ", chargeState: " + chargeState;
		details += ", precursorMZ: " + precursorMZ;

		return details;
	}
	
	public double getPrecursorMHP()
	{ return (precursorMZ - MassInfo.proton) * chargeState + MassInfo.proton; }
	
	public String toReport()
	{ return scanNumber + "\t" + chargeState + "\t" + precursorMZ + "\t" + getPrecursorMHP(); }
}