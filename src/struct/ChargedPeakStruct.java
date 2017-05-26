// data structure for m/z peak with specified charge state
public class ChargedPeakStruct
{
	public final int peakID;
	public final int chargeState;
	public final double massWithOneCharge;

	public ChargedPeakStruct(int peakID, int chargeState, double mz) // input m/z
	{
		this.peakID = peakID;
		this.chargeState = chargeState;
		this.massWithOneCharge = mz * chargeState - (chargeState - 1) * MassInfo.proton; // convert to mass with one charges
	}
}