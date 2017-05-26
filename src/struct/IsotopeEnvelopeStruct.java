import java.util.ArrayList;

// data structure for candidate isotopic envelope
public class IsotopeEnvelopeStruct
{
	public final int chargeState;
	public final ArrayList<Integer> peaks;

	public IsotopeEnvelopeStruct(int chargeState)
	{
		this.chargeState = chargeState;
		peaks = new ArrayList<Integer>();
	}

	public IsotopeEnvelopeStruct(int chargeState, int firstPeak)
	{
		this.chargeState = chargeState;
		peaks = new ArrayList<Integer>();
		peaks.add(new Integer(firstPeak));
	}

	// create new 'IsotopeEnvelopeStruct' with partial elements [startIndex, ..., endIndex] inclusive
	public IsotopeEnvelopeStruct(IsotopeEnvelopeStruct source, int startIndex, int endIndex)
	{
		chargeState = source.chargeState;
		peaks = new ArrayList<Integer>(endIndex - startIndex + 1);

		for (int i = startIndex; i < endIndex + 1; i++)
			peaks.add(source.get(i));
	}

	public void add(int nextPeak)
	{ peaks.add(new Integer(nextPeak)); }

	public int size()
	{ return peaks.size(); }

	public int get(int location)
	{ return peaks.get(location).intValue(); }

	public IsotopeEnvelopeStruct[] split(int location)
	{
		IsotopeEnvelopeStruct[] result = new IsotopeEnvelopeStruct[2];

		if (location > 1) // we have at least 2 peaks before 'location'
			result[0] = new IsotopeEnvelopeStruct(this, 0, location - 1);
		else
			result[0] = null;

		if (size() - location - 1 > 1) // we have at least 2 peaks after 'location'
			result[1] = new IsotopeEnvelopeStruct(this, location + 1, size() - 1);
		else
			result[1] = null;

		return result;
	}
}