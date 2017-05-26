import java.util.ArrayList;
import java.util.Comparator;
import java.util.Collections;

// data structure for MS/MS spectrum
public class SpectrumStruct
{
	public final String rawfileName;
	public final PrecursorStruct precursor;
	public final int parentScanNumber;
	public final double retentionTime, basePeakIntensity, precursorIntensity;
	public double[][] spectrum; // can be changed via de-isotoping

	public SpectrumStruct(String rawfileName, PrecursorStruct precursor, int parentScanNumber, double retentionTime, double basePeakIntensity, double precursorIntensity, double[][] spectrum)
	{
		this.rawfileName = rawfileName;
		this.precursor = precursor;
		this.parentScanNumber = parentScanNumber;
		this.retentionTime = retentionTime;
		this.basePeakIntensity = basePeakIntensity;
		this.precursorIntensity = precursorIntensity;
		this.spectrum = spectrum;
	}

	public void setSpectrum(double[][] spectrum)
	{ this.spectrum = spectrum; }

	// size of spectrum
	public int length()
	{ return spectrum.length; }
	
	// return peak intensity
	public double getIntensity(int position)
	{ return spectrum[position][1]; }
	
	// return peak count and total intensity for scoring IntCov and PeakCov
	public double[] getTotalPeakAndInt(ParamStruct param)
	{
		double[] output = {0, 0};
		
		for (int i = 0; i < spectrum.length; i++)
		{
			if (spectrum[i][1] > param.minPeakIntToScore)
			{
				output[0] += 1;
				output[1] += spectrum[i][1];
			}
		}
		
		return output;
	}

	// return charge peaks
	// use minimum peak intensity threshold from 'ParamStruct'
	public ArrayList<ChargedPeakStruct> getSortedChargedPeaks(ParamStruct param)
	{
		ArrayList<ChargedPeakStruct> sortedChargedPeaks = new ArrayList<ChargedPeakStruct>();

		if (spectrum[0].length == 2) // no charge state determined
		{
			for (int i = 0; i < spectrum.length; i++)
			{
				if (spectrum[i][1] > param.minPeakInt)
				{
					if (spectrum[i][0] < precursor.precursorMZ)
					{
						for (int z = 1; z <= precursor.chargeState; z++)
							sortedChargedPeaks.add(new ChargedPeakStruct(i, z, spectrum[i][0]));
					}

					else // m/z higher than precursor's
					{
						for (int z = 1; z < precursor.chargeState; z++)
							sortedChargedPeaks.add(new ChargedPeakStruct(i, z, spectrum[i][0]));
					}
				}
			}
		}

		else
		{
			for (int i = 0; i < spectrum.length; i++)
			{
				if (spectrum[i][1] > param.minPeakInt)
				{
					if (spectrum[i][2] == 0) // unknown charge state
					{
						if (spectrum[i][0] < precursor.precursorMZ)
						{
							for (int z = 1; z <= precursor.chargeState; z++)
								sortedChargedPeaks.add(new ChargedPeakStruct(i, z, spectrum[i][0]));
						}

						else // m/z higher than precursor's
						{
							for (int z = 1; z < precursor.chargeState; z++)
								sortedChargedPeaks.add(new ChargedPeakStruct(i, z, spectrum[i][0]));
						}
					}

					else
						sortedChargedPeaks.add(new ChargedPeakStruct(i, (int) spectrum[i][2], spectrum[i][0]));
				}
			}
		}

		Collections.sort(sortedChargedPeaks, new Comparator<ChargedPeakStruct>()
		{
			public int compare(ChargedPeakStruct entry1, ChargedPeakStruct entry2) // compare based on mass
    		{
    			if (entry1.massWithOneCharge > entry2.massWithOneCharge)
    				return 1;
    			if (entry1.massWithOneCharge < entry2.massWithOneCharge)
    				return -1;

    			return 0;
    		}
		});

		return sortedChargedPeaks;
	}

	public String toString()
	{
		String details = precursor.toString();
		details += ", parentScanNumber: " + parentScanNumber;
		details += ", retentionTime: " + retentionTime;
		details += ", basePeakIntensity: " + basePeakIntensity;
		details += ", precursorIntensity: " + precursorIntensity + "\n";

		if (spectrum[0].length == 2)
		{
			for (int i = 0; i < spectrum.length; i++)
				details += "m/z: " + spectrum[i][0] + ", intensity: " + spectrum[i][1] + "\n";
		}

		else
		{
			for (int i = 0; i < spectrum.length; i++)
				details += "m/z: " + spectrum[i][0] + ", intensity: " + spectrum[i][1] + ", charge: " + spectrum[i][2] + "\n";
		}

		return details;
	}
	
	// report format for text output via SpectrumMatcher
	public String toReport()
	{ return rawfileName + "\t" + precursor.toReport() + "\t" + parentScanNumber + "\t" + retentionTime + "\t" + basePeakIntensity + "\t" + precursorIntensity + "\t" + length(); }
}