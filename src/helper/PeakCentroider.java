import java.util.ArrayList;

// convert from profile to centroided data
public class PeakCentroider
{
	static final double zeroIntensity = 1; // numerical bound for 'practically-zero' intensity valuein full scan
	
	public PeakCentroider()
	{}
	
	public static double[][] centroided(double[][] spectrum)
	{
		ArrayList<double[]> localList = new ArrayList<double[]>();
		
		int currentStart = 0, currentEnd, currentIntPeak;
		
		while (currentStart < spectrum.length)
		{
			while (currentStart < spectrum.length && spectrum[currentStart][1] < zeroIntensity) // move across zero intensity peaks
				currentStart++;
			
			if (currentStart == spectrum.length) // terminate
				break;
			
			// all set with start position
			currentEnd = currentStart + 1;
			
			while (currentEnd < spectrum.length && spectrum[currentEnd][1] >= spectrum[currentEnd - 1][1]) // keep climbing a peak
				currentEnd++;
			
			currentIntPeak = currentEnd - 1; // location of the highest intensity peak
			
			if (currentEnd == spectrum.length) // record peak and terminate
			{
				localList.add(getCentroidData(spectrum, currentStart, currentEnd, currentEnd - 1));				
				break;
			}
			
			while (currentEnd < spectrum.length && spectrum[currentEnd][1] <= spectrum[currentEnd - 1][1] && spectrum[currentEnd][1] > zeroIntensity) // keep going down hill
				currentEnd++;

			localList.add(getCentroidData(spectrum, currentStart, currentEnd, currentIntPeak));
			
			if (currentEnd == spectrum.length) // terminate
				break;
			else
				currentStart = currentEnd; // continue
		}
		
		double[][] output = new double[localList.size()][2];
		
		for (int i = 0; i < localList.size(); i++)
			output[i] = localList.get(i);
		
		return output;
	}
	
	public static double[] getCentroidData(double[][] spectrum, int start, int end, int peak)
	{
		double[] result = new double[2];
		
		if (end - start == 1)
		{
			result[0] = spectrum[start][0];
			result[1] = spectrum[start][1];
			
			return result;
		}
		
		if (end - start == 2)
		{
			result[0] = (spectrum[start][0] * spectrum[start][1] + spectrum[start + 1][0] * spectrum[start + 1][1]) / (spectrum[start][1] + spectrum[start + 1][1]); // weighted average
			result[1] = spectrum[start][1] + spectrum[start + 1][1]; // sum
			
			return result;
		}
		
		// 3 or more data points
		int[] centers = new int[3];
		
		if (peak == spectrum.length - 1) // peak is at the end of the data
		{
			centers[0] = peak - 2;
			centers[1] = peak - 1;
			centers[2] = peak;
		}
		
		else
		{
			centers[0] = peak - 1;
			centers[1] = peak;
			centers[2] = peak + 1;
		}
		
		result[1] = 0;
		
		for (int i = start; i < end; i++)
			result[1] += spectrum[i][1]; // sum intensity
		
		double[] coeff = {Math.log(spectrum[centers[1]][1]) - Math.log(spectrum[centers[2]][1]), 
				Math.log(spectrum[centers[2]][1]) - Math.log(spectrum[centers[0]][1]), 
				Math.log(spectrum[centers[0]][1]) - Math.log(spectrum[centers[1]][1])};
		
		result[0] = 0.5 * (coeff[0] * Math.pow(spectrum[centers[0]][0], 2) + coeff[1] * Math.pow(spectrum[centers[1]][0], 2) + coeff[2] * Math.pow(spectrum[centers[2]][0], 2)) 
				/ (coeff[0] * spectrum[centers[0]][0] + coeff[1] * spectrum[centers[1]][0] + coeff[2] * spectrum[centers[2]][0]); // gaussian fit
		
		return result;
	}
}
