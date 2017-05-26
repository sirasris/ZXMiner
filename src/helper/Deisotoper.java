import java.util.ArrayList;
import java.util.TreeSet;
import java.util.Arrays;
import java.util.Comparator;

// de-isotope algorithm similar to that described in Senko M. et al.
public class Deisotoper
{
	static final double[] averagineRatio = {4.9384, 7.7583, 1.3577, 1.4773, 0.0417}; // ratio for C, H, N, O, S

	// final double[] chiSquareCutoff = {0.06, 0.45, 1.01, 1.65, 2.34, 3.07, 3.82, 4.59, 5.38, 6.18}; // cutoff at prob chance = 0.8
	// static final double[] chiSquareCutoff = {0.15, 0.71, 1.42, 2.20, 3.00, 3.83, 4.67, 5.53, 6.39, 7.27}; // cutoff at prob chance = 0.7
	// final double[] chiSquareCutoff = {0.27, 1.02, 1.87, 2.75, 3.66, 4.57, 5.49, 6.42, 7.36, 8.30}; // cutoff at prob chance = 0.6
	// final double[] chiSquareCutoff = {0.46, 1.39, 2.37, 3.36, 4.35, 5.35, 6.35, 7.34, 8.34, 9.34}; // cutoff at prob chance = 0.5

	// static final int maxEnvelopeLength = chiSquareCutoff.length;
	static final int maxEnvelopeLength = 10;
	static final int minEnvelopeLength = 2;
	static final double minAveragineCorr = 0.6; // quality criteria for matching to averagine's isotope pattern
	
	static final double averagineWeight = averagineRatio[0] * MassInfo.getMass("C", false) + averagineRatio[1] * MassInfo.getMass("H", false) +
    									  averagineRatio[2] * MassInfo.getMass("N", false) + averagineRatio[3] * MassInfo.getMass("O", false) +
    									  averagineRatio[4] * MassInfo.getMass("S", false); // averagine molecular weight

	public Deisotoper()
	{}

	public static double[][] deisotope(double[][] spectrum, int precursorCharge, ParamStruct param)
	{
		double currentMZ, nextMZ;
    	int nextPeak;
    	boolean[] assigned = new boolean[spectrum.length]; // keep track of peaks that are already assigned to an envelope at each charge state level
    	boolean foundEnvelope;
    	IsotopeEnvelopeStruct currentCandidate;

    	TreeSet<IsotopeEnvelopeStruct> candidateEnvelopes = new TreeSet<IsotopeEnvelopeStruct>(new IsotopeEnvelopeComparator()); // sorted list of candidate envelopes

    	// for each charge state, we generate candidate envelopes of the form "chargeState,peak1,peak2,...,peakN"
    	for (int charge = precursorCharge; charge > 0; charge--)
    	{
    		// reset "assigned" indicator
    		for (int i = 0; i < spectrum.length; i++)
    			assigned[i] = false;

    		// for each starting peak
    		for (int currentPeak = 0; currentPeak < spectrum.length - 1; currentPeak++)
    		{
    			foundEnvelope = false; // starting a new search for candidate envelope

    			// continue only if the current peak hasn't been used for this charge state
    			// no need to update status for 'currentPeak' since we move forward only
    			if (!assigned[currentPeak])
    			{
    				// initialize the candidate string format and m/z values
    				currentCandidate = new IsotopeEnvelopeStruct(charge, currentPeak);
    				currentMZ = spectrum[currentPeak][0];
    				nextPeak = currentPeak + 1;
    				nextMZ = spectrum[nextPeak][0];

    				// repeat until the envelope can't be extended (no m/z within window or end of spectrum table)
    				while (currentCandidate.size() < maxEnvelopeLength && (nextMZ - currentMZ < (MassInfo.neutron / charge) + (currentMZ * param.isotopeWindowppm * 0.000001)))
					{
						// if nextMZ differ from currentMZ within reasonable error from (MassInfo.neutron / charge), then we extend the envelope
						if (nextMZ - currentMZ > (MassInfo.neutron / charge) - (currentMZ * param.isotopeWindowppm * 0.000001))
						{
							// move currentMZ along and update envelope information
							foundEnvelope = true;
							assigned[nextPeak] = true;
							currentCandidate.add(nextPeak);
							currentMZ = nextMZ;
						}

						nextPeak++;

						// if we hit the end of the spectrum, terminate
						if (nextPeak == spectrum.length)
							break;

						nextMZ = spectrum[nextPeak][0];
					}

					// record the candidate envelope
					if (foundEnvelope && currentCandidate.size() >= minEnvelopeLength)
						candidateEnvelopes.add(currentCandidate);
    			}
    		}
    	}

    	ArrayList<IsotopeEnvelopeStruct> validatedEnvelopes = validateEnvelope(candidateEnvelopes, spectrum); // list of validated envelopes + unused peaks at the end
    	int numPeak = validatedEnvelopes.size() - 1 + validatedEnvelopes.get(validatedEnvelopes.size() - 1).size();
    	double intensity;
    	IsotopeEnvelopeStruct envelope;

    	double[][] deisotoped = new double[numPeak][3];

    	for (int i = 0; i < validatedEnvelopes.size() - 1; i++) // extract envelopes
    	{
			envelope = validatedEnvelopes.get(i);
			intensity = 0;

			for (int j = 0; j < envelope.size(); j++) // sum intensity
				intensity += spectrum[envelope.get(j)][1];

			deisotoped[i][0] = spectrum[envelope.get(0)][0]; // monoisotopic mz
			deisotoped[i][1] = intensity;
			deisotoped[i][2] = envelope.chargeState;
    	}

    	// peaks with no charge state
    	envelope = validatedEnvelopes.get(validatedEnvelopes.size() - 1);

    	for (int i = 0; i < envelope.size(); i++)
    	{
    		deisotoped[validatedEnvelopes.size() + i - 1][0] = spectrum[envelope.get(i)][0];
    		deisotoped[validatedEnvelopes.size() + i - 1][1] = spectrum[envelope.get(i)][1];
    	}

    	Arrays.sort(deisotoped, new Comparator<double[]>()
    	{
    		public int compare(double[] entry1, double[] entry2) // compare based on mz value
    		{
    			if (entry1[0] > entry2[0])
    				return 1;
    			if (entry1[0] < entry2[0])
    				return -1;

    			return 0;
    		}
    	});

		return deisotoped;
	}

	// validate candidate envelopes from largest to smallest
	// check against averagine intensity distribution
	// also make sure not to use the same peak multiple time
	public static ArrayList<IsotopeEnvelopeStruct> validateEnvelope(TreeSet<IsotopeEnvelopeStruct> candidateEnvelopes, double[][] spectrum)
	{
		ArrayList<IsotopeEnvelopeStruct> validatedEnvelopes = new ArrayList<IsotopeEnvelopeStruct>();
		double[] observedIntDistribution;
		double[] theoIntDistribution;
		double temp;
		boolean[] assigned = new boolean[spectrum.length]; // track peak usage

		IsotopeEnvelopeStruct envelope = candidateEnvelopes.pollFirst();
		IsotopeEnvelopeStruct[] tempenvelope;
		boolean broken;

		while (envelope != null) // repeat until all candidates were analyzed
		{
			broken = false;

			for (int i = 0; i < envelope.size(); i++) // check whether all peaks are available
			{
				if (assigned[envelope.get(i)])
				{
					tempenvelope = envelope.split(i);
					broken = true;

					if (tempenvelope[0] != null && tempenvelope[0].size() >= minEnvelopeLength)
						candidateEnvelopes.add(tempenvelope[0]);
					if (tempenvelope[1] != null && tempenvelope[0].size() >= minEnvelopeLength)
						candidateEnvelopes.add(tempenvelope[1]); // add smaller envelopes back in

					break;
				}
			}

			if (!broken) // continue only if passed
			{
				theoIntDistribution = averagineIntensity(envelope, spectrum); // obtain averagine distribution
				observedIntDistribution = new double[envelope.size()];
				temp = 0;

				for (int i = 0; i < envelope.size(); i++) // extract intensity
				{
					observedIntDistribution[i] = spectrum[envelope.get(i)][0];
					temp += observedIntDistribution[i];
				}

				for (int i = 0; i < envelope.size(); i++) // normalize
					observedIntDistribution[i] /= temp;

				//if (chiSquareValue(observedIntDistribution, theoIntDistribution) < chiSquareCutoff[envelope.size() - 1]) // passed
				if (correlation(observedIntDistribution, theoIntDistribution) >= minAveragineCorr)
				{
					validatedEnvelopes.add(envelope);

					for (int i = 0; i < envelope.size(); i++) // update peak status as 'used'
						assigned[envelope.get(i)] = true;
				}

				else // try either removing the first or last peaks
				{
					if (envelope.size() > minEnvelopeLength)
					{
						candidateEnvelopes.add(new IsotopeEnvelopeStruct(envelope, 0, envelope.size() - 2)); // remove last peak
						candidateEnvelopes.add(new IsotopeEnvelopeStruct(envelope, 1, envelope.size() - 1)); // remove first peak
					}
				}
			}
			
			envelope = candidateEnvelopes.pollFirst(); // retrieve next element
		}

		// add unused peaks as a zero-charge group
		envelope = new IsotopeEnvelopeStruct(0);

		for (int i = 0; i < spectrum.length; i++)
		{
			if (!assigned[i])
				envelope.add(i);
		}

		validatedEnvelopes.add(envelope);

		return validatedEnvelopes;
	}

	// compute chi-square score
	public static double chiSquareValue(double[] real, double[] expect)
	{
		double sum = 0;

		for (int i = 0; i < real.length; i++)
			sum += Math.pow(real[i] - expect[i], 2) / expect[i];

		return sum;
	}

	// compute ideal envelope's intensity ratio based on averagine
	public static double[] averagineIntensity(IsotopeEnvelopeStruct envelope, double[][] spectrum)
	{
		double mass = spectrum[envelope.get(0)][0] * envelope.chargeState - envelope.chargeState * MassInfo.proton;
		double fold = mass / averagineWeight;
		int[] newRatio = new int[5]; // rounded [C, H, N, O, S] ratio

		for (int i = 0; i < 5; i++)
			newRatio[i] = (int) Math.round(averagineRatio[i] * fold);

		double[] distribution = new double[envelope.size()];
		double sum = 0;

		for (int i = 0; i < envelope.size(); i++)
		{
			distribution[i] = isotopeProb(newRatio, i);
			sum += distribution[i];
		}

		for (int i = 0; i < envelope.size(); i++)
			distribution[i] /= sum;

		return distribution; // sum to 1
	}

	// compute P(isotope = +n neutrons) given [C, H, N, O, S] composition
	// consider only C13 (+1), N15 (+1), and O18 (+2)
	public static double isotopeProb(int[] comp, int n)
	{
		// if we don't have enough atom to generate the isotope, return zero
		if (comp[0] + comp[2] + 2*comp[3] < n)
			return 0;

		double[] isoProbs = {0.01109, 0.00366, 0.00201}; // C13, N15, and O18 existences in nature
		double prob = 0, temp;

		// iterate over all possible n-object distribution using "star & bar" technique
		for (int n1 = 1; n1 < n + 2; n1++)
		for (int n2 = n1 + 1; n2 < n + 3; n2++)
		{
			// there will be (n1 - 1) of C13, (n2 - n1 - 1) of N15, and (n + 3 - n2 - 1)/2 of O18 for a grand total of exactly n neutrons
			// proceed only if we have an integer number of O18 and does not exceed any atomic composition
			if ((n - n2) % 2 == 0 && n1 < comp[0] + 2 && n2 - n1 < comp[2] + 2 && (n - n2) / 2 < comp[3])
			{
				temp = logBinomial(comp[0], n1 - 1, isoProbs[0]) + logBinomial(comp[2], n2 - n1 - 1, isoProbs[1]) +
					logBinomial(comp[3], (n - n2) / 2 + 1, isoProbs[2]);

				prob += Math.exp(temp);
			}
		}

		return prob;
	}

	// compute log bionomial term (n choose k) (1 - a)^(n-k) a^k [x^k]
	public static double logBinomial(int n, int k, double a)
	{
		double result = 0;

		result += Math.log(1 - a) * (n - k) + Math.log(a) * k;

		if (k == 0)
			return result;

		for (int i = 1; i < k + 1; i++)
			result += Math.log(n + 1 - i) - Math.log(i);

		return result;
	}
	
	// mean data = 1 for both inputs
	// compute standard Pearson correlation
	public static double correlation(double[] observed, double[] expected)
	{
		double[][] temp = new double[observed.length][2];
		
		for (int i = 0; i < temp.length; i++)
		{
			temp[i][0] = observed[i] - (1.0 / temp.length);
			temp[i][1] = expected[i] - (1.0 / temp.length);
		}
		
		double num = temp[0][0] * temp[0][1];
		double denom1 = temp[0][0] * temp[0][0];
		double denom2 = temp[0][1] * temp[0][1];
		
		for (int i = 1; i < temp.length; i++)
		{
			num += temp[i][0] * temp[i][1];
			denom1 += temp[i][0] * temp[i][0];
			denom2 += temp[i][1] * temp[i][1];
		}
		
		return num / Math.sqrt(denom1 * denom2);
	}
}