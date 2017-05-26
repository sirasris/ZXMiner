import java.util.TreeMap;
import java.util.HashSet;
import java.util.Stack;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.Collections;

// match spectrum to theoretical peptide
// called as multi-thread via 'DataGrabber'
public class SpectrumMatcher implements Runnable
{
	double[] overallScores = new double[10]; // peak coverage, intensity coverage, ion coverage peptideA b-ion, ion coverage peptideA y-ion, ion coverage peptideA
											 // ion coverage peptideB b-ion, ion coverage peptideB y-ion, ion coverage peptideB, ion coverage, GM score
	double[] globalSpectrumStat; // get the total peak counts and total intensity
	int[] peptideLength = new int[2]; // pre-compute peptide length
	int[][] crosslinkSiteBound = new int[2][2]; // boundary of crosslink sites
	boolean isCrosslink;
	double precursorMassError;
	ParamStruct param;
	PeptideStruct peptide = null;
	CrosslinkStruct crosslink = null;
	Stack<FragmentIonStruct> fragmentedIons; // list of major b and y fragmented ions, only used once, can be edited safely
	SpectrumStruct currentSpectrum;
	TreeMap<Integer, ArrayList<MatchedFragmentStruct>> matchedMap = new TreeMap<Integer, ArrayList<MatchedFragmentStruct>>();

	public SpectrumMatcher(Stack<FragmentIonStruct> fragmentedIons, SpectrumStruct currentSpectrum, Object generalizedPeptide, ParamStruct param)
	{
		this.fragmentedIons = fragmentedIons;
		this.currentSpectrum = currentSpectrum;
		this.param = param;
		globalSpectrumStat = currentSpectrum.getTotalPeakAndInt(param);
		
		// HelperFunctions.debug("currentFragment", fragmentedIons);

		if (generalizedPeptide instanceof PeptideStruct)
		{
			peptide = (PeptideStruct) generalizedPeptide;
			peptideLength[0] = peptide.getSequenceLength(false);
			peptideLength[1] = -1;
			precursorMassError = HelperFunctions.getppmError(MassInfo.getMass(peptide) + MassInfo.proton, currentSpectrum.precursor.getPrecursorMHP());
			isCrosslink = false;
			// HelperFunctions.debug("peptide", peptide);
		}
		
		else
		{
			crosslink = (CrosslinkStruct) generalizedPeptide;
			peptideLength[0] = crosslink.peptideA.getSequenceLength(false);
			peptideLength[1] = crosslink.peptideB.getSequenceLength(false);
			precursorMassError = HelperFunctions.getppmError(MassInfo.getMass(crosslink, param) + MassInfo.proton, currentSpectrum.precursor.getPrecursorMHP());
			isCrosslink = true;
			// HelperFunctions.debug("crosslink", crosslink);
		}
		
		// HelperFunctions.debug("peptide length", peptideLength);
	}

	// report as target list, i.e. [m/z, RT]
	public String toTargetListEntry()
	{ return currentSpectrum.precursor.precursorMZ + "\t" + currentSpectrum.retentionTime + "\t" + precursorMassError; }
	
	// generate all neutral losses form of a major ion
	public void generateAllNeutralLosses(int step, FragmentIonStruct original, ArrayList<FragmentIonStruct> target)
	{
		boolean[] tempFlag = original.neutralLossFlag;

		if (step < tempFlag.length) // valid step
		{
			if (tempFlag[step])
			{
				String neutralLossName = NeutralLossGenerator.neutralLossName[step];
				int currentLength = target.size();

				if (neutralLossName.equals("H2O") || neutralLossName.equals("NH3")) // repeat according to 'ParamStruct'
				{
					int repeatCount = param.numNeutralLoss;

					if (original.isPrecursor())
						repeatCount = repeatCount * 2; // double for precursor

					for (int n = 0; n < repeatCount; n++)
					{
						for (int i = currentLength * n; i < currentLength * (n + 1); i++)
							target.add(new FragmentIonStruct(target.get(i), neutralLossName));

						// currentLength = target.size(); // update length of 'target'
					}
				}

				else // add only once for each entry in 'target'
				{
					for (int i = 0; i < currentLength; i++)
						target.add(new FragmentIonStruct(target.get(i), neutralLossName));
				}
			}

			generateAllNeutralLosses(step + 1, original, target); // continue
		}
	}

	// check whether a fragment ion agree with crosslink site
	public boolean checkFragmentWithLinkSite(FragmentIonStruct ion, int siteA, int siteB)
	{
		if (ion.ionType == 'p') // precursor
			return true;
		
		if (ion.fromSecondPeptide) // peptideB
		{
			if (ion.ionType == 'b') // b-ion
			{
				if (ion.ionID < siteB && !ion.isCrosslink) // non-linked
					return true;
				if (ion.ionID >= siteB && ion.isCrosslink) // linked
					return true;
			}
			
			else // y-ion
			{
				if (ion.ionID <= peptideLength[1] - siteB && !ion.isCrosslink) // non-linked
					return true;
				if (ion.ionID > peptideLength[1] - siteB && ion.isCrosslink) // linked
					return true;
			}
		}
		
		else // peptideA
		{
			if (ion.ionType == 'b') // b-ion
			{
				if (ion.ionID < siteA && !ion.isCrosslink) // non-linked
					return true;
				if (ion.ionID >= siteA && ion.isCrosslink) // linked
					return true;
			}
			
			else // y-ion
			{
				if (ion.ionID <= peptideLength[0] - siteA && !ion.isCrosslink) // non-linked
					return true;
				if (ion.ionID > peptideLength[0] - siteA && ion.isCrosslink) // linked
					return true;
			}
		}
		
		return false; // otherwise
	}
	
	// check whether a fragment ion agree with crosslink bounary
	public boolean checkFragmentWithLinkSite(FragmentIonStruct ion, int[][] boundary)
	{
		if (ion.ionType == 'p') // precursor
			return true;
		
		if (ion.fromSecondPeptide) // peptideB
		{
			if (ion.ionType == 'b') // b-ion
			{
				if (ion.ionID < boundary[1][0] && ion.isCrosslink) // linked
					return false;
				if (ion.ionID >= boundary[1][1] && !ion.isCrosslink) // non-linked
					return false;
			}
			
			else // y-ion
			{
				if (ion.ionID > peptideLength[1] - boundary[1][0] && !ion.isCrosslink) // non-linked
					return false;
				if (ion.ionID <= peptideLength[1] - boundary[1][1] && ion.isCrosslink) // linked
					return false;
			}
		}
		
		else // peptideA
		{
			if (ion.ionType == 'b') // b-ion
			{
				if (ion.ionID < boundary[0][0] && ion.isCrosslink) // linked
					return false;
				if (ion.ionID >= boundary[0][1] && !ion.isCrosslink) // non-linked
					return false;
			}
			
			else // y-ion
			{
				if (ion.ionID > peptideLength[0] - boundary[0][0] && !ion.isCrosslink) // non-linked
					return false;
				if (ion.ionID <= peptideLength[0] - boundary[0][1] && ion.isCrosslink) // linked
					return false;
			}
		}
		
		return true; // otherwise
	}
	
	// compute scores based on specific crosslink site
	public double[] computeScore(int siteA, int siteB)
	{
		double[] matchedStat = new double[2]; // count number of matched peak and total matched intensity
		double[] localScore = new double[overallScores.length];
		
		HashSet<String> majorMatchA = new HashSet<String>();
		HashSet<String> majorMatchB = new HashSet<String>(); // record major ion matches for each peptide
		
		ArrayList<MatchedFragmentStruct> tempMatchList;
		MatchedFragmentStruct bestMatch, bestMajorMatch;
		Integer peakID;
		
		for (Iterator<Integer> iter = matchedMap.keySet().iterator(); iter.hasNext();)
		{
			peakID = iter.next();
			tempMatchList = matchedMap.get(peakID);
			bestMatch = null;
			bestMajorMatch = null;
			
			for (int i = 0; i < tempMatchList.size(); i++) // look for the best match which agree with crosslink site interpretation
			{
				if (checkFragmentWithLinkSite(tempMatchList.get(i).ion, siteA, siteB))
				{
					bestMatch = tempMatchList.get(i);
					// HelperFunctions.debug("bestMatch" + peakID.toString(), bestMatch);
					break;  // only need the best one
				}
			}
			
			for (int i = 0; i < tempMatchList.size(); i++) // look for the best major ion match which agree with crosslink site interpretation
			{
				if (tempMatchList.get(i).isMajorIon() && checkFragmentWithLinkSite(tempMatchList.get(i).ion, siteA, siteB))
				{
					bestMajorMatch = tempMatchList.get(i);
					// HelperFunctions.debug("bestMatch" + peakID.toString(), bestMatch);
					break;  // only need the best one
				}
			}
			
			if (bestMatch != null)
			{
				if (currentSpectrum.getIntensity(peakID.intValue()) >= param.minPeakIntToScore)
				{
					matchedStat[0] += 1; // add matched peak count
					matchedStat[1] += currentSpectrum.getIntensity(peakID.intValue());
				}
			
				if (bestMajorMatch != null) // only count major ion for ion coverage
				{
					if (bestMajorMatch.ion.fromSecondPeptide) // from peptideB
						majorMatchB.add(bestMajorMatch.ion.getMajorTag());
					else
						majorMatchA.add(bestMajorMatch.ion.getMajorTag());
				}
			}
		}
		
		int[] bIonCount = new int[2];
		
		// count number of b-ions
		for (Iterator<String> iter = majorMatchA.iterator(); iter.hasNext();)
			if (iter.next().charAt(1) == 'b')
				bIonCount[0]++;
		
		for (Iterator<String> iter = majorMatchB.iterator(); iter.hasNext();)
			if (iter.next().charAt(1) == 'b')
				bIonCount[1]++;
		
		// compute coverage percentage
		localScore[0] = matchedStat[0] / globalSpectrumStat[0]; // peak coverage
		localScore[1] = matchedStat[1] / globalSpectrumStat[1]; // intensity coverage
		
		localScore[2] = bIonCount[0] * 0.5 / (peptideLength[0] - 1); // ion coverage peptideA b-ion
		localScore[3] = (majorMatchA.size() - bIonCount[0]) * 0.5 / (peptideLength[0] - 1); // ion coverage peptideA y-ion
		localScore[4] = majorMatchA.size() * 0.5 / (peptideLength[0] - 1); // ion coverage peptideA
		
		localScore[5] = bIonCount[1] * 0.5 / (peptideLength[1] - 1); // ion coverage peptideB b-ion
		localScore[6] = (majorMatchB.size() - bIonCount[1]) * 0.5 / (peptideLength[1] - 1); // ion coverage peptideB y-ion	
		localScore[7] = majorMatchB.size() * 0.5 / (peptideLength[1] - 1); // ion coverage peptideB
		
		localScore[8] = Math.sqrt(localScore[4] * localScore[7]); // update ion coverage
		localScore[9] = Math.pow(localScore[0] * localScore[1] * localScore[8], 1.0 / 3); // GM score
		
		// HelperFunctions.debug("majorMatchA", majorMatchA);
		// HelperFunctions.debug("majorMatchB", majorMatchB);
		
		// HelperFunctions.debug("matched count", matchedStat);		
		// HelperFunctions.debug("matching score", localScore);
		
		return localScore;
	}
	
	// compute scores based solely on best match of each peak
	public void computeOverallScores()
	{		
		if (isCrosslink) // obtain boundary of crosslink sites
			crosslinkSiteBound = CrosslinkSiteIdentifier.crosslinkSiteBound(crosslink, param.crosslinker);
		
		double[] matchedStat = new double[2]; // count number of matched peak and total matched intensity
		
		HashSet<String> majorMatchA = new HashSet<String>();
		HashSet<String> majorMatchB = new HashSet<String>(); // record major ion matches for each peptide
		
		ArrayList<MatchedFragmentStruct> tempMatchList;
		MatchedFragmentStruct bestMatch, bestMajorMatch;
		Integer peakID;
		
		for (Iterator<Integer> iter = matchedMap.keySet().iterator(); iter.hasNext();)
		{
			peakID = iter.next();
			bestMatch = null;
			bestMajorMatch = null;
			tempMatchList = matchedMap.get(peakID);
			
			/*
			if (currentSpectrum.precursor.scanNumber == 1047)
			{
				HelperFunctions.debug("MSMSpeak", currentSpectrum.spectrum[peakID.intValue()]);
				HelperFunctions.debug("tempMatchList", tempMatchList);
			}
			*/
			
			if (!isCrosslink)
			{
				bestMatch = tempMatchList.get(0); // retrieve best match for this peak
				
				for (int i = 0; i < tempMatchList.size(); i++)
				{
					if (tempMatchList.get(i).isMajorIon()) // retrieve best major ion match
					{
						bestMajorMatch = tempMatchList.get(i);
						break;
					}
				}
			}
			
			else
			{
				for (int i = 0; i < tempMatchList.size(); i++) // look for the best match which agree with crosslink site interpretation
				{
					if (checkFragmentWithLinkSite(tempMatchList.get(i).ion, crosslinkSiteBound))
					{
						bestMatch = tempMatchList.get(i);
						break;  // only need the best one
					}
				}
				
				for (int i = 0; i < tempMatchList.size(); i++) // look for the best major ion match which agree with crosslink site interpretation
				{
					if (tempMatchList.get(i).isMajorIon() && checkFragmentWithLinkSite(tempMatchList.get(i).ion, crosslinkSiteBound))
					{
						bestMajorMatch = tempMatchList.get(i);
						break;  // only need the best one
					}
				}
			}
			
			if (bestMatch != null)
			{
				/*
				if (currentSpectrum.precursor.scanNumber == 1047)
				{
					HelperFunctions.debug("MSMSpeak", currentSpectrum.spectrum[peakID.intValue()]);
					HelperFunctions.debug("bestMatch", bestMatch);
				}
				*/
				
				if (currentSpectrum.getIntensity(peakID.intValue()) >= param.minPeakIntToScore)
				{
					matchedStat[0] += 1; // add matched peak count
					matchedStat[1] += currentSpectrum.getIntensity(peakID.intValue());
				}
				
				if (bestMajorMatch != null) // only count major ion for ion coverage
				{
					if (bestMajorMatch.ion.fromSecondPeptide) // from peptideB
						majorMatchB.add(bestMajorMatch.ion.getMajorTag());
					else
						majorMatchA.add(bestMajorMatch.ion.getMajorTag());
				}
			}
		}
		
		int[] bIonCount = new int[2];
		
		// count number of b-ions
		for (Iterator<String> iter = majorMatchA.iterator(); iter.hasNext();)
			if (iter.next().charAt(1) == 'b')
				bIonCount[0]++;
		
		// compute coverage percentage
		overallScores[0] = matchedStat[0] / globalSpectrumStat[0]; // peak coverage
		overallScores[1] = matchedStat[1] / globalSpectrumStat[1]; // intensity coverage
		
		overallScores[2] = bIonCount[0] * 0.5 / (peptideLength[0] - 1); // ion coverage peptideA b-ion
		overallScores[3] = (majorMatchA.size() - bIonCount[0]) * 0.5 / (peptideLength[0] - 1); // ion coverage peptideA y-ion
		overallScores[4] = majorMatchA.size() * 0.5 / (peptideLength[0] - 1); // ion coverage peptideA
		
		overallScores[8] = overallScores[4]; // ion coverage
		
		if (isCrosslink)
		{
			for (Iterator<String> iter = majorMatchB.iterator(); iter.hasNext();)
				if (iter.next().charAt(1) == 'b')
					bIonCount[1]++;
			
			overallScores[5] = bIonCount[1] * 0.5 / (peptideLength[1] - 1); // ion coverage peptideB b-ion
			overallScores[6] = (majorMatchB.size() - bIonCount[1]) * 0.5 / (peptideLength[1] - 1); // ion coverage peptideB y-ion	
			overallScores[7] = majorMatchB.size() * 0.5 / (peptideLength[1] - 1); // ion coverage peptideB
			
			overallScores[8] = Math.sqrt(overallScores[4] * overallScores[7]); // update ion coverage
		}
		
		overallScores[9] = Math.pow(overallScores[0] * overallScores[1] * overallScores[8], 1.0 / 3); // GM score
		
		/*
		if (currentSpectrum.precursor.scanNumber == 1047)
		{
			HelperFunctions.debug("majorMatchA", majorMatchA);
			HelperFunctions.debug("majorMatchB", majorMatchB);
		}
		*/
		
		// HelperFunctions.debug("source count", globalStat);
		// HelperFunctions.debug("matched count", matchedStat);		
		// HelperFunctions.debug("matching score", overallScores);
	}

	// output text result via synchronized method in 'HelperFunctions.appendToFile()'
	public void outputText()
	{
		String[] tempoutput;
		String[] output;
		
		if (!isCrosslink) // for linear peptide, just output
		{
			output = new String[1];
			output[0] = currentSpectrum.toReport() + "\t" + peptide.toReport() + "\tn/a\t" + HelperFunctions.dummyPeptideReport + "\tn/a\t" + precursorMassError;
			
			for (int i = 0; i < overallScores.length; i++)
				output[0] += "\t" + overallScores[i];
			
			HelperFunctions.appendToFile(".\\output\\" + param.outputFileName + "_allScores.out", output);
			// HelperFunctions.appendToFile(".\\temp\\scores-" + Thread.currentThread().getName() + ".temp", output);
		}

		else // for crosslink, output the global score followed by individual crosslink site possibility
		{
			String report = currentSpectrum.toReport();
			
			if (crosslink.isForward()) // only compute site-specific data for forward crosslinks
			{
				int[] siteA_A = CrosslinkSiteIdentifier.crosslinkSites(crosslink.peptideA, 'A', param.crosslinker);
				int[] siteA_B = CrosslinkSiteIdentifier.crosslinkSites(crosslink.peptideA, 'B', param.crosslinker);
				int[] siteB_A = CrosslinkSiteIdentifier.crosslinkSites(crosslink.peptideB, 'A', param.crosslinker);
				int[] siteB_B = CrosslinkSiteIdentifier.crosslinkSites(crosslink.peptideB, 'B', param.crosslinker);
				
				int numPossibility = siteA_A.length * siteB_B.length + siteA_B.length * siteB_A.length; // number of different crosslink site pairs
				
				tempoutput = new String[numPossibility + 1];
				tempoutput[0] = report + "\t" + crosslink.toReport() + "\t" + precursorMassError;
				
				for (int i = 0; i < overallScores.length; i++)
					tempoutput[0] += "\t" + overallScores[i];
				
				int current = 1;
				
				double[] localScore = new double[overallScores.length];
				
				for (int i = 0; i < siteA_A.length; i++) // A_A to B_B
				for (int j = 0; j < siteB_B.length; j++)
				{
					localScore = computeScore(siteA_A[i], siteB_B[j]); // compute GM score
				
					if (localScore[localScore.length - 1] >= param.minScoreSite)
					{
						tempoutput[current] = report + "\t" + crosslink.toReport(siteA_A[i], siteB_B[j]) + "\t" + precursorMassError;
				
						for (int k = 0; k < overallScores.length; k++)
							tempoutput[current] += "\t" + localScore[k];
				
						current++;
					}
				}
			
				for (int i = 0; i < siteA_B.length; i++) // A_B to B_A
				for (int j = 0; j < siteB_A.length; j++)
				{
					localScore = computeScore(siteA_B[i], siteB_A[j]);
					
					if (localScore[localScore.length - 1] >= param.minScoreSite)
					{
						tempoutput[current] = report + "\t" + crosslink.toReport(siteA_B[i], siteB_A[j]) + "\t" + precursorMassError;
					
					
						for (int k = 0; k < overallScores.length; k++)
							tempoutput[current] += "\t" + localScore[k];
					
						current++;
					}
				}
				
				output = new String[current];
				
				for (int i = 0 ; i < current; i++)
					output[i] = tempoutput[i];
			}
			
			else
			{
				output = new String[1];
				output[0] = report + "\t" + crosslink.toReport() + "\t" + precursorMassError;
				
				for (int i = 0; i < overallScores.length; i++)
					output[0] += "\t" + overallScores[i];
			}
			
			HelperFunctions.appendToFile(".\\output\\" + param.outputFileName + "_allScores.out", output);
			// HelperFunctions.appendToFile(".\\temp\\scores-" + Thread.currentThread().getName() + ".temp", output);
		}
	}

	public void run()
	{
		FragmentIonStruct currentFragment = null;
		ArrayList<FragmentIonStruct> tempFragments = null;
		MatchedFragmentStruct tempMatch = null;
		ArrayList<MatchedFragmentStruct> tempMatchList = null;
		ArrayList<Integer> matchedID = null;
		ArrayList<ChargedPeakStruct> sortedChargedPeaks = currentSpectrum.getSortedChargedPeaks(param);
		ChargedPeakStruct tempChargedPeak;
		Double[] sortedMasses = null;
		double massTolerance = param.fragmentTolerance, massError;
		boolean found;
		int[] chargeRange;
		
		ArrayList<Integer[]> sortedMassedIDAll = new ArrayList<Integer[]>();
		ArrayList<Double[]> sortedMassesAll = new ArrayList<Double[]>();
		Double[] tempSortedMasses = null;
		Integer[] tempSortedID = null;
		
		// System.out.println(sortedChargedPeaks.size());
		
		// debug
		/*
		if (currentSpectrum.precursor.scanNumber == 6239 && isCrosslink)
		{
			HelperFunctions.debug("crosslink", crosslink);
			HelperFunctions.debug("currentFragment", fragmentedIons);
		}
		*/
		
		try
		{
			if (param.isHighResFragmentTolerance) // ppm, charge state doesn't matter
			{
				sortedMasses = HelperFunctions.getMassArray(sortedChargedPeaks);

				while (!fragmentedIons.empty()) // repeat until the stack is empty
				{
					currentFragment = fragmentedIons.pop();
					chargeRange = HelperFunctions.getChargeRange(currentFragment, currentSpectrum.precursor.chargeState);
					matchedID = MassMatcher.match(sortedMasses, currentFragment.mass + MassInfo.proton, massTolerance);
					found = false;

					if (matchedID.size() > 0)
					{
						for (int i = 0; i < matchedID.size(); i++) // add match result to global map
						{
							tempChargedPeak = sortedChargedPeaks.get(matchedID.get(i).intValue());

							if (tempChargedPeak.chargeState >= chargeRange[0] && tempChargedPeak.chargeState <= chargeRange[1]) // within allowable charge state range
							{
								massError = HelperFunctions.getppmError(currentFragment.mass + MassInfo.proton, tempChargedPeak.massWithOneCharge);
								tempMatch = new MatchedFragmentStruct(massError, currentFragment, tempChargedPeak.chargeState);
								
								if (matchedMap.containsKey(new Integer(tempChargedPeak.peakID)))
									matchedMap.get(new Integer(tempChargedPeak.peakID)).add(tempMatch);
								else
								{
									tempMatchList = new ArrayList<MatchedFragmentStruct>();
									tempMatchList.add(tempMatch);
									matchedMap.put(new Integer(tempChargedPeak.peakID), tempMatchList);
								}
																
								found = true;
							}
						}
					}

					// add all neutral losses and a-ion to applicable major fragments
					// only allow neutral losses to be considered for identified major ions
					// always allow precursor
					if ((found && currentFragment.isUnmodifiedABYIon()) || currentFragment.isPrecursor())
					{
						// a-ion
						if (currentFragment.canHaveAion(param))
							fragmentedIons.push(new FragmentIonStruct(currentFragment, "CO"));

						// others
						tempFragments = new ArrayList<FragmentIonStruct>();
						tempFragments.add(currentFragment); // starting point
						generateAllNeutralLosses(0, currentFragment, tempFragments);

						for (int i = 1; i < tempFragments.size(); i++) // skip the source
							fragmentedIons.push(tempFragments.get(i));
					}
				}
				
				// debug
				/*
				if (currentSpectrum.precursor.scanNumber == 6239 && isCrosslink)
				{
					HelperFunctions.debug("crosslink", crosslink);
					HelperFunctions.debug("currentSpectrum", currentSpectrum);
				}
				*/
			}

			else // Da mass tolerance, charge state matters
			{
				for (int z = 1; z <= currentSpectrum.precursor.chargeState; z++) // split masses according to charge state
					HelperFunctions.getMassArray(sortedChargedPeaks, z, sortedMassedIDAll, sortedMassesAll);

				// System.out.println(sortedMassedIDAll);
				// System.out.println(sortedMassesAll);

				while (!fragmentedIons.empty()) // repeat until the stack is empty
				{
					currentFragment = fragmentedIons.pop();
					// System.out.println(currentFragment);
					chargeRange = HelperFunctions.getChargeRange(currentFragment, currentSpectrum.precursor.chargeState);
					found = false;

					for (int z = chargeRange[0]; z <= chargeRange[1]; z++) // only consider possible charge states
					{
						tempSortedID = sortedMassedIDAll.get(z - 1);
						tempSortedMasses = sortedMassesAll.get(z - 1);
						// HelperFunctions.debug("tempSortedMasses", tempSortedMasses);
						// HelperFunctions.debug("tempSortedMasses", tempSortedID);
						massTolerance = param.fragmentTolerance * z * 1000000 / (currentFragment.mass + MassInfo.proton);
						matchedID = MassMatcher.match(tempSortedMasses, currentFragment.mass + MassInfo.proton, massTolerance);
						// System.out.println(matchedID.size());
						
						/*
						if (currentSpectrum.precursor.scanNumber == 1047 && currentFragment.getMajorTag().equals("Bb-2") && currentFragment.isMajorIon() && !currentFragment.isCrosslink)
						{
							HelperFunctions.debug("currentFragment", currentFragment);
							HelperFunctions.debug("charge: " + z);
							HelperFunctions.debug("masses", tempSortedMasses);
						}
						*/

						if (matchedID.size() > 0)
						{
							for (int i = 0; i < matchedID.size(); i++) // add match result to global map
							{
								tempChargedPeak = sortedChargedPeaks.get(tempSortedID[matchedID.get(i).intValue()].intValue());
								massError = HelperFunctions.getppmError(currentFragment.mass + MassInfo.proton, tempChargedPeak.massWithOneCharge);
								tempMatch = new MatchedFragmentStruct(massError, currentFragment, z);
								
								if (matchedMap.containsKey(new Integer(tempChargedPeak.peakID)))
									matchedMap.get(new Integer(tempChargedPeak.peakID)).add(tempMatch);
								else
								{
									tempMatchList = new ArrayList<MatchedFragmentStruct>();
									tempMatchList.add(tempMatch);
									matchedMap.put(new Integer(tempChargedPeak.peakID), tempMatchList);
								}
								
								found = true;
							}
						}
					}

					// add all neutral losses and a-ion to applicable major fragments
					// only allow neutral losses to be considered for identified major ions
					// always allow precursor
					if ((found && currentFragment.isUnmodifiedABYIon()) || currentFragment.isPrecursor())
					{
						// System.out.println(currentFragment);
						if (currentFragment.isMajorIon()) // add C13 peak for low-res major b or y ion
							fragmentedIons.push(new FragmentIonStruct(currentFragment, "C13"));

						// a-ion
						if (currentFragment.canHaveAion(param))
							fragmentedIons.push(new FragmentIonStruct(currentFragment, "CO"));

						// others
						tempFragments = new ArrayList<FragmentIonStruct>();
						tempFragments.add(currentFragment); // starting point
						generateAllNeutralLosses(0, currentFragment, tempFragments);

						for (int i = 1; i < tempFragments.size(); i++) // skip the source
						 	fragmentedIons.push(tempFragments.get(i));
					}
				}
			}
			
			// sort all matches for each peak, favor major ion, then sort by mass error
			for (Iterator<Integer> iter = matchedMap.keySet().iterator(); iter.hasNext();)
				Collections.sort(matchedMap.get(iter.next()), new MatchedFragmentComparator());
			
			// System.out.println(matchedMap);
			
			computeOverallScores(); // compute global scores
			
			if (overallScores[overallScores.length - 1] >= param.minScoreTargetList && isCrosslink && param.outputTargetList) // pass score threshold
				HelperFunctions.appendToFile(".\\temp\\" + param.outputFileName + "_targetList.temp", toTargetListEntry());
			
			if (overallScores[overallScores.length - 1] > 0 && !param.outputTargetList)
				outputText(); // output text file
		}

		catch (Exception e)
		{
			HelperFunctions.debug("SpectrumMatcher", HelperFunctions.getStackTrace(e));
		}
	}
}