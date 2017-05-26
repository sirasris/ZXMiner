import java.util.ArrayList;

// multi-threading class for PrecursorMatcher
public class PrecursorMatcherThread implements Runnable
{
	SynchronizedTreeMap globalMatches; // global precursor match list
	int peptideThreadID;
	boolean isCrosslink;
	double massTolerance; // mass tolerance
	Double[] peptideMasses; // list of peptide masses, sorted
	Integer[] indexedPeptideMasses; // map to original peptide ID
	Double[] observedMasses; // list of observed masses, sorted
	Integer[] indexedObservedMasses; // map to original precursor ID

	public PrecursorMatcherThread(SynchronizedTreeMap globalMatches, int peptideThreadID, boolean isCrosslink, Double[] peptideMasses, Integer[] indexedPeptideMasses, Double[] observedMasses, Integer[] indexedObservedMasses, double massTolerance)
	{
		this.globalMatches = globalMatches;
		this.peptideThreadID = peptideThreadID;
		this.isCrosslink = isCrosslink;
		this.peptideMasses = peptideMasses;
		this.indexedPeptideMasses = indexedPeptideMasses;
		this.observedMasses = observedMasses;
		this.indexedObservedMasses = indexedObservedMasses;
		this.massTolerance = massTolerance;
	}

	public void run()
	{
		try
		{
			ArrayList<Integer> tempmatches;
			ArrayList<ThreadPeptideStruct> matches;

			for (int i = 0; i < observedMasses.length; i++)
			{
				tempmatches = MassMatcher.match(peptideMasses, observedMasses[i].doubleValue(), massTolerance);

				if (tempmatches.size() > 0) // found some matches
				{
					matches = new ArrayList<ThreadPeptideStruct>();

					for (int j = 0; j < tempmatches.size(); j++) // update index
						matches.add(new ThreadPeptideStruct(peptideThreadID, isCrosslink, indexedPeptideMasses[tempmatches.get(j).intValue()]));
					
					globalMatches.addToValueList(indexedObservedMasses[i], matches);
				}
			}
		}

		catch (Exception e)
		{
			HelperFunctions.debug("PrecursorMatcherThread", HelperFunctions.getStackTrace(e));
			HelperFunctions.debug("peptideMasses", peptideMasses.length);
			HelperFunctions.debug("indexedPeptideMasses", indexedPeptideMasses.length);
		}
	}
}