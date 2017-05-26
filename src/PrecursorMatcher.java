import java.util.ArrayList;
// import java.util.concurrent.ExecutorService;
// import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledThreadPoolExecutor;

// match observed precursor masses to theoretical crosslinks and linear peptides
public class PrecursorMatcher
{
	Double[] observedMasses; // list of observed masses, sorted
	Integer[] indexedObservedMasses; // contain original ID of sorted observed masses
	public SynchronizedTreeMap globalPrecursorMatches = new SynchronizedTreeMap(); // mapping between observed mass ID to matched peptides

	public PrecursorMatcher()
	{}

	// import observed precursors
	public void importObservedMasses(ArrayList<PrecursorStruct> candidates)
	{
		observedMasses = new Double[candidates.size()];
		PrecursorStruct temp;

		// extract m/z and charge state
		for (int i = 0; i < candidates.size(); i++)
		{
			temp = candidates.get(i);
			observedMasses[i] = new Double(temp.chargeState * temp.precursorMZ - temp.chargeState * MassInfo.proton); // compute mass wihtout and charges
		}

		indexedObservedMasses = HelperFunctions.getIndexArray(observedMasses); // obtain index of the sorted observedMasses
		// HelperFunctions.debug("precursor masses", observedMasses);
	}

	// match precursor masses and add result to 'globalMatches'
	// multi-thread
	public void matchPrecursors(SynchronizedTreeMap[] globalPeptideMasses, ParamStruct param)
	{
		// ExecutorService executor = Executors.newFixedThreadPool(param.numCPU, new IDThreadFactory("thread"));
		ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(param.numCPU, new IDThreadFactory("thread"));
		Double[] peptideMasses;
		Integer[] indexedPeptideMasses; 

		for (int i = 0; i < param.numCPU; i++)
		{
			peptideMasses = ((ArrayList<Double>) globalPeptideMasses[0].get("thread-" + i)).toArray(new Double[0]);
			
			if (peptideMasses.length > 0)
			{
				indexedPeptideMasses = HelperFunctions.getIndexArray(peptideMasses); // sort and obtain index array
				executor.execute(new PrecursorMatcherThread(globalPrecursorMatches, i, false, peptideMasses, indexedPeptideMasses, observedMasses, indexedObservedMasses, param.precursorTolerance)); // linear peptide + dead-end
			}
			
			peptideMasses = ((ArrayList<Double>) globalPeptideMasses[1].get("thread-" + i)).toArray(new Double[0]);
			
			if (peptideMasses.length > 0)
			{
				indexedPeptideMasses = HelperFunctions.getIndexArray(peptideMasses); // sort and obtain index array
				executor.execute(new PrecursorMatcherThread(globalPrecursorMatches, i, true, peptideMasses, indexedPeptideMasses, observedMasses, indexedObservedMasses, param.precursorTolerance)); // crosslink + loop
			}
			
			while (executor.getQueue().size() > 2 * param.numCPU) {} // wait if the queue size is too high
		}

		executor.shutdown();
        while (!executor.isTerminated()) {} // wait
	}
}