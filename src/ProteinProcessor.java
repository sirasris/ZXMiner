import java.util.*;
import java.io.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledThreadPoolExecutor;

// process protein sequence --> in silico digestion and zero-length crosslinking
public class ProteinProcessor
{
	public ArrayList<ProteinStruct> proteins; // list of all proteins
	ArrayList<ThreadLimitStruct> queue; // queue for generating all peptides
	public SynchronizedTreeMap[] globalPeptideMasses = new SynchronizedTreeMap[2]; // contain peptide and crosslink mass lists for all threads

	public ProteinProcessor()
	{
		proteins = new ArrayList<ProteinStruct>(); // reset each time
		queue = new ArrayList<ThreadLimitStruct>();
		globalPeptideMasses[0] = new SynchronizedTreeMap(); // peptide, deadend
		globalPeptideMasses[1] = new SynchronizedTreeMap(); // crosslink, loop
	}
	
	// clear memory
	public void clearResource()
	{
		globalPeptideMasses[0].clear();
		globalPeptideMasses[1].clear();
		queue.clear();
	}

	// load proteins from fasta file, expected to be in '.\sequence\'
	// 'flag' indicate forward or decoy database
	public void loadProteins(String fname, boolean isForward)
	{
		String seq, head, st;
		String[] tempst;
		ProteinStruct protein;

		try
		{
			BufferedReader infile = new BufferedReader(new FileReader(".\\sequence\\" + fname));
			head = infile.readLine();

			while (head != null && head.length() > 0)
			{
				st = infile.readLine();
				seq = "";

				while(st!= null && st.length() > 0 && st.charAt(0) != '>') // in case where sequence span multiple lines
				{
					seq += st;
					st = infile.readLine();
				}

				tempst = head.split(" "); // split protein ID and description at the first space

				protein = new ProteinStruct(proteins.size(), tempst[0], head.substring(tempst[0].length(), head.length()), seq.toUpperCase(), isForward);
				// HelperFunctions.debug("protein", protein);
				// System.out.println(seq.toUpperCase());
				proteins.add(protein);
				head = st;
			}
			
			infile.close();
		}

		catch (Exception e)
		{ HelperFunctions.debug("ProteinProcessor::loadProteins", HelperFunctions.getStackTrace(e)); }
	}

	// digest all proteins
	public void digest(ProteaseStruct protease)
	{
		for (Iterator<ProteinStruct> iter = proteins.iterator(); iter.hasNext();)
		{
			ProteinStruct temp = iter.next();
			temp.computeCleaveSites(protease);
		}
	}

	// generate all peptides
	// multi-thread
	public void generatePeptides(ParamStruct param)
	{
		// ExecutorService executor = Executors.newFixedThreadPool(param.numCPU, new IDThreadFactory("thread"));
		ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(param.numCPU, new IDThreadFactory("thread"));
		int numStartSite, chunkSize;

		// add to queue
		// for each starting protein
		for (int i = 0; i < proteins.size(); i++)
		{
			numStartSite = proteins.get(i).getNumStartSites();
			chunkSize = numStartSite / param.numCPU;

			for (int j = 0; j < param.numCPU - 1; j++)
			{
				queue.add(new ThreadLimitStruct(i, chunkSize * j, chunkSize * (j + 1) - 1)); // sourceID2 = '-1' for linear peptides

				for (int k = i; k < proteins.size(); k++)
					queue.add(new ThreadLimitStruct(i, k, chunkSize * j, chunkSize * (j + 1) - 1)); // crosslinks
			}

			queue.add(new ThreadLimitStruct(i, chunkSize * (param.numCPU - 1), numStartSite - 1)); // sourceID2 = '-1' for linear peptides, for last chunk

			for (int k = i; k < proteins.size(); k++)
				queue.add(new ThreadLimitStruct(i, k, chunkSize * (param.numCPU - 1), numStartSite - 1)); // crosslinks, for last chunk
		}

		// initialize global mass lists
		for (int i = 0; i < param.numCPU; i++)
		{
			globalPeptideMasses[0].put("thread-" + i, new ArrayList<Double>());
			globalPeptideMasses[1].put("thread-" + i, new ArrayList<Double>());
		}

		// HelperFunctions.debug("peptide queue", queue);
		for (Iterator<ThreadLimitStruct> iter = queue.iterator(); iter.hasNext();)
		{
			executor.execute(new ProteinProcessorThread(globalPeptideMasses, proteins, iter.next(), param));
			
			while (executor.getQueue().size() > 2 * param.numCPU) {} // wait if the queue size is too high
		}

		executor.shutdown();
        while (!executor.isTerminated()) {} // wait
	}
	
	// generate all peptides
	// multi-thread
	public void generatePeptidesAndMatch(PrecursorMatcher matcher, ParamStruct param)
	{
		ExecutorService executor = Executors.newFixedThreadPool(param.numCPU, new IDThreadFactory("thread"));
		int numStartSite, chunkSize;

		// add to queue
		// for each starting protein
		for (int i = 0; i < proteins.size(); i++)
		{
			numStartSite = proteins.get(i).getNumStartSites();
			chunkSize = numStartSite / param.numCPU;

			for (int j = 0; j < param.numCPU - 1; j++)
			{
				queue.add(new ThreadLimitStruct(i, chunkSize * j, chunkSize * (j + 1) - 1)); // sourceID2 = '-1' for linear peptides

				for (int k = i; k < proteins.size(); k++)
					queue.add(new ThreadLimitStruct(i, k, chunkSize * j, chunkSize * (j + 1) - 1)); // crosslinks
			}

			queue.add(new ThreadLimitStruct(i, chunkSize * (param.numCPU - 1), numStartSite - 1)); // sourceID2 = '-1' for linear peptides, for last chunk

			for (int k = i; k < proteins.size(); k++)
				queue.add(new ThreadLimitStruct(i, k, chunkSize * (param.numCPU - 1), numStartSite - 1)); // crosslinks, for last chunk
		}

		// initialize global mass lists
		for (int i = 0; i < param.numCPU; i++)
		{
			globalPeptideMasses[0].put("thread-" + i, new ArrayList<Double>());
			globalPeptideMasses[1].put("thread-" + i, new ArrayList<Double>());
		}

		// HelperFunctions.debug("peptide queue", queue);
		for (Iterator<ThreadLimitStruct> iter = queue.iterator(); iter.hasNext();)
			executor.execute(new ProcessAndMatchPrecursor(globalPeptideMasses, proteins, iter.next(), param, matcher.globalPrecursorMatches, matcher.observedMasses, matcher.indexedObservedMasses));

		executor.shutdown();
        while (!executor.isTerminated()) {} // wait
	}

	// return number of proteins
	public int numProtein()
	{ return proteins.size(); }
}