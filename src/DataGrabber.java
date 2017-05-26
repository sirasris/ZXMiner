import java.util.ArrayList;
import java.util.TreeMap;
import java.util.Iterator;
import java.util.Stack;
import java.io.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.regex.Pattern;

// extract data from mzXML and peptide files and feed to SpectrumMatcher and FragmentGenerator
public class DataGrabber
{
	public SynchronizedTreeMap globalSpectralMatch;
	TreeMap<Integer, TreeMap<Integer, Long>> indexedSpectraOffset;
	RandomAccessFile[] mzXMLReaders, peptideReaders, crosslinkReaders;

	String mzXMLPath;
	final Pattern indexPattern = Pattern.compile("\\D*(\\d+)\\D*(\\d+)\\D*");

	public DataGrabber()
	{}
	
	public DataGrabber(String mzXMLPath)
	{
		this.mzXMLPath = mzXMLPath;

		if (!mzXMLPath.endsWith("\\")) // make sure the path ends with '\' symbol
			this.mzXMLPath += "\\";
	}

	// initialize spectrum-matching output file with header line
	public void initSpectrumMatchOutputFile(ParamStruct param)
	{
		try
		{
			BufferedWriter outfile = new BufferedWriter(new FileWriter(".\\output\\" + param.outputFileName + "_allScores.out"));
			String[] header = HelperFunctions.tempScoreFileHeader;
			
			for (int i = 0; i < header.length - 1; i++)
				outfile.write(header[i] + "\t");
			
			outfile.write(header[header.length - 1] + "\n");

			// outfile.flush();
			outfile.close();
		}

		catch (Exception e)
		{ HelperFunctions.debug("DataGrabber::initSpectrumMatchOutput", HelperFunctions.getStackTrace(e)); }
	}
	
	/*
	// initialize spectrum-matching output file with header line
	public void initTargetOutputFile(ParamStruct param)
	{
		try
		{
			BufferedWriter outfile = new BufferedWriter(new FileWriter(".\\output\\" + param.outputFileName + "_targetList.out"));
			outfile.close();
		}

		catch (Exception e)
		{ HelperFunctions.debug("DataGrabber::initTargetOutputFile", HelperFunctions.getStackTrace(e)); }
	}
	*/

	// obtain index of all scan number
	// multi-thread
	public void initmzXMLIndex(ArrayList<String> mzXMLFiles, ParamStruct param)
	{
		indexedSpectraOffset = new TreeMap<Integer, TreeMap<Integer, Long>>();
		ArrayList<Future<TreeMap<Integer, Long>>> future = new ArrayList<Future<TreeMap<Integer, Long>>>(mzXMLFiles.size());
		ExecutorService executor = Executors.newFixedThreadPool(param.numCPU, new IDThreadFactory("thread"));

		for (int i = 0; i < mzXMLFiles.size(); i++)
			future.add((Future<TreeMap<Integer, Long>>) executor.submit(new IndexGrabberThread(i, mzXMLPath + mzXMLFiles.get(i), indexPattern))); // reuse pattern

		executor.shutdown();
        while (!executor.isTerminated()) {} // wait

        for (int i = 0; i < mzXMLFiles.size(); i++)
        {
        	try
        	{
        		indexedSpectraOffset.put(new Integer(i), future.get(i).get()); // add to the true list
        	}

        	catch (Exception e)
        	{ HelperFunctions.debug("DataGrabber::initmzXMLIndex", HelperFunctions.getStackTrace(e)); }
        }

		// HelperFunctions.debug("mzXML spectra index", indexedSpectraOffset);
	}

	// initialize access to all mzXML files
	public void initmzXMLReaders(ArrayList<String> mzXMLFiles)
	{
		mzXMLReaders = new RandomAccessFile[mzXMLFiles.size()];

		try
		{
			for (int i = 0; i < mzXMLFiles.size(); i++)
				mzXMLReaders[i] = new RandomAccessFile(mzXMLPath + mzXMLFiles.get(i) + ".mzXML", "r");
		}

		catch (Exception e)
		{ HelperFunctions.debug("DataGrabber::initmzXMLReader", HelperFunctions.getStackTrace(e)); }
	}

	// initialize access to all peptide files
	public void initPeptideReaders(ParamStruct param)
	{
		peptideReaders = new RandomAccessFile[param.numCPU];
		crosslinkReaders = new RandomAccessFile[param.numCPU];

		try
		{
			for (int i = 0; i < param.numCPU; i++)
			{
				peptideReaders[i] = new RandomAccessFile(".\\temp\\peptide-thread-" + i + ".temp", "r");
				crosslinkReaders[i] = new RandomAccessFile(".\\temp\\crosslink-thread-" + i + ".temp", "r");
			}
		}

		catch (Exception e)
		{ HelperFunctions.debug("DataGrabber::initPeptideReaders", HelperFunctions.getStackTrace(e)); }
	}
	
	// close all input streams
	public void closeAllFileReaders()
	{
		try
		{
			for (int i = 0; i < mzXMLReaders.length; i++)
				mzXMLReaders[i].close();
			
			for (int i = 0; i < peptideReaders.length; i++)
			{
				peptideReaders[i].close();
				crosslinkReaders[i].close();
			}
		}

		catch (Exception e)
		{ HelperFunctions.debug("DataGrabber::initmzXMLReader", HelperFunctions.getStackTrace(e)); }
	}

	// grab peptide information
	public PeptideStruct grabPeptide(RandomAccessFile reader, int peptideID, ArrayList<ProteinStruct> proteins, ParamStruct param)
	{
		try
		{
			reader.seek((long) peptideID * param.peptideStructBase64Length);
			byte[] temp = new byte[param.peptideStructBase64Length];
			reader.read(temp);
			// HelperFunctions.debug("byte array", temp);
			// System.out.println(new String(temp));
			return new PeptideStruct(proteins, Base64Parser.base64toInteger(new String(temp)));
		}

		catch (Exception e)
		{ HelperFunctions.debug("DataGrabber::grabPeptide", HelperFunctions.getStackTrace(e)); }

		return null; // fails
	}

	// grab crosslink information
	public CrosslinkStruct grabCrosslink(RandomAccessFile reader, int crosslinkID, ArrayList<ProteinStruct> proteins, ParamStruct param)
	{
		try
		{
			reader.seek((long) crosslinkID * param.crosslinkStructBase64Length);
			byte[] tempA = new byte[param.peptideStructBase64Length];
			reader.read(tempA);

			reader.seek((long) crosslinkID * param.crosslinkStructBase64Length + param.peptideStructBase64Length);

			byte[] tempB = new byte[param.peptideStructBase64Length];
			reader.read(tempB);

			return new CrosslinkStruct(proteins, Base64Parser.base64toInteger(new String(tempA)), Base64Parser.base64toInteger(new String(tempB)));
		}

		catch (Exception e)
		{ HelperFunctions.debug("DataGrabber::grabPeptide", HelperFunctions.getStackTrace(e)); }

		return null; // fails
	}

	// grab spectrum information
	public SpectrumStruct grabSpectrum(RandomAccessFile reader, PrecursorStruct precursor, String rawfileName)
	{
		try
		{
			long offset = indexedSpectraOffset.get(new Integer(precursor.rawFileID)).get(new Integer(precursor.scanNumber)).longValue();
			reader.seek(offset);
			String st = reader.readLine().trim();
			String[] tempst;

			int parentScanNumber = -1, precision = 32;
			double retentionTime = -1, basePeakIntensity = -1, precursorIntensity = -1;
			boolean centroided = true;

			while (!st.startsWith("<precursorMz"))
			{
				if (st.startsWith("centroided"))
				{
					st = st.split("=")[1];
					st = st.substring(1, st.length() - 1); // remove "..."
					
					if (!st.equals("1"))
						centroided = false;
				}
				
				if (st.startsWith("retentionTime"))
				{
					st = st.split("=")[1];
					st = st.substring(3, st.length() - 2); // remove "PT...S"
					retentionTime = Double.valueOf(st) / 60.0; // convert to minutes
				}

				if (st.startsWith("basePeakIntensity"))
				{
					st = st.split("=")[1];
					st = st.substring(1, st.length() - 1); // remove "..."
					basePeakIntensity = Double.valueOf(st); // convert to minutes
				}

				st = reader.readLine().trim();
			}

			// reached precursorMz
			tempst = st.split(" ");

			for (int i = 0; i < tempst.length; i++)
			{
				if (tempst[i].startsWith("precursorScanNum"))
				{
					st = tempst[i].split("=")[1];
					st = st.substring(1, st.length() - 1); // remove "..."
					parentScanNumber = Integer.valueOf(st);
				}

				if (tempst[i].startsWith("precursorIntensity"))
				{
					st = tempst[i].split("=")[1];
					st = st.substring(1, st.length() - 1); // remove "..."
					precursorIntensity = Double.valueOf(st);
				}
			}

			// last entry are peaks
			st = reader.readLine().trim();

			while (!st.startsWith("contentType"))
			{
				if (st.startsWith("precision"))
				{
					st = st.split("=")[1];
					st = st.substring(1, st.length() - 1); // remove "..."
					precision = Integer.valueOf(st);
				}

				st = reader.readLine().trim();
			}

			// reached 'contentType=...'
			double[][] spectrum = Base64Parser.base64toDouble(st.substring(st.indexOf(">") + 1, st.length() - 8), precision);
			
			if (!centroided)
				spectrum = PeakCentroider.centroided(spectrum); // centroided data
			
			/*
			if (precursor.scanNumber == 13844)
			{
				for (int i = 0 ; i < spectrum.length; i++)
					System.out.println(spectrum[i][0] + "\t" + spectrum[i][1]);
			}
			*/

			return new SpectrumStruct(rawfileName, precursor, parentScanNumber, retentionTime, basePeakIntensity, precursorIntensity, spectrum);
		}

		catch (Exception e)
		{ HelperFunctions.debug("DataGrabber::grabSpectrum", HelperFunctions.getStackTrace(e)); }

		return null; // fails
	}

	// grab all matched pairs
	// multi-thread
	public void grabAll(SynchronizedTreeMap precursorMatches, ArrayList<PrecursorStruct> precursorInfo, ArrayList<String> mzXMLFiles, ArrayList<ProteinStruct> proteins, ParamStruct param)
	{
		globalSpectralMatch = new SynchronizedTreeMap();
		initmzXMLIndex(mzXMLFiles, param);

		initmzXMLReaders(mzXMLFiles); // prepare RandomAccessFile instances
		initPeptideReaders(param);
		
		if (!param.outputTargetList)
			initSpectrumMatchOutputFile(param); // initialize score file 
		
		/*
		if (param.outputTargetList) // initialize target list file
			initTargetOutputFile(param);
		*/
		
		// ExecutorService executor = Executors.newFixedThreadPool(param.numCPU, new IDThreadFactory("thread"));
		ScheduledThreadPoolExecutor executor = new ScheduledThreadPoolExecutor(param.numCPU, new IDThreadFactory("thread"));
		ArrayList<ThreadPeptideStruct> matches;
		Stack<FragmentIonStruct> fragmentedIons;
		Integer precursorID;
		ThreadPeptideStruct temppeptide;
		PrecursorStruct currentPrecursor;
		SpectrumStruct currentSpectrum;
		PeptideStruct currentPeptide;
		CrosslinkStruct currentCrosslink;

		for (Iterator iter = precursorMatches.keySet().iterator(); iter.hasNext();) // for each precursor that matched
		{
			precursorID = (Integer) iter.next();
			currentPrecursor = precursorInfo.get(precursorID.intValue());
			matches = precursorMatches.get(precursorID);

			currentSpectrum = grabSpectrum(mzXMLReaders[currentPrecursor.rawFileID], currentPrecursor, mzXMLFiles.get(currentPrecursor.rawFileID));

			if (param.performDeisotope && param.isHighResFragmentTolerance) // de-isotope spectrum
				currentSpectrum.setSpectrum(Deisotoper_ChiSquare.deisotope(currentSpectrum.spectrum, currentPrecursor.chargeState, param));

			// HelperFunctions.debug("spectrum", currentSpectrum);

			for (Iterator<ThreadPeptideStruct> iter2 = matches.iterator(); iter2.hasNext();) // for each matched peptide
			{
				temppeptide = iter2.next();

				if (temppeptide.isCrosslink) // crosslink
				{
					currentCrosslink = grabCrosslink(crosslinkReaders[temppeptide.threadID], temppeptide.peptideID, proteins, param);
					// fragmentedIons = FragmentGenerator.fragment(currentCrosslink, param);
					fragmentedIons = FragmentGenerator.fragmentManual(currentCrosslink, param);
					executor.execute(new SpectrumMatcher(fragmentedIons, currentSpectrum, currentCrosslink, param));
				}

				else // linear
				{
					currentPeptide = grabPeptide(peptideReaders[temppeptide.threadID], temppeptide.peptideID, proteins, param);
					// fragmentedIons = FragmentGenerator.fragment(currentPeptide);
					fragmentedIons = FragmentGenerator.fragmentManual(currentPeptide);
					executor.execute(new SpectrumMatcher(fragmentedIons, currentSpectrum, currentPeptide, param));
				}
				
				while (executor.getQueue().size() > 2 * param.numCPU) {} // wait if the queue size is too high
			}
		}

		executor.shutdown();
        while (!executor.isTerminated()) {} // wait
        
        closeAllFileReaders(); // close all input streams
	}
}