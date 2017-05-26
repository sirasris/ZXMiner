import java.util.*;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;

// top-level class which dictate the flow of the whole pipeline
public class XlinkMiner
{
	ArrayList<PrecursorStruct> precursorInfo;
	ArrayList<String> mzXMLFiles;

	public XlinkMiner(ParamStruct parameter)
	{
		long startTime = System.currentTimeMillis();
		long currentTime = startTime;

		ParamStruct param;

	// LOAD PARAMETERS
		if (parameter == null)
			param = new ParamStruct(); // default parameters
		else
			param = parameter; // load parameters
		
		MassInfo.clearMemory();
		MassInfo.addModifications(param);
		currentTime = HelperFunctions.reportStatus("SETUP PARAMETERS", currentTime);
		
		if (param.outputFileName.endsWith("_allScores.out"))
		{
			HelperFunctions.debug("Existing allScores.out file specified. Skip to filtering result.");
			param.autosetOutputFileName(); // format outputFileName for compatibility with CandidatePicker
			
			CandidatePicker.pickSpectra(param);
			currentTime = HelperFunctions.reportStatus("SUMMARIZE SPECTRA", currentTime);
			
			CandidatePicker.pickPeptide(param);
			currentTime = HelperFunctions.reportStatus("SUMMARIZE PEPTIDES", currentTime);
			
			CandidatePicker.pickHighConfPeptide(param);
			currentTime = HelperFunctions.reportStatus("SUMMARIZE CROSSLINKED SITES", currentTime);
		}
		
		else
		{
			// LOAD PRECURSOR INFORMATION
			readCandidates(param.candidateFilePath); // load precursor information
			currentTime = HelperFunctions.reportStatus("LOAD CANDIDATES", currentTime);

		// LOAD PROTEIN DATABASE
			ProteinProcessor proteinproc = new ProteinProcessor();
			proteinproc.loadProteins(param.forwardFastaFileName, true); // load forward database

			if (param.decoyFastaFileName.length() > 0)
				proteinproc.loadProteins(param.decoyFastaFileName, false); // optionally load decoy database
			
			currentTime = HelperFunctions.reportStatus("LOAD PROTEIN DATABASES", currentTime);

		// WRITE ONLY MATCHED PEPTIDES
		    proteinproc.digest(param.protease); // compute cleavage site
			PrecursorMatcher precursormatch = new PrecursorMatcher();
			precursormatch.importObservedMasses(precursorInfo); // process precursor masses
			
			proteinproc.digest(param.protease); // compute cleavage site
			proteinproc.generatePeptidesAndMatch(precursormatch, param); // generate peptides and match precursor mass at once
			currentTime = HelperFunctions.reportStatus("CROSSLINKING AND PRECURSOR MATCHING", currentTime);
		
		// SPECTRUM MATCHING
			proteinproc.clearResource(); // clear resource
			DataGrabber grabber = new DataGrabber(param.mzXMLFilePath); // setup path to mzXML files
			grabber.grabAll(precursormatch.globalPrecursorMatches, precursorInfo, mzXMLFiles, proteinproc.proteins, param); // grab data and perform spectral match
			currentTime = HelperFunctions.reportStatus("SPECTRUM MATCHING", currentTime);
			
			if (param.outputTargetList)
			{
				GenerateTargetList.generate(param);
				currentTime = HelperFunctions.reportStatus("GENERATE TARGET LIST", currentTime);
			}
			
			else
			{
				CandidatePicker.pickSpectra(param);
				currentTime = HelperFunctions.reportStatus("SUMMARIZE SPECTRA", currentTime);
				
				CandidatePicker.pickPeptide(param);
				currentTime = HelperFunctions.reportStatus("SUMMARIZE PEPTIDES", currentTime);
				
				CandidatePicker.pickHighConfPeptide(param);
				currentTime = HelperFunctions.reportStatus("SUMMARIZE CROSSLINKED SITES", currentTime);
			}
			
			clearTempFolder(param);
			currentTime = HelperFunctions.reportStatus("CLEANING TEMP FOLDER", currentTime);
			HelperFunctions.reportStatus("DONE", startTime);
		}

	// CODE FOR DEBUGGING BASE64 CONVERSION
		// String binstr = HelperFunctions.decodeBase64("AA");
		// System.out.println(binstr);
		// System.out.println(HelperFunctions.encodeBase64(binstr));

	// CODE FOR DEBUGGING RANDOMACCESS TO FILE
		/*
    	try
    	{
    		RandomAccessFile testRead = new RandomAccessFile(".\\mzXML\\example.mzXML", "r");
    		System.out.println(testRead.length());
    		long offset = 1772;
	    	testRead.seek(offset);
	    	System.out.println(testRead.readLine().trim());
	    	System.out.println(testRead.readLine().trim());
	    	System.out.println(testRead.readLine().trim());

	    	offset = 69154;
	    	testRead.seek(offset);
			System.out.println(testRead.readLine().trim());
			System.out.println(testRead.readLine().trim());
			System.out.println(testRead.readLine().trim());

			testRead.close();
    	}

    	catch (Exception e)
    	{ HelperFunctions.debug("RAF", HelperFunctions.getStackTrace(e)); }
    	*/

    // CODE FOR DEBUGGING INDEX EXTRACTION
    	/*
    	DataGrabber grabber = new DataGrabber();
    	ArrayList<String> test = new ArrayList<String>(1);
    	test.add("example");
    	test.add("QE_20130222_GST_0C_Xlinks_12");
		grabber.initmzXMLIndex(test, param);
		*/

	// CODE FOR DEBUGGING THREADPEPTIDE EXTRACTION
		/*
		DataGrabber grabber = new DataGrabber();
		try
		{
			RandomAccessFile reader = new RandomAccessFile(".\\temp\\peptide-thread-0.temp", "r");
			PeptideStruct peptide = grabber.grabPeptide(reader, 0, proteinproc.proteins, param);
			System.out.println(peptide);
		}

		catch (Exception e)
    	{ HelperFunctions.debug("RAF", HelperFunctions.getStackTrace(e)); }
    	*/
		
	// CODE FOR DEBUGGING CONCURRENCY OF PROTEIN PROCESSOR
		/*
		DataGrabber grabber = new DataGrabber();
		try
		{
			for (int i = 0; i < param.numCPU; i++)
			{
				RandomAccessFile reader = new RandomAccessFile(".\\temp\\peptide-thread-" + i + ".temp", "r");
				ArrayList<Double> masses = (ArrayList<Double>) proteinproc.globalPeptideMasses[0].get("thread-" + i);
				
				for (int j = 0; j < masses.size(); j++)
				{
					PeptideStruct peptide = grabber.grabPeptide(reader, j, proteinproc.proteins, param);
					
					if (Math.abs(MassInfo.getMass(peptide) - masses.get(j).doubleValue()) > 0.1)
						System.out.println(MassInfo.getMass(peptide) + ", " + masses.get(j));
				}
				
				reader.close();
				
				reader = new RandomAccessFile(".\\temp\\crosslink-thread-" + i + ".temp", "r");
				masses = (ArrayList<Double>) proteinproc.globalPeptideMasses[1].get("thread-" + i);
				
				for (int j = 0; j < masses.size(); j++)
				{
					CrosslinkStruct crosslink = grabber.grabCrosslink(reader, j, proteinproc.proteins, param);
					
					if (Math.abs(MassInfo.getMass(crosslink, param) - masses.get(j).doubleValue()) > 0.1)
						System.out.println(MassInfo.getMass(crosslink, param) + ", " + masses.get(j));
				}
				
				reader.close();
			}
		}

		catch (Exception e)
    	{ HelperFunctions.debug("RAF", HelperFunctions.getStackTrace(e)); }
    	*/
	}

	// read a candidate file and create ID table
	public void readCandidates(String fpath)
	{
		precursorInfo = new ArrayList<PrecursorStruct>();
		mzXMLFiles = new ArrayList<String>();
		String st;
		String[] tempst;
		int temp; // ID starts from 0

		try
		{
			BufferedReader infile = new BufferedReader(new FileReader(fpath));
			st = infile.readLine();

			// repeat until end of file
			while (st != null)
			{
				if (st.length() > 0)
				{
					tempst = st.split("\t");
					temp = mzXMLFiles.indexOf(tempst[0]); // find raw file ID

					if (temp > -1) // found
						precursorInfo.add(new PrecursorStruct(temp, tempst));
					else // add new raw file
					{
						precursorInfo.add(new PrecursorStruct(mzXMLFiles.size(), tempst));
						mzXMLFiles.add(tempst[0]);
					}
				}

				st = infile.readLine();
			}

			// HelperFunctions.debug("mzXML file list", mzXMLFiles);
			infile.close();
		}

		catch (Exception e)
		{ HelperFunctions.debug("XlinkMiner::readCandidates", HelperFunctions.getStackTrace(e)); }
	}

	// clear temporary files
	public void clearTempFolder(ParamStruct param)
	{
		try
		{
			for (int i = 0; i < param.numCPU; i++)
			{
				Files.deleteIfExists(Paths.get(".\\temp\\peptide-thread-" + i + ".temp"));
				Files.deleteIfExists(Paths.get(".\\temp\\crosslink-thread-" + i + ".temp"));
			}
			
			if (param.outputTargetList)
				Files.deleteIfExists(Paths.get(".\\temp\\" + param.outputFileName + "_targetList.temp"));
		}
		
		catch (Exception e)
		{ HelperFunctions.debug("XlinkMiner::clearTempFolder", HelperFunctions.getStackTrace(e)); }
	}
	
	/*
	// expect at least 3 inputs: [path to mzXML files folder, path to candidate file, and forward database] + (decoy database and/or parameter file)
	public static void main(String[] args)
	{
		new XlinkMiner(null);		
	}
	*/
}