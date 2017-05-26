import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.TreeMap;

// data structure for search parameters
public class ParamStruct
{
	// DIGESTION
	public final ProteaseStruct protease;
	public final int maxMissedCleave;
	public final int minPeptideLength;
	public final int maxPeptideLength;

	// CROSSLINKING
	public final CrosslinkerStruct crosslinker;
	public final boolean allowAdjacent;

	// MODIFICATION
	public final String fixedModification;
	public final String varModification;
	public final int maxVarModPerPeptide;

	// MASS TOLERANCE
	public final double precursorTolerance; // always ppm
	public final double fragmentTolerance;
	public final boolean isHighResFragmentTolerance; // flag for Da or ppm tolerance

	// MS/MS PREPROCESSING
	public final double minPeakInt; // minimum peak intensity to take into account when computing matches
	public final double minPeakIntToScore; // minimum peak intensity to take into account when computing IntCov and PeakCov scores
	public final double isotopeWindowppm;
	public final boolean performDeisotope;

	// MS/MS MATCHING
	public final int numAion;
	public final int numNeutralLoss; // number of neutral losses allowed each for H2O and NH3, doubled for precursor peak
	public final boolean enableMetOxLoss; // consider loss of (OH)S(CH3) from oxidized methionine

	// FILE NAMES AND PATHS
	public final String forwardFastaFileName;
	public final String decoyFastaFileName;
	public final String candidateFilePath;
	public final String mzXMLFilePath;
	public String outputFileName; // output file name
	
	// MISC
	public final int numCPU;
	public final int peptideStructBase64Length;
	public final int crosslinkStructBase64Length; // length of base64 code
	
	// OUTPUT
	public final double minScoreTargetList; // minimum GM score to output to target list 
	public final double targetListRTHalfWindow; // half the width of retention time window for the target list, in minutes
	public final boolean outputTargetList; // whether to output target list for high-res MS/MS
	public final double maximumFDR; // maximum FDR level to output detailed crosslink information, in decimal, e.g. 0.05 = 5% FDR
	public final double minScoreSite; // minimum GM score to output crosslinked site-level information
	public final double minIonCov; // minimum Ion Coverage score to output crosslink ID
	public final double minDeltaScore; // minimum dCn score to output crosslink ID

	// default settings
	public ParamStruct()
	{
		// DATABASES
		forwardFastaFileName = "GST.fasta";
		decoyFastaFileName = "GST_reverse.fasta";
		
		// FILES
		candidateFilePath = ".\\input\\GST Oct 2011 Sample3 targetedOBOB precursor 071513.txt";
		mzXMLFilePath = "C:\\Research\\Crosslinking Projects\\Standard Proteins\\GST\\October 2011\\Raw Files\\Targeted OBOB\\";
		outputFileName = "GST Oct2011 Sample3 targetedOBOB reanalysis 100913";
		
		// DIGESTION
		protease = new ProteaseStruct(false, "KR"); // trypsin w/o ?P restriction
		maxMissedCleave = 3;
		minPeptideLength = 7;
		maxPeptideLength = 50;
		fixedModification = "C,carboxyamidomethyl,57.02146";
		varModification = "M,oxidation,15.99491";
		maxVarModPerPeptide = 3;
		
		// CROSSLINKING
		crosslinker = new CrosslinkerStruct(false, "K[", "DE]", -18.01528); // EDC
		allowAdjacent = false;
		
		// MASS TOLERANCES
		precursorTolerance = 10; // 10 for 60K full scans
		fragmentTolerance = 15;
		isHighResFragmentTolerance = true;
		
		// PEAK INTENSITY
		minPeakInt = 1000;
		minPeakIntToScore = 5000;
		
		// DE-ISOTOPING
		isotopeWindowppm = 20;
		performDeisotope = isHighResFragmentTolerance;
		
		// NEUTRAL LOSSES
		numAion = 5;
		numNeutralLoss = 2;
		enableMetOxLoss = true;
		
		// # OF CPU CORES
		numCPU = 2;
		
		// OUTPUT
		minScoreSite = 0.0;
		minIonCov = 0.2;
		minDeltaScore = 0.1;
		maximumFDR = 0.05;
		
		// OUTPUT - TARGET LIST
		targetListRTHalfWindow = 4;
		minScoreTargetList = 0.5;
		outputTargetList = false;

		// AUTOMATICALLY SET
		peptideStructBase64Length = (32 * (3 + 2 * maxMissedCleave) + 2 * (maxMissedCleave % 3)) / 6;
		crosslinkStructBase64Length = peptideStructBase64Length * 2;
	}

	// load
	public ParamStruct(TreeMap<String, String> setupValue)
	{
		if (setupValue.get("paramFilePath") != null && new File(setupValue.get("paramFilePath")).exists())
		{
			// DATABASES
			forwardFastaFileName = setupValue.get("forwardFastaFileName");
			decoyFastaFileName = setupValue.get("decoyFastaFileName");
			
			// FILES
			candidateFilePath = setupValue.get("candidateFilePath");
			mzXMLFilePath = setupValue.get("mzXMLFilePath");
			outputFileName = setupValue.get("outputFileName");
			
			TreeMap<String, String> paramValue = new TreeMap<String, String>();
			
			try
			{
				BufferedReader infile = new BufferedReader(new FileReader(setupValue.get("paramFilePath")));
				String st = infile.readLine();
				String[] tempst;
				
				while (st != null)
				{
					if (st.trim().length() > 0 && !st.startsWith("//")) // ignore comments and empty line
					{
						st = st.split("//")[0].trim(); // remove comments and empty spaces
						tempst = st.split("=");
						paramValue.put(tempst[0].trim(), tempst[1].trim());
					}
					
					st = infile.readLine();
				}
				
				
				infile.close();
			}
			
			catch (Exception e)
			{ HelperFunctions.debug("ParamStruct", "Error loading parameters: " + HelperFunctions.getStackTrace(e)); }
			
			// DIGESTION
			
			if (paramValue.get("cleavageDirection").equals("C-term"))
			{
				if (paramValue.get("cleavageRestriction").equals("none"))
					protease = new ProteaseStruct(false, paramValue.get("cleavageSite"));
				else
					protease = new ProteaseStruct(false, paramValue.get("cleavageRestriction"), paramValue.get("cleavageSite"));
			}
			
			else
			{
				if (paramValue.get("cleavageRestriction").equals("none"))
					protease = new ProteaseStruct(true, paramValue.get("cleavageSite"));
				else
					protease = new ProteaseStruct(true, paramValue.get("cleavageRestriction"), paramValue.get("cleavageSite"));
			}
			
			maxMissedCleave = Integer.parseInt(paramValue.get("maxMissedCleave"));
			minPeptideLength = Integer.parseInt(paramValue.get("minPeptideLength"));
			maxPeptideLength = Integer.parseInt(paramValue.get("maxPeptideLength"));
			fixedModification = paramValue.get("fixedModification");
			varModification = paramValue.get("varModification");
			maxVarModPerPeptide = Integer.parseInt(paramValue.get("maxVarModPerPeptide"));
			
			// CROSSLINKING
			// System.out.println(paramValue.get("crosslinkerDeltaMass"));
			crosslinker = new CrosslinkerStruct(false, paramValue.get("crosslinkerSiteA"), paramValue.get("crosslinkerSiteB"), Double.parseDouble(paramValue.get("crosslinkerDeltaMass"))); // EDC
			allowAdjacent = paramValue.get("allowAdjacent").equals("true");
			
			// MASS TOLERANCES
			precursorTolerance = Double.parseDouble(paramValue.get("precursorTolerance")); // 10 for 60K full scans
			fragmentTolerance = Double.parseDouble(paramValue.get("fragmentTolerance"));
			isHighResFragmentTolerance = paramValue.get("isHighResFragmentTolerance").equals("true");
			
			// PEAK INTENSITY
			minPeakInt = Double.parseDouble(paramValue.get("minPeakInt"));
			minPeakIntToScore = Double.parseDouble(paramValue.get("minPeakIntToScore"));
			
			// DE-ISOTOPING
			isotopeWindowppm = Double.parseDouble(paramValue.get("isotopeWindowppm"));
			performDeisotope = isHighResFragmentTolerance;
			
			// NEUTRAL LOSSES
			numAion = Integer.parseInt(paramValue.get("numAion"));
			numNeutralLoss = Integer.parseInt(paramValue.get("numNeutralLoss"));
			enableMetOxLoss = paramValue.get("enableMetOxLoss").equals("true");
			
			// # OF CPU CORES
			numCPU = Integer.parseInt(setupValue.get("numCPU"));
			
			// OUTPUT
			minScoreSite = Double.parseDouble(paramValue.get("minScoreSite"));
			
			// OUTPUT - TARGET LIST
			outputTargetList = setupValue.get("outputTargetList").equals("true");
			
			if (outputTargetList)
			{
				targetListRTHalfWindow = Double.parseDouble(setupValue.get("targetListRTHalfWindow"));
				minScoreTargetList = Double.parseDouble(setupValue.get("minScoreTargetList"));
				
				minIonCov = 0;
				minDeltaScore = 0;
				maximumFDR = 1;
			}
			
			else
			{
				targetListRTHalfWindow = 0;
				minScoreTargetList = 0;
				
				minIonCov = Double.parseDouble(setupValue.get("minIonCov"));
				minDeltaScore = Double.parseDouble(setupValue.get("minDeltaScore"));
				maximumFDR = Double.parseDouble(setupValue.get("maximumFDR"));
			}
			
			// AUTOMATICALLY SET
			peptideStructBase64Length = (32 * (3 + 2 * maxMissedCleave) + 2 * (maxMissedCleave % 3)) / 6;
			crosslinkStructBase64Length = peptideStructBase64Length * 2;
		}
		
		else // this must be a re-processing attempt
		{
			// DATABASES
			forwardFastaFileName = "";
			decoyFastaFileName = "";
			
			// FILES
			candidateFilePath = "";
			mzXMLFilePath = "";
			outputFileName = setupValue.get("outputFileName");
			
			// DIGESTION
			protease = new ProteaseStruct(false, "KR"); // trypsin w/o ?P restriction
			maxMissedCleave = 3;
			minPeptideLength = 7;
			maxPeptideLength = 50;
			fixedModification = "C,carboxyamidomethyl,57.02146";
			varModification = "M,oxidation,15.99491";
			maxVarModPerPeptide = 3;
			
			// CROSSLINKING
			crosslinker = new CrosslinkerStruct(false, "K[", "DE]", -18.01528); // EDC
			allowAdjacent = false;
			
			// MASS TOLERANCES
			precursorTolerance = 10; // 10 for 60K full scans
			fragmentTolerance = 15;
			isHighResFragmentTolerance = true;
			
			// PEAK INTENSITY
			minPeakInt = 1000;
			minPeakIntToScore = 5000;
			
			// DE-ISOTOPING
			isotopeWindowppm = 20;
			performDeisotope = isHighResFragmentTolerance;
			
			// NEUTRAL LOSSES
			numAion = 5;
			numNeutralLoss = 2;
			enableMetOxLoss = true;
			
			// # OF CPU CORES
			numCPU = 2;
			
			// OUTPUT
			minScoreSite = 0.0;
			minIonCov = Double.parseDouble(setupValue.get("minIonCov"));
			minDeltaScore = Double.parseDouble(setupValue.get("minDeltaScore"));
			maximumFDR = Double.parseDouble(setupValue.get("maximumFDR")); // the only new input are these numbers
			
			// OUTPUT - TARGET LIST
			targetListRTHalfWindow = 4;
			minScoreTargetList = 0.5;
			outputTargetList = false;

			// AUTOMATICALLY SET
			peptideStructBase64Length = (32 * (3 + 2 * maxMissedCleave) + 2 * (maxMissedCleave % 3)) / 6;
			crosslinkStructBase64Length = peptideStructBase64Length * 2;
		}
	}

	// return whether the crosslinker can generate dead-ends
	public boolean hasDeadEnd()
	{ return crosslinker.hasDeadEnd; }
	
	// for filtering existing allScores.out file
	public void autosetOutputFileName()
	{
		String[] temp = outputFileName.split("\\\\");
		outputFileName = temp[temp.length - 1].split("_allScores")[0];
	}
}