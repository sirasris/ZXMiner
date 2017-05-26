# ZXMiner
Source code for ZXMiner, a zero-length crosslinks identification tool originally published in http://pubs.acs.org/doi/abs/10.1021/pr400953w.

###############################################################

ZXMiner - Software tool for identification of zero-length chemical crosslinks

Sira Sriswasdi
PI: David W. Speicher, Ph.D.
The Wistar Institute & University of Pennsylvania, PA, USA

###############################################################

0. INTRODUCTION

1. REQUIREMENT

2. USER GUIDE
	2.1 installation
	2.2 preparing input files
	2.3 configuration
	2.4 types of output
	2.5 reprocessing of allscores.out
	2.6 running multiple instance of ZXMiner

3. SOURCE CODE & CITATION

###############################################################

0. INTRODUCTION
ZXMiner is a software tool for facilitaing identification of zero-length crosslinks.
The software is developed primarily based on crosslinking dataset using 
1-Ethyl-3-[3-dimethylaminopropyl]carbodiimide hydrochloride (EDC) and Thermo's LTQ Orbitrap XL 
mass spectrometer.

Because all of ZXMiner's scoring criteria do not involve any parameters that could be 
specific to the usage of EDC, ZXMiner should be able to handle data from other zero-length 
crosslinker just as well.

However, because different types of mass spectrometer do produce considerable variations in the 
resulting MS/MS spectra of the same crosslinked peptide, please be careful when using ZXMiner to analyze 
data from other instruments.

1. REQUIREMENT
Java 7. ZXMiner was developed using Java 7 without any external libraries.

2. USER GUIDE
2.1 installation
ZXMiner comes as a combination of ZXMiner.jar file and 5 main folders: config, input,
output, sequence, and temp. "config", "input", and "sequence" folders are provided just for the sake of
being organized. However, "output" folder is where all output will be generated and "temp" folder is 
where all temporary data files will be dynamically written and deleted. Both of the latters are essential.

To be able to see error messages from ZXMiner, please either utilize the 
accompanying "execute.bat" or run ZXMiner.jar through the command line. We plan to direct all error
messages to the GUI in future versions.

2.2 preparing input files
ZXMiner requires 3 types of inputs: a tab-delimited text files containing candidate precursors, 
a folder of mzXML files, and protein sequence files in fasta format. Input files should be put into the
"input" folder but this is not necessary. The same holds for protein sequence files. A few examples are 
provided in these folders. 

Basically, the tab-delimited text file containing candidate precursors will
have 4 columns: mzXML file name (w/o .mzXML extension), scan number, charge state, and precursor m/z. 
These information can be obtained from the MS/MS header information in mzXML files or from any 
label-free analysis software.

mzXML files can be obtained by converting MS raw file using publicly available tools such as
proteoWizard [http://proteowizard.sourceforge.net/]. Be careful to uncheck the "using zlib compression"
option.

2.3 configuration
Examples of configuration file are provided in the "config" folder. Most of the parameters 
are given names whose meaning can be recognized by eye and comments following each "//" symbol provide 
further explanations. The following is from "default_OBOB_config.txt":

// DIGESTION
cleavageDirection = C-term
cleavageSite = KR
cleavageRestriction = none				                  // default = trypsin without Proline restriction

maxMissedCleave = 3

minPeptideLength = 7
maxPeptideLength = 50 					                    // limit on peptide size (individual peptide)

fixedModification = C,carboxyamidomethyl,57.02146
varModification = M,oxidation,15.99491			        // format is "residue1,name1,mass1;residue1,name2,mass2;..."
maxVarModPerPeptide = 3
		
// CROSSLINKING
crosslinkerSiteA = K[
crosslinkerSiteB = DE]					                    // [ and ] represents the protein N- and C-terminal residues
crosslinkerDeltaMass = -18.01528			              // default = EDC

allowAdjacent = false					                      // whether to consider crosslinks of adjacent peptides
		
// MASS TOLERANCES
precursorTolerance = 10					                    // ppm unit, 10 for 60K full scans, 15 for 15K
fragmentTolerance = 15 					                    // ppm unit if isHighResFragmentTolerance = true, Da unit otherwise
isHighResFragmentTolerance = true
		
// PEAK INTENSITY
minPeakInt = 1000
minPeakIntToScore = 5000
		
// DE-ISOTOPING
isotopeWindowppm = 20 					                    // ppm unit, for determining isotopic envelopes
		
// NEUTRAL LOSSES
numAion = 5 						                            // number of a-ion to consider
numNeutralLoss = 2 					                        // number of water and ammonia loss to consider (each)
enableMetOxLoss = true 					                    // whether to consider loss of oxidized methionine from precursor ions

// OUTPUT
minScoreSite = 0.2 					                        // for reducing the output file size, only retain information for candidates with GM score at least minScoreSite 


The first series of parameters - cleavageDirection, cleavageSite, and cleavageRestriction - concern the digestion enzyme, e.g. trypsin.
"cleavageDirection" defines the direction of cleavage relative to the cleavage site.
"cleavageSite" defines the residues which are the cleavage site.
"cleavageRestriction" defines residues that can block cleavage.

Next are for generating peptides.
"maxMissedCleave" defines the maxmimum number of missed cleavage in each peptide generated.
"minPeptideLength" and "maxPeptideLength" defines the sizes of peptide to be generated.

Next are amino acid modicfications.
"fixedModification" defines static modifications. Example: C,carboxyamidomethyl,57.02146
"varModification" defines variable modifications. Example: M,oxidation,15.99491
"maxVarModPerPeptide" defines maximum number of variable modifications on each peptide.

The format for specifying modification is "residue1,name1,mass1" separated by ";", i.e. "residue1,name1,mass1;residue2,name2,mass2"
		
The parameters for crosslinker are as followed:
"crosslinkerSiteA" defines residues to be used as the first crosslinkable site.
"crosslinkerSiteB" defines residues to be used as the second crosslinkable site.
"crosslinkerDeltaMass" defines the mass lost during the crosslinking reaction, i.e. -18.01528 for the loss of water molecule.

Use "[" and "]" to designate the protein N- and C-terminal.

"allowAdjacent" defines whether cross-links involving adjacent peptides on the same protein to be output. This is recommended to be
turned off for EDC since the masses of these cross-links and the alternative explanation of a long contiguous peptides are exactly
the same and there will be only a few fragmented ions specific to the cross-linked peptides.
		
Next are for mass tolerance.
"precursorTolerance" defines mass tolerance for precursor m/z (in ppm). Please adjust this according to the resolution of full MS data.
"fragmentTolerance" defines mass tolerance for MS/MS m/z (in amu or ppm, as indicated by the next parameter).
"isHighResFragmentTolerance" indicates whether the mass tolerance provided for "fragmentTolerance" is in amu or ppm units.
		
Next are for filtering MS/MS peak intensities.
"minPeakInt" defines minimal MS/MS peak intensity to retain for analysis. All peaks with less intensities will be removed prior to all MS/MS processing steps.
"minPeakIntToScore" defines minimal MS/MS peak intensity to be included in the scoring of peak coverage and intensity coverage for each spectrum.

Peaks with intensities between these two settings will be present in the MS/MS processing step for the purpose of identifying isotopic envelopes but will 
not affect the scores.
		
"isotopeWindowppm" defines the ppm mass tolerance for adjacent peaks in the same isotopic envelope. This is used to determine preliminary isotopic envelope
candidates which will then be confirmed using the profile of peak intensities.
		
Next are for neutral losses.
"numAion" defines the number of a-ion to consider, e.g. use 5 for considering a-1, a-2, ..., upto a-5 ions.
"numNeutralLoss" defines the maxmimal number of neutral losses allowed on each relevant fragmented ions. This limit is for water and ammonia individually,
e.g. use 2 for allowing upto 2 water losses and 2 ammonia losses (for a maximum of 4 losses on some ions).
"enableMetOxLoss" defines whether the loss of oxidized methionine side chain (CH3 SOH, about 64 Dalton) from the precursor ions will be considered.

The last parameter is for controlling the size of "allScores.out" file.
"minScoreSite" defines the mimimal GM scores of candidate cross-links that will be output to the allScores.out file. This setting is a balance between
file size and reprocessibility. Setting this to zero can result in 10+ GB file size but will allow complete reprocessing without having to re-run 
ZXMiner again - unless there are other drastic changes in analysis parameters.

2.4 types of output
To facilitate follow-up targeted MS analysis, ZXMiner can output a list of target m/z and retention time to be incorporated into the 
setting on the mass spectrometer instrument. To enable this mode, please check the "Output target list only" option and specify retention window and minimal GM score.

Otherwise, ZXMiner output 4 files for each analysis - allScores, spectraScores, peptideScores, and crosslinkSiteScores. "allScores" contain all raw 
spectrum-sequence match data which can be reprocessed again through ZXMiner (by selecting the "allScores" file for "output file name"). "spectraScores" is
obtained from "allScores" by determining the best peptide for each spectrum using GM scores. "peptideScores" is obtained from "spectraScores" by determining the best
spectrum for each peptide identified. Finally, "crosslinkSiteScores" contain detailed score information for each peptides passing the defined false discovery rate (FDR)
and include individual scores for every possible cross-link sites for peptides with multiple cross-linkable sites (accept those that were filtered out from "allScores" 
by the "minScoreSite" setting in configuration file).

	2.5 reprocessing of allscores.out
		This option is provided so that the FDR can be re-calculated under different score filters without the need to re-run the entire process in ZXMiner.
This mode is activated by selecting the "allScores" file for "output file name". Right now, there are two filters available for reprocessing: ion coverage and deltaGM.
We will include more control in the future versions of ZXMiner.

2.6 running multiple instance of ZXMiner
More CPU cores can be allocated to ZXMiner to speed up the process. However, if situations come up where running two ZXMiner processes in parallel is more desirable,
because ZXMiner writes files with fixed names to the "temp" folder dynamically, multiple processes of ZXMiner will interfere with each other. To bypass this limitation, 
please copy ZXMiner into two locations so that they can be run independently.

3. SOURCE CODE & CITATION
To allow for wider application of ZXMiner, the source code is provided under the MIT license so that others 
can optimize the MS/MS processing step, peptide scoring step, etc. to suit their circumstances.

To cite ZXMiner, please use:
Sriswasdi S, Harper SL, Tang H-Y, Speicher DW. Enhanced identification of zero-length chemical 
cross-links using label-free quantitation and high-resolution fragment ion spectra. J Proteome Res 
13(2):898-914 (2014).

and

Rivera-Santiago, R.F., Sriswasdi, S., Harper, S.L., and Speicher, D.W. Probing structures of large 
protein complexes using zero-length cross-linking. Methods. 89:99-111 (2015).
