import java.io.FileWriter;
import java.util.ArrayList;

// do in silico digestion and precursor matching at the same time
public class ProcessAndMatchPrecursor implements Runnable
{
	SynchronizedTreeMap globalMatches; // global precursor match list
	SynchronizedTreeMap[] globalResult;
	ArrayList<Double> peptideResult; // peptide mass lists of this thread name
	ArrayList<Double> crosslinkResult; // crosslink mass lists of this thread name
	ArrayList<ProteinStruct> proteins;
	ThreadLimitStruct scope;
	ParamStruct param;
	double massTolerance; // mass tolerance
	Double[] observedMasses; // list of observed masses, sorted
	Integer[] indexedObservedMasses; // map to original precursor ID
	int threadID;
	
	public ProcessAndMatchPrecursor(SynchronizedTreeMap[] globalResult, ArrayList<ProteinStruct> proteins, ThreadLimitStruct scope, ParamStruct param, SynchronizedTreeMap globalMatches, Double[] observedMasses, Integer[] indexedObservedMasses)
	{
		this.globalResult = globalResult;
   		this.proteins = proteins;
   		this.scope = scope;
		this.param = param;
		this.globalMatches = globalMatches;
		this.observedMasses = observedMasses;
		this.indexedObservedMasses = indexedObservedMasses;
   		massTolerance = param.precursorTolerance;
	}
	
    // check whether two peptides are adjacent
    public boolean isAdjacent(PeptideStruct peptideA, PeptideStruct peptideB)
    {
    	if (peptideA.parent.entryID == peptideB.parent.entryID) // same parent
    	{
    		if (peptideA.startPos == peptideB.endPos + 1) // A after B
    			return true;
    		if (peptideB.startPos == peptideA.endPos + 1) // B after A
    			return true;
    	}
    	
    	return false;
    }
    
   	// iterate and generate all modified form of a peptide
   	// consider residue at 'currentPos'
   	public void generateAllModForms(int currentPos, PeptideStruct original, ArrayList<PeptideStruct> list)
	{
		if (original.getNumVarMod() == param.maxVarModPerPeptide) // maximum number of variable modifications reached
		{
			list.add(original);
			return;
		}

		if (currentPos == original.getSequenceLength(true)) // we complete one construction
			list.add(original);
		else
		{
			if (MassInfo.hasFixedMod(original.getResidue(currentPos))) // skip residues with fixed modification
				generateAllModForms(currentPos + 1, original, list);

			else
			{
				ArrayList<ModificationStruct> modList = MassInfo.getAllVarMod(original.getResidue(currentPos));

				if (modList != null) // this residue can be modified
				{
					PeptideStruct newPeptide;

					for (int i = 0; i < modList.size(); i++)
					{
						newPeptide = original.clone(); // create a deep copy
						newPeptide.addModification(currentPos, modList.get(i)); // add new modification
						generateAllModForms(currentPos + 1, newPeptide, list); // continue the recursion
					}
				}

				generateAllModForms(currentPos + 1, original, list); // also always proceed without adding anything
			}
		}
	}

   	
   	public void run()
    {
   		ProteinStruct protein1 = null, protein2 = null;
   		String threadName = Thread.currentThread().getName();
   		threadID = Integer.valueOf(threadName.split("-")[1]);
   		// System.out.print(threadName + ", ");
   		// HelperFunctions.debug("limit", scope);
   		
        try
        {
        	FileWriter outfilePeptide, outfileCrosslink;
        	peptideResult = new ArrayList<Double>();
        	crosslinkResult = new ArrayList<Double>();

			if (globalResult[0].get(threadName).size() == 0) // create new file for the first time
				outfilePeptide = new FileWriter(".\\temp\\peptide-" + threadName + ".temp"); // outfilePeptide = new BufferedWriter(new FileWriter(".\\temp\\peptide-" + threadName + ".temp"));
			else // append otherwise
				outfilePeptide = new FileWriter(".\\temp\\peptide-" + threadName + ".temp", true); // outfilePeptide = new BufferedWriter(new FileWriter(".\\temp\\peptide-" + threadName + ".temp", true));

			if (globalResult[1].get(threadName).size() == 0) // create new file for the first time
				outfileCrosslink = new FileWriter(".\\temp\\crosslink-" + threadName + ".temp"); // outfileCrosslink = new BufferedWriter(new FileWriter(".\\temp\\crosslink-" + threadName + ".temp"));
			else // append otherwise
				outfileCrosslink = new FileWriter(".\\temp\\crosslink-" + threadName + ".temp", true); // outfileCrosslink = new BufferedWriter(new FileWriter(".\\temp\\crosslink-" + threadName + ".temp", true));

		// CODE FOR DEBUGGING THREADING
			// PeptideStruct peptide = new PeptideStruct(proteins.get(0), 1, 5);
			// localResult.add(new Double(localResult.size()));
			// outfile.write(HelperFunctions.encodeBase64(peptide, param.maxVarModPerPeptide));

			protein1 = proteins.get(scope.sourceID1);
			ArrayList<Integer> cleaveSites1 = protein1.getCleaveSites();
			PeptideStruct temppeptide1;
			ArrayList<PeptideStruct> templist1;
			// String combinedOutput = "";
			int tempint;
			double tempmass;

			if (scope.sourceID2 == -1) // linear peptide, dead-end, or loop
			{
				for (int i = scope.startID; i < scope.endID + 1; i++) // for each starting position
				for (int j = 0; j < param.maxMissedCleave; j++) // for each level of missed cleavage
				{
					if (i + j + 1 < cleaveSites1.size()) // make sure to not exceed the end of sequence
					{
						// System.out.println(threadName + ":" + i + "," + j);
						tempint = protein1.getSequenceLength(false, cleaveSites1.get(i) + 1, cleaveSites1.get(i + j + 1));

						if (tempint >= param.minPeptideLength && tempint <= param.maxPeptideLength) // check sequence length
						{
							temppeptide1 = new PeptideStruct(protein1, cleaveSites1.get(i) + 1, cleaveSites1.get(i + j + 1));
							// System.out.println(threadName + ":" + temppeptidHelperFunctions.getStackTrace(e));

							templist1 = new ArrayList<PeptideStruct>();
							generateAllModForms(0, temppeptide1, templist1); // populate variable modifications

							for (int k = 0; k < templist1.size(); k++) // for each form of this peptide
							{
								// peptideResult.add(new Double(MassInfo.getMass(templist1.get(k)))); // record mass
								// outfilePeptide.write(Base64Parser.encodeBase64(templist1.get(k), param.maxVarModPerPeptide)); // output peptide to file
								// globalResult[0].addAndWrite(threadName, new Double(MassInfo.getMass(templist1.get(k))), outfilePeptide, Base64Parser.encodeBase64(templist1.get(k), param.maxVarModPerPeptide));
								tempmass = MassInfo.getMass(templist1.get(k));
								globalResult[0].checkMassAddAndWrite(threadName, new Double(tempmass), outfilePeptide, Base64Parser.encodeBase64(templist1.get(k), param.maxVarModPerPeptide), observedMasses, indexedObservedMasses, tempmass, massTolerance, globalMatches, threadID, false);
								
								if (param.hasDeadEnd()) // generate loop // dead-ends were taken care of as variable modifications
								{
								}
							}
						}
					}
				}
				
				// globalResult[0].addToValueList(threadName, peptideResult);
				// outfilePeptide.flush();
				outfilePeptide.close();
				outfileCrosslink.close();
				return;
			}

			PeptideStruct temppeptide2;
			boolean[] crosslinkSites;
			ArrayList<PeptideStruct> templist2;
			CrosslinkStruct crosslink;

			if (scope.sourceID2 == scope.sourceID1) // same-protein crosslink
			{
				for (int i1 = scope.startID; i1 < scope.endID + 1; i1++) // for each starting position
				for (int j1 = 0; j1 < param.maxMissedCleave; j1++) // for each level of missed cleavage
				{
					if (i1 + j1 + 1 < cleaveSites1.size()) // make sure to not exceed the end of sequence
					{
						// System.out.println(threadName + ":" + i1 + "," + j1);
						tempint = protein1.getSequenceLength(false, cleaveSites1.get(i1) + 1, cleaveSites1.get(i1 + j1 + 1));

						if (tempint >= param.minPeptideLength && tempint <= param.maxPeptideLength) // check sequence length
						{
							temppeptide1 = new PeptideStruct(protein1, cleaveSites1.get(i1) + 1, cleaveSites1.get(i1 + j1 + 1));
							// System.out.println(threadName + ":" + temppeptidHelperFunctions.getStackTrace(e));

							templist1 = new ArrayList<PeptideStruct>();
							generateAllModForms(0, temppeptide1, templist1); // populate variable modifications

							for (int k1 = 0; k1 < templist1.size(); k1++) // for each form of this peptide
							{
								crosslinkSites = CrosslinkSiteIdentifier.containsCrosslinkSites(templist1.get(k1), param.crosslinker); // preprocess crosslink site
								
								for (int j2 = j1 + 1; j2 < param.maxMissedCleave; j2++) // same start site, different missed cleavage
								{
									if (i1 + j2 + 1 < cleaveSites1.size()) // make sure to not exceed the end of sequence
									{
										tempint = protein1.getSequenceLength(false, cleaveSites1.get(i1) + 1, cleaveSites1.get(i1 + j2 + 1));

										if (tempint >= param.minPeptideLength && tempint <= param.maxPeptideLength) // check sequence length
										{
											temppeptide2 = new PeptideStruct(protein1, cleaveSites1.get(i1) + 1, cleaveSites1.get(i1 + j2 + 1));
											// System.out.println(threadName + ":" + temppeptide2.toString());

											if (param.allowAdjacent || !isAdjacent(temppeptide1, temppeptide2)) // check adjacent peptides
											{
												if (CrosslinkSiteIdentifier.isCrosslinkable(crosslinkSites, temppeptide2, param.crosslinker)) // original peptide2 crosslinkable
												{
													templist2 = new ArrayList<PeptideStruct>();
													generateAllModForms(0, temppeptide2, templist2); // populate variable modifications

													for (int k2 = 0; k2 < templist2.size(); k2++) // for each form of this peptide
													{
														if (CrosslinkSiteIdentifier.isCrosslinkable(crosslinkSites, templist2.get(k2), param.crosslinker)) // final peptide crosslinkable
														{
															crosslink = new CrosslinkStruct(templist1.get(k1), templist2.get(k2)); // create crosslink
															// crosslinkResult.add(new Double(MassInfo.getMass(crosslink, param))); // record mass
															// outfileCrosslink.write(Base64Parser.encodeBase64(crosslink, param.maxVarModPerPeptide)); // output peptide to file
															// globalResult[1].addAndWrite(threadName, new Double(MassInfo.getMass(crosslink, param)), outfileCrosslink, Base64Parser.encodeBase64(crosslink, param.maxVarModPerPeptide));
															tempmass = MassInfo.getMass(crosslink, param);
															globalResult[1].checkMassAddAndWrite(threadName, new Double(tempmass), outfileCrosslink, Base64Parser.encodeBase64(crosslink, param.maxVarModPerPeptide), observedMasses, indexedObservedMasses, tempmass, massTolerance, globalMatches, threadID, true);
														}
													}
												}
											}
										}
									}
								}

								for (int i2 = i1 + 1; i2 < cleaveSites1.size() - 1; i2++) // for peptide starting after current startID
								for (int j2 = 0; j2 < param.maxMissedCleave; j2++) // for each level of missed cleavage
								{
									if (i2 + j2 + 1 < cleaveSites1.size()) // make sure to not exceed the end of sequence
									{
										tempint = protein1.getSequenceLength(false, cleaveSites1.get(i2) + 1, cleaveSites1.get(i2 + j2 + 1));

										if (tempint >= param.minPeptideLength && tempint <= param.maxPeptideLength) // check sequence length
										{
											temppeptide2 = new PeptideStruct(protein1, cleaveSites1.get(i2) + 1, cleaveSites1.get(i2 + j2 + 1));
											// System.out.println(threadName + ":" + temppeptide2.toString());
											
											if (param.allowAdjacent || !isAdjacent(temppeptide1, temppeptide2)) // check adjacent peptides
											{
												if (CrosslinkSiteIdentifier.isCrosslinkable(crosslinkSites, temppeptide2, param.crosslinker)) // original peptide2 crosslinkable
												{
													templist2 = new ArrayList<PeptideStruct>();
													generateAllModForms(0, temppeptide2, templist2); // populate variable modifications
	
													for (int k2 = 0; k2 < templist2.size(); k2++) // for each form of this peptide
													{
														if (CrosslinkSiteIdentifier.isCrosslinkable(crosslinkSites, templist2.get(k2), param.crosslinker)) // final peptide crosslinkable
														{
															crosslink = new CrosslinkStruct(templist1.get(k1), templist2.get(k2)); // create crosslink
															// crosslinkResult.add(new Double(MassInfo.getMass(crosslink, param))); // record mass
															// outfileCrosslink.write(Base64Parser.encodeBase64(crosslink, param.maxVarModPerPeptide)); // output peptide to file
															// globalResult[1].addAndWrite(threadName, new Double(MassInfo.getMass(crosslink, param)), outfileCrosslink, Base64Parser.encodeBase64(crosslink, param.maxVarModPerPeptide));
															tempmass = MassInfo.getMass(crosslink, param);
															globalResult[1].checkMassAddAndWrite(threadName, new Double(tempmass), outfileCrosslink, Base64Parser.encodeBase64(crosslink, param.maxVarModPerPeptide), observedMasses, indexedObservedMasses, tempmass, massTolerance, globalMatches, threadID, true);

														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
				
				// globalResult[1].addToValueList(threadName, crosslinkResult);
				// outfileCrosslink.flush();
				outfilePeptide.close();
				outfileCrosslink.close();
				return;
			}

			// inter-protein crosslink
			protein2 = proteins.get(scope.sourceID2);
			ArrayList<Integer> cleaveSites2 = protein2.getCleaveSites();

			for (int i1 = scope.startID; i1 < scope.endID + 1; i1++) // for each starting position
			for (int j1 = 0; j1 < param.maxMissedCleave; j1++) // for each level of missed cleavage
			{
				if (i1 + j1 + 1 < cleaveSites1.size()) // make sure to not exceed the end of sequence
				{
					// System.out.println(threadName + ":" + i1 + "," + j1);
					tempint = protein1.getSequenceLength(false, cleaveSites1.get(i1) + 1, cleaveSites1.get(i1 + j1 + 1));

					if (tempint >= param.minPeptideLength && tempint <= param.maxPeptideLength) // check sequence length
					{
						temppeptide1 = new PeptideStruct(protein1, cleaveSites1.get(i1) + 1, cleaveSites1.get(i1 + j1 + 1));
						// System.out.println(threadName + ":" + temppeptidHelperFunctions.getStackTrace(e));

						templist1 = new ArrayList<PeptideStruct>();
						generateAllModForms(0, temppeptide1, templist1); // populate variable modifications

						for (int k1 = 0; k1 < templist1.size(); k1++) // for each form of this peptide
						{
							crosslinkSites = CrosslinkSiteIdentifier.containsCrosslinkSites(templist1.get(k1), param.crosslinker); // preprocess crosslink site

							for (int i2 = 0; i2 < cleaveSites2.size() - 1; i2++) // for each starting position on second protein
							for (int j2 = 0; j2 < param.maxMissedCleave; j2++) // for each level of missed cleavage
							{
								if (i2 + j2 + 1 < cleaveSites2.size()) // make sure to not exceed the end of sequence
								{
									tempint = protein2.getSequenceLength(false, cleaveSites2.get(i2) + 1, cleaveSites2.get(i2 + j2 + 1));

									if (tempint >= param.minPeptideLength && tempint <= param.maxPeptideLength) // check sequence length
									{
										temppeptide2 = new PeptideStruct(protein2, cleaveSites2.get(i2) + 1, cleaveSites2.get(i2 + j2 + 1));
										// System.out.println(threadName + ":" + temppeptide2.toString());

										if (CrosslinkSiteIdentifier.isCrosslinkable(crosslinkSites, temppeptide2, param.crosslinker)) // original peptide2 crosslinkable
										{
											templist2 = new ArrayList<PeptideStruct>();
											generateAllModForms(0, temppeptide2, templist2); // populate variable modifications

											for (int k2 = 0; k2 < templist2.size(); k2++) // for each form of this peptide
											{
												if (CrosslinkSiteIdentifier.isCrosslinkable(crosslinkSites, templist2.get(k2), param.crosslinker)) // final peptide crosslinkable
												{
													crosslink = new CrosslinkStruct(templist1.get(k1), templist2.get(k2)); // create crosslink
													// crosslinkResult.add(new Double(MassInfo.getMass(crosslink, param))); // record mass
													// outfileCrosslink.write(Base64Parser.encodeBase64(crosslink, param.maxVarModPerPeptide)); // output peptide to file
													// globalResult[1].addAndWrite(threadName, new Double(MassInfo.getMass(crosslink, param)), outfileCrosslink, Base64Parser.encodeBase64(crosslink, param.maxVarModPerPeptide));
													tempmass = MassInfo.getMass(crosslink, param);
													globalResult[1].checkMassAddAndWrite(threadName, new Double(tempmass), outfileCrosslink, Base64Parser.encodeBase64(crosslink, param.maxVarModPerPeptide), observedMasses, indexedObservedMasses, tempmass, massTolerance, globalMatches, threadID, true);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			
			// globalResult[1].addToValueList(threadName, crosslinkResult);
			// outfileCrosslink.flush();
			outfilePeptide.close();
			outfileCrosslink.close();
        }

        catch (Exception e)
        {
        	HelperFunctions.debug("ProteinProcessorThread", HelperFunctions.getStackTrace(e));
        	HelperFunctions.debug("protein1", protein1);
        	HelperFunctions.debug("protein1", protein2);
        }
    }

}
