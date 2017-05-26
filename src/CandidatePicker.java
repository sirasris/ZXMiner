import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeMap;

// pick the best explanation for each scan
public class CandidatePicker
{
	public CandidatePicker()
	{}
	
	// select best peptide for each spectrum
	public static void pickSpectra(ParamStruct param)
	{
		try
		{
			TreeMap<String, String> candidateInfo = new TreeMap<String, String>();
			TreeMap<String, Double> bestGMScore = new TreeMap<String, Double>();
			TreeMap<String, Double> bestDCNScore = new TreeMap<String, Double>();
			String[] tempst;
			String key, header, st;
			double dCn;
			
			BufferedReader infile = new BufferedReader(new FileReader(".\\output\\" + param.outputFileName + "_allScores.out"));
			header = infile.readLine(); // header
			st = infile.readLine();
			
			while (st != null)
			{
				tempst = st.split("\t");
				
				if (tempst[17].equals("n/a") || tempst[18].equals("n/a")) // filter for generic crosslink or linear peptide
				{
					key = tempst[0] + "::" + tempst[1];
					
					if (!candidateInfo.containsKey(key)) // new entry
					{
						candidateInfo.put(key, st);
						bestGMScore.put(key, new Double(tempst[tempst.length - 1]));
						bestDCNScore.put(key, new Double(tempst[tempst.length - 1]));
					}
					
					else // check GM score
					{
						dCn = Double.valueOf(tempst[tempst.length - 1]) - bestGMScore.get(key).doubleValue();
						
						if (dCn >= 0) // update
						{
							candidateInfo.put(key, st);
							bestDCNScore.put(key, new Double(dCn));
							bestGMScore.put(key, new Double(tempst[tempst.length - 1]));
						}
						
						else // still check if dCn should be update
						{
							if (-dCn < bestDCNScore.get(key).doubleValue())
							{
								bestDCNScore.put(key, new Double(-dCn));
							}
						}
					}
				}
				
				st = infile.readLine();
			}
			
			infile.close();
			
			// now, output to another file
			FileWriter outfile = new FileWriter(".\\output\\" + param.outputFileName + "_spectraScores.out");
			outfile.write(header + "\tdeltaGM\n");
			
			for (Iterator<String> iter = candidateInfo.keySet().iterator(); iter.hasNext();)
			{
				key = iter.next();
				outfile.write(candidateInfo.get(key) + "\t" + bestDCNScore.get(key) + "\n");
			}
			
			outfile.close();
		}
		
		catch (Exception e)
		{ HelperFunctions.debug("CandidatePicker::pickSpectra", HelperFunctions.getStackTrace(e)); }
	}
	
	// select best spectrum for each crosslinked peptide
	public static void pickPeptide(ParamStruct param)
	{
		try
		{
			TreeMap<String, String> candidateInfo = new TreeMap<String, String>();
			TreeMap<String, Double> bestGMScore = new TreeMap<String, Double>();
			TreeMap<String, String> entryType = new TreeMap<String, String>();
			String[] tempst;
			String key, header, st;
			
			BufferedReader infile = new BufferedReader(new FileReader(".\\output\\" + param.outputFileName + "_spectraScores.out"));
			header = infile.readLine(); // header
			st = infile.readLine();
			
			while (st != null)
			{
				tempst = st.split("\t");
				
				if (!tempst[18].equals("n/a") && Double.parseDouble(tempst[31]) > param.minIonCov 
						&& Double.parseDouble(tempst[34]) > param.minIonCov && Double.parseDouble(tempst[37]) > param.minDeltaScore) // filter for crosslink only and apply score filters
				{
					key = tempst[2] + "::" + tempst[10] + "::" + tempst[13] + "::" + tempst[14] + "::" + tempst[16];
					key = key + "::" + tempst[18] + "::" + tempst[21] + "::" + tempst[22] + "::" + tempst[24];
					
					if (!candidateInfo.containsKey(key)) // new entry
					{
						entryType.put(key, tempst[12] + "::" + tempst[20]);
						candidateInfo.put(key, st);
						bestGMScore.put(key, Double.valueOf(tempst[36]));
					}
					
					else // check GM score
					{
						if (Double.parseDouble(tempst[36]) > bestGMScore.get(key).doubleValue())
						{						
							candidateInfo.put(key, st);
							bestGMScore.put(key, Double.valueOf(tempst[36]));
						}
					}
				}
				
				st = infile.readLine();
			}
			
			infile.close();
			
			// compute FDR
			TreeMap<String, Double> FDR = computeFDR(bestGMScore, entryType);
			
			// now, output to another file
			FileWriter outfile = new FileWriter(".\\output\\" + param.outputFileName + "_peptideScores.out");
			outfile.write(header + "\tClass\tFDR\n");
			
			for (Iterator<String> iter = candidateInfo.keySet().iterator(); iter.hasNext();)
			{
				key = iter.next();
				outfile.write(candidateInfo.get(key));
				
				switch (entryType.get(key))
				{
					case "true::true": outfile.write("\tForward\t"); break;
					case "true::false": outfile.write("\tHybrid\t"); break;
					case "false::true": outfile.write("\tHybrid\t"); break;
					case "false::false": outfile.write("\tReversed\t"); break;
					default: outfile.write("\tn/a\n"); break;
				}
				
				outfile.write(FDR.get(key).toString() + "\n");
			}
			
			outfile.close();
		}
		
		catch (Exception e)
		{ HelperFunctions.debug("CandidatePicker::pickPeptide", HelperFunctions.getStackTrace(e)); }
	}
	
	// report all crosslink sites for high-confidence crosslinked peptides
	public static void pickHighConfPeptide(ParamStruct param)
	{
		try
		{
			TreeMap<String, Integer> highconfCrosslinks = new TreeMap<String, Integer>();
			String[] tempst;
			String key, header, st;
			int count = 1;
			
			BufferedReader infile = new BufferedReader(new FileReader(".\\output\\" + param.outputFileName + "_peptideScores.out"));
			header = infile.readLine(); // header	
			header = header.split("deltaGM")[0];
			st = infile.readLine();
			
			while (st != null)
			{
				tempst = st.split("\t");
				
				if (tempst[tempst.length - 2].equals("Forward") && Double.valueOf(tempst[tempst.length - 1]) <= 100 * param.maximumFDR) // filter for low FDR and forward sequence
				{
					key = tempst[0] + "::" + tempst[1];
					key = key + tempst[2] + "::" + tempst[10] + "::" + tempst[13] + "::" + tempst[14] + "::" + tempst[16];
					key = key + "::" + tempst[18] + "::" + tempst[21] + "::" + tempst[22] + "::" + tempst[24];
					
					highconfCrosslinks.put(key, new Integer(count));
					count++;
				}
				
				st = infile.readLine();
			}
			
			infile.close();
			
			// System.out.println(highconfCrosslinks);
			
			// extract all crosslinked site information
			infile = new BufferedReader(new FileReader(".\\output\\" + param.outputFileName + "_allScores.out"));
			st = infile.readLine();
			st = infile.readLine();
			
			TreeMap<Integer, ArrayList<String>> allSites = new TreeMap<Integer, ArrayList<String>>();
			ArrayList<String> templist;
			Integer keyID;
			
			while (st != null)
			{
				tempst = st.split("\t");
				
				key = tempst[0] + "::" + tempst[1];
				key = key + tempst[2] + "::" + tempst[10] + "::" + tempst[13] + "::" + tempst[14] + "::" + tempst[16];
				key = key + "::" + tempst[18] + "::" + tempst[21] + "::" + tempst[22] + "::" + tempst[24];
				
				if (tempst[12].equals("true") && tempst[20].equals("true") && highconfCrosslinks.containsKey(key))
				{			
					keyID = highconfCrosslinks.get(key);
					// System.out.println(highconfCrosslinks.containsKey(key));
					// System.out.println(keyID);
					
					if (allSites.containsKey(keyID)) // add to an existing entry
					{
						templist = allSites.get(keyID);
						templist.add(st);
					}
					
					else
					{
						templist = new ArrayList<String>();
						templist.add(st);
						allSites.put(keyID, templist);
					}
				}
				
				st = infile.readLine();
			}
			
			infile.close();
			
			// System.out.println(allSites);
			
			// now, output to another file
			FileWriter outfile = new FileWriter(".\\output\\" + param.outputFileName + "_crosslinkSiteScores.out");
			outfile.write("CrosslinkID\t" + header + "\n");
			
			for (Iterator<String> iter = highconfCrosslinks.keySet().iterator(); iter.hasNext();)
			{
				keyID = highconfCrosslinks.get(iter.next());
				templist = allSites.get(keyID);
				
				for (Iterator<String> iter2 = templist.iterator(); iter2.hasNext();)
					outfile.write(keyID.intValue() + "\t" + iter2.next() + "\n");
			}
			
			outfile.close();
		}
		
		catch (Exception e)
		{ HelperFunctions.debug("CandidatePicker::pickHighConfPeptide", HelperFunctions.getStackTrace(e)); }
	}
	
	public static TreeMap<String, Double> computeFDR(TreeMap<String, Double> scores, TreeMap<String, String> types)
	{	
		String[] keys = scores.keySet().toArray(new String[0]); // extract data into array
		Double[] score = new Double[scores.size()];
		
		for (int i = 0; i < scores.size(); i++)
			score[i] = scores.get(keys[i]);
		
		String[] type = new String[types.size()];
		
		for (int i = 0; i < types.size(); i++)
			type[i] = types.get(keys[i]);
		
		Integer[] index = HelperFunctions.getIndexArray(score); // sort score and obtain indexes
		
		String[] sortedTypes = new String[types.size()];
		
		for (int i = 0; i < types.size(); i++) // sort types
			sortedTypes[i] = type[index[i].intValue()];
		
		String[] sortedPeptides = new String[scores.size()];
		
		for (int i = 0; i < scores.size(); i++) // sort peptides
			sortedPeptides[i] = keys[index[i].intValue()];
		
		int[] count = new int[3]; // count the number of forward, hybrid, and reversed entries
		TreeMap<String, Double> FDR = new TreeMap<String, Double>();
		
		for (int i = scores.size() - 1; i > -1; i--) // compute FDR
		{
			// System.out.println(sortedTypes[i] + "\t" + score[i] + "\t" + sortedPeptides[i]);
			
			switch (sortedTypes[i])
			{
				case "true::true": count[0]++; break;
				case "true::false": count[1]++; break;
				case "false::true": count[1]++; break;
				case "false::false": count[2]++; break;
				default: HelperFunctions.debug("CandidatePicker::computeFDR", "Error! Unknown class of crosslink!"); break;
			}
			
			FDR.put(sortedPeptides[i], new Double((count[1] - count[2]) * 100.0 / (scores.size() - i)));
		}
		
		return FDR;
	}
}
