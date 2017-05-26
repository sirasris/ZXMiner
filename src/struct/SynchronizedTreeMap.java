import java.util.TreeMap;
import java.util.ArrayList;
import java.util.Set;
import java.io.FileWriter;

// tree map version for multi-threading
// every value must be ArrayList
public class SynchronizedTreeMap
{
	private TreeMap data;

	public SynchronizedTreeMap()
	{ data = new TreeMap(); }

	public synchronized void put(Object key, ArrayList value)
	{ data.put(key, value); }

	public synchronized int size(Object key)
	{
		if (data.containsKey(key))
			return ((ArrayList) data.get(key)).size();
		else
			return -1;
	}

	public synchronized Set keySet()
	{ return data.keySet(); }

	public synchronized ArrayList get(Object key)
	{
		if (data.containsKey(key))
			return ((ArrayList) data.get(key));
		else
			return null;
	}

	public synchronized void addToValueList(Object key, Object value)
	{
		if (data.containsKey(key))
			((ArrayList) data.get(key)).add(value);
		else
		{
			ArrayList temp = new ArrayList();
			temp.add(value);
			data.put(key, temp);
		}
	}
	
	// add to list and write to output file at the same time to maintain perfect order
	public synchronized void addAndWrite(Object key, Object value, FileWriter writer, String line)
	{
		if (data.containsKey(key))
			((ArrayList) data.get(key)).add(value);
		else
		{
			ArrayList temp = new ArrayList();
			temp.add(value);
			data.put(key, temp);
		}
		
		try
		{
			writer.write(line);
		}
		
		catch (Exception e)
		{ HelperFunctions.debug("SynchronizedTreeMap", HelperFunctions.getStackTrace(e)); }
	}

	public synchronized void addToValueList(Object key, ArrayList values)
	{
		if (data.containsKey(key))
			((ArrayList) data.get(key)).addAll(values);
		else
			data.put(key, values);
	}

	public synchronized String toString()
	{ return data.toString(); }
	
	// check precursor mass first
	public synchronized void checkMassAddAndWrite(Object key, Object value, FileWriter writer, String line, Double[] observedMasses, Integer[] indexedObservedMasses, double peptideMass, double massTolerance, SynchronizedTreeMap globalMatches, int threadID, boolean isCrosslink)
	{
		ArrayList<Integer> tempmatches = MassMatcher.match(observedMasses, peptideMass, massTolerance);
		ThreadPeptideStruct tempThreadPeptide;

		if (tempmatches.size() > 0) // found some matches
		{
			if (data.containsKey(key)) // add value to map
			{
				((ArrayList) data.get(key)).add(value);
				tempThreadPeptide = new ThreadPeptideStruct(threadID, isCrosslink, ((ArrayList) data.get(key)).size() - 1);
			}
			
			else
			{
				ArrayList tempList = new ArrayList();
				tempList.add(value);
				data.put(key, tempList);
				tempThreadPeptide = new ThreadPeptideStruct(threadID, isCrosslink, 0);
			}
			
			for (int j = 0; j < tempmatches.size(); j++) // update global matches
				globalMatches.addToValueList(indexedObservedMasses[tempmatches.get(j).intValue()], tempThreadPeptide);			
			
			try // only write mached peptides
			{
				writer.write(line);
			}
			
			catch (Exception e)
			{ HelperFunctions.debug("SynchronizedTreeMap::checkMassAddAndWrite", HelperFunctions.getStackTrace(e)); }
		}
	}
	
	// clear resource
	public void clear()
	{ data.clear(); }
}