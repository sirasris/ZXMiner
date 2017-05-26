import java.util.TreeMap;
import java.util.concurrent.*;
import java.io.RandomAccessFile;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// multi-threading for extracting indices from mzXML file
public class IndexGrabberThread implements Callable<TreeMap<Integer, Long>>
{
	int mzXMLFileID;
	String mzXMLPath;
	Pattern indexPattern;

	public IndexGrabberThread(int mzXMLFileID, String mzXMLPath, Pattern indexPattern)
	{
		this.mzXMLFileID = mzXMLFileID;
		this.mzXMLPath = mzXMLPath;
		this.indexPattern = indexPattern;
	}

	public TreeMap<Integer, Long> call()
	{
		try
		{
			TreeMap<Integer, Long> tempmap = new TreeMap<Integer, Long>();
			String st;
			Matcher match;

			RandomAccessFile reader = new RandomAccessFile(mzXMLPath + ".mzXML", "r");
			reader.seek((long) (reader.length() * 0.999)); // heuristic seek

			st = reader.readLine().trim();

			while (!st.equals("<index name=\"scan\">")) // seek for index info
			{
				// System.out.println(st);
				st = reader.readLine().trim();

				if (st.startsWith("<indexOffset>")) // heuristic seek is too deep
				{
					long offset = Long.valueOf(st.split(">")[1].split("<")[0]);
					reader.seek(offset);
					st = reader.readLine().trim(); // this should be '<index name="scan">'
					break;
				}
			}

			st = reader.readLine().trim();

			while (!st.equals("</index>"))
			{
				// System.out.println(st);
				match = indexPattern.matcher(st);
				st = reader.readLine().trim();

				if (match.find())
				{
					tempmap.put(new Integer(match.group(1)), new Long(match.group(2)));
					// System.out.println(match.group(1));
					// System.out.println(match.group(2));
				}

				else
					HelperFunctions.debug("Unexpected pattern in mzXML index!");
			}

			reader.close();
			return tempmap;
		}

		catch (Exception e)
		{ HelperFunctions.debug("IndexGrabberThread", HelperFunctions.getStackTrace(e)); }

		return null; // fails
	}
}