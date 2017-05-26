// generalized data structure define the scope of each thread
public class ThreadLimitStruct
{
	public final int sourceID1, sourceID2; // entry ID of source data (eg. 'raw file' or 'protein')
										   // sourceID2 is used for crosslinks and can be '-1' in other cases
	public final int startID, endID; // range of entry IDs of subdata (eg. 'scan number' or 'cleaveage site')

	public ThreadLimitStruct(int sourceID1, int sourceID2, int startID, int endID)
	{
		this.sourceID1 = sourceID1;
		this.sourceID2 = sourceID2;
		this.startID = startID;
		this.endID = endID;
	}

	public ThreadLimitStruct(int sourceID1, int startID, int endID)
	{
		this.sourceID1 = sourceID1;
		this.sourceID2 = -1;
		this.startID = startID;
		this.endID = endID;
	}

	public String toString()
	{
		String details = "sourceID1: " + sourceID1;
		details += ", sourceID2: " + sourceID2;
		details += ", startID: " + startID;
		details += ", endID: " + endID;

		return details;
	}
}