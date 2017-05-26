// data structure for thread-specific peptide data
public class ThreadPeptideStruct
{
	public final int threadID;
	public final boolean isCrosslink;
	public final int peptideID;

	public ThreadPeptideStruct(int threadID, boolean isCrosslink, int peptideID)
	{
		this.threadID = threadID;
		this.isCrosslink = isCrosslink;
		this.peptideID = peptideID;
	}
}